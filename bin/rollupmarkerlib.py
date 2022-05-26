# Purpose: to encapsulate the rollup rules for rolling up genotype-level
#	MP and disease annotations to markers, eliminating paths where we
#	cannot determine a single causative marker.
# Notes: We will want to carefully consider memory management in this module,
#	as we are dealing with over 270,000 annotations for more than 60,000
#	genotypes as of August 2014.  Add in evidence records and notes, and
#	memory needs could be substantial.  Also, we'll really want to do most
#	batch processing in SQL, rather than iterating in code (where possible)
#	to reduce the number of repetitive SQL commands that would be needed.

import gc
import os
import db
import copy
import mgi_utils

###--- globals ---###

db.setTrace()
Error = 'rollupmarkerlib.Error'
DEBUG = False

MAX_ANNOTATIONS = 5000		# maximum number of annotations to cache in

CURRENT_ANNOT_TYPE = None	# either DO_GENOTYPE or MP_GENOTYPE
DO_GENOTYPE = 1020		# DO/Genotype annotation type
MP_GENOTYPE = 1002		# Mammalian Phenotype/Genotype annotation type

NO_PHENOTYPIC_ANALYSIS = 293594	# term key for 'no phenotypic analysis' term
SOURCE_ANNOT_KEY = None		# term key for _SourceAnnot_key property

GT_ROSA = 37270			# marker key for Gt(ROSA)26Sor marker
HPRT = 9936			# marker key for Hprt marker
COL1A1 = 1092			# marker key for Col1a1 marker

DOCKING_SITES = [ GT_ROSA, HPRT, COL1A1 ]	# loci that can be generally knocked into without causing a phenotype

INITIALIZED = False		# have we finished initializing this module?

ANNOTATION_COUNTS = {}		# marker key -> count of rolled-up annotations

MARKER_KEYS = []		# ordered list of marker keys

LAST_MARKER_KEY_INDEX = None	# index into MARKER_KEYS of last marker key which had its details loaded

MARKERS_TO_DO = []		# list of markers loaded and waiting to be processed

TERM_MAP = None			# KeyMap for term key -> term ID
MARKER_MAP = None		# KeyMap for marker key -> marker ID
JNUM_MAP = None			# KeyMap for refs key -> J: number
EVIDENCE_MAP = None		# KeyMap for evidence key -> evidence abbrev.
QUALIFIER_MAP = None		# KeyMap for qualifier key -> qualifier term
USER_MAP = None			# KeyMap for user key -> user
PROPERTY_MAP = None		# KeyMap for property key -> property name

# no testing
testSQL = ""

# rule #1/crm159
#testSQL = '''
#and exists (select 1 from ACC_Accession testg 
#where gag._genotype_key = testg._object_key and testg._mgitype_key = 12 
#and testg.accid in (
#'MGI:4946295',
#'MGI:2173405',
#'MGI:3776495',
#'MGI:5449569',
#'MGI:5293787',
#'MGI:4948663',
#'MGI:5297696',
#'MGI:3720108',
#'MGI:6414594',
#'MGI:5567826'
#)
#)
#'''

# rule #2/crm160
#testSQL = '''
#and exists (select 1 from ACC_Accession testg 
#where gag._genotype_key = testg._object_key and testg._mgitype_key = 12 
#and testg.accid in (
#'MGI:4948663',
#'MGI:5439284',
#'MGI:6294154',
#'MGI:5755142',
#'MGI:6414594',
#'MGI:5312862',
#'MGI:4830769',
#'MGI:3653173',
#'MGI:5521544',
#'MGI:5448443',
#'MGI:5523279',
#'MGI:5312862',
#'MGI:4830769',
#'MGI:3712071',
#'MGI:4361923',
#'MGI:5312862',
#'MGI:4830769',
#'MGI:5297696',
#'MGI:5567826',
#'MGI:6414594',
#'MGI:5003460',
#'MGI:3721552',
#'MGI:5288490',
#'MGI:5003501',
#'MGI:5430308'
#)
#)
#'''

# rule #3/crm192
#testSQL = '''
#and exists (select 1 from ACC_Accession testg 
#where gag._genotype_key = testg._object_key and testg._mgitype_key = 12 
#and testg.accid in (
#'MGI:5297696',
#'MGI:3721552',
#'MGI:5430308',
#'MGI:4822407',
#'MGI:2388127',
#'MGI:5610007',
#'MGI:5634906',
#'MGI:3717464',
#'MGI:5521546',
#'MGI:4415690',
#'MGI:5288490',
#'MGI:6693445',
#'MGI:5487451',
#'MGI:6710974',
#'MGI:3036838',
#'MGI:3052475',
#'MGI:3805456',
#'MGI:2663960',
#'MGI:5527455'
#)
#)
#'''

# rule #9/crm204
#testSQL = '''
#and exists (select 1 from ACC_Accession testg 
#where gag._genotype_key = testg._object_key and testg._mgitype_key = 12 
#and testg.accid in (
#'MGI:5567826',
#'MGI:6693445',
#'MGI:5529093',
#'MGI:6377632',
#'MGI:6192446',
#'MGI:6403447',
#'MGI:5605719',
#'MGI:3623489',
#'MGI:3810360',
#'MGI:5757706'
#)
#)
#'''

###--- classes ---###

class KeyMap:
        # Is: a memory-based cache of keys and values
        # Has: a mapping from keys to their values
        # Does: populates the mapping from the database and provides easy
        #	access to look up the value for each key

        def __init__ (self, query, keyField, valueField):
                # instantiate a new key mapping using the given SQL 'query',
                # using values for 'keyField' as the keys and values for
                # 'valueField' for their corresponding values

                self.mapping = {}

                for row in db.sql(query, 'auto'):
                        self.mapping[row[keyField]] = row[valueField]
                return

        def __len__ (self):
                # return the number of key/value pairs cached
                return len(self.mapping)

        def get (self, key):
                # return the value corresponding to the given key, or None
                # if there is no value for that 'key'

                if key in self.mapping:
                        return self.mapping[key]
                return None

class Marker:
        # Is: a single marker from the database
        # Has: a marker key and sets of annotations, evidence, evidence
        #	properties, and evidence notes that can be rolled up to that
        #	marker.
        # Does: converts primary keys and annotation types so the data can
        #	be extracted from the object for loading into the database
        #	with no primary key conflicts
        # Notes: Once you call a get*() method to extract data from this
        #	object, you can no longer call any set*() methods.  This is
        #	due to the need to prepare the data (as noted above).

        def __init__ (self, markerKey):
                # constructor; initializes object for marker with given key

                self.finalized = False

                # input data to be populated from outside this object:

                self.markerKey = markerKey

                # list of annotation rows
                self.annotations = []

                # annotation key -> list of evidence rows
                self.evidence = {}

                # evidence key -> list of property rows
                self.evidenceProperties = {}

                # evidence key -> { note key : { note record } }
                self.notes = {}

                # output data to be computed by finalize() method, in columns
                # expected by annotload product 
                self.finalAnnotations = []
                return

        def _concatenateNotes (self, evidenceKey):
                # concatenate the various notes together from noteDict
                # for the given 'evidenceKey'

                notes = ''
                if evidenceKey in self.notes:
                        for (noteKey, noteRow) in list(self.notes[evidenceKey].items()):
                                if notes:
                                        notes = notes.strip() + ' '
                                notes = notes + noteRow['note']

                return notes.replace('\n', ' ').replace('\t', ' ').strip()

        def _buildPropertiesValue (self, evidenceKey):
                # build a str.to encapsulate the various properties in the
                # manner expected by the annotload

                stanzas = []

                if evidenceKey in self.evidenceProperties:
                        lastStanza = None

                        for row in self.evidenceProperties[evidenceKey]:
                                stanza = row['stanza']

                                # if we have a change in stanza, then add a
                                # new stanza to our list of stanzas

                                if lastStanza != stanza:
                                        lastStanza = stanza
                                        stanzas.append([])

                                # compose our clause and add it to the most
                                # recent stanza

                                clause = '%s&=&%s' % ( PROPERTY_MAP.get(row['_PropertyTerm_key']), row['value'] )
                                stanzas[-1].append(clause)

                # finally, use &==& to separate clauses within a stanza, and
                # use &===& to separate the stanzas

                x = []
                for stanza in stanzas:
                        x.append('&==&'.join(stanza))

                return '&===&'.join(x)

        def finalize (self):
                # The finalize() method takes this marker object from its
                # original state (old primary keys, old annotation types) and
                # converts them to rows appropriate to be loaded as new
                # records by the annotation loader.  This can be called
                # multiple times, as it will be a no-op if the marker is
                # already finalized.

                if self.finalized:
                        return

                # The input file for the annotation loader allows up to eleven
                # fields per line.  We will generate rows in that format for
                # self.finalAnnotations, specifically for our data set:
                #    1. vocab term ID
                #    2. marker ID
                #    3. J: num
                #    4. evidence code abbreviation
                #    5. inferred from
                #    6. qualifier
                #    7. username
                #    8. empty -- use the current date by default
                #    9. notes -- can be empty
                #   10. empty -- defaults to MGI IDs for markers
                #   11. properties -- can be empty

                self.finalAnnotations = []

                for annotRow in self.annotations:
                        annotKey = annotRow['_Annot_key']

                        if annotRow['_Term_key'] == NO_PHENOTYPIC_ANALYSIS:
                                # skip annotations to this term
                                continue

                        termID = TERM_MAP.get(annotRow['_Term_key'])
                        markerID = MARKER_MAP.get(annotRow['_Marker_key'])
                        qualifier = QUALIFIER_MAP.get( annotRow['_Qualifier_key'])

                        if qualifier == None:
                                qualifier = ''

                        if annotKey not in self.evidence:
                                continue

                        for evidRow in self.evidence[annotKey]:
                                evidKey = evidRow['_AnnotEvidence_key']
                                inferredFrom = evidRow['inferredFrom']
                                evidenceCode = EVIDENCE_MAP.get( evidRow['_EvidenceTerm_key'])
                                jnumID = JNUM_MAP.get( evidRow['_Refs_key'])
                                user = USER_MAP.get(evidRow['_ModifiedBy_key'])

                                if inferredFrom == None:
                                        inferredFrom = ''

                                if qualifier == None:
                                        qualifier = ''

                                # concatenate all notes together into a single
                                # mega-note, as the annotation loader can only
                                # load a single General note and these are
                                # only for searching anyway.

                                notes = self._concatenateNotes(evidKey)

                                # build properties str.to include any
                                # properties

                                properties = self._buildPropertiesValue(evidKey)

                                if DEBUG:
                                        row = [
                                                termID,
                                                markerID
                                                ]
                                else:
                                        row = [
                                                termID,
                                                markerID,
                                                jnumID,
                                                evidenceCode,
                                                inferredFrom,
                                                qualifier,
                                                user,
                                                '',
                                                notes,
                                                '',
                                                properties
                                                ]

                                self.finalAnnotations.append(row) 

                # finished translating old records to new records
                self.finalized = True
                return

        def checkWriteable (self, setWhat):
                if self.finalized:
                        s = 'Cannot set %s on finalized marker: %d' % (setWhat, self.markerKey)
                        raise Error(s)
                return

        def setAnnotations (self, annotations):
                self.checkWriteable('annotations')
                self.annotations = annotations
                return

        def setEvidence (self, evidence):
                self.checkWriteable('evidence')
                self.evidence = evidence
                return

        def setProperties (self, properties):
                self.checkWriteable('properties')
                self.evidenceProperties = properties
                return

        def setNotes (self, notes):
                self.checkWriteable('notes')
                self.notes = notes
                return

        def getAnnotations (self):
                # gets a list of annotation rows, suitable for loading by the
                # annotation loader
                self.finalize()
                return self.finalAnnotations

###--- private functions ---###

def _stamp (s):
        # log message 's' to the global profiler
        #print(mgi_utils.date(), s)
        print(s)
        return

def _stampResults(s, order = ''):
        if DEBUG:
                results = db.sql('select * from %s %s' % (s, order), 'auto')
                for r in results:
                        _stamp(r)
                _stamp('\n')

def _addKeeper():
        if DEBUG:
                _stamp('after add to genotype_keepers')
                _stampResults('genotype_keepers', 'order by gaccid, maccid')

def _deleteKeeper():
        if DEBUG:
                _stamp('after delete from genotype_keepers')
                _stampResults('genotype_keepers', 'order by gaccid, maccid')

def _deleteScratchpad():
        if DEBUG:
                _stamp('after delete from scratchpad')
                _stampResults('scratchpad', 'order by gaccid, maccid')

        # re-create genotype_pair_counts from scratchpad
        _countAllelePairsPerGenotype('scratchpad')

def _getCount(table):
        results = db.sql('select count(1) as get_count from %s' % table, 'auto')
        if not results:
                return 0
        return results[0]['get_count']

def _identifyExpressesComponent():
        # builds a has_expresses_component temp table with (genotype key, # allele key) pairs 
        # for those genotypes and alleles with 'expresses  component' relationships.  
        # For cases where this temp table is used, it doesn't matter if the relationship is to a mouse gene or to an orthologous gene.

        _stamp('\n1:_identifyExpressesComponent/has_expresses_component')

        cmd = '''
                select distinct aa.accid, gag._Genotype_key, gag._Allele_key
                into temp table has_expresses_component
                from VOC_Annot a, GXD_AlleleGenotype gag, MGI_Relationship mr, ACC_Accession aa
                where a._AnnotType_key = %d
                        and a._Term_key != %d
                        and a._Object_key = gag._Genotype_key
                        and gag._Allele_key = mr._Object_key_1
                        and mr._Category_key = 1004
                        and gag._Genotype_key = aa._Object_key
                        and aa._MGIType_key = 12
                        and aa._logicaldb_key = 1
                        and aa.prefixpart = 'MGI:'
                        %s
                ''' % (CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS, testSQL)

        db.sql(cmd, None)
        db.sql('create index hec1 on has_expresses_component (_Genotype_key)', None)
        db.sql('create index hec2 on has_expresses_component (_Allele_key)', None)
        _stamp('1:has_expresses_component : %d\n' % (_getCount('has_expresses_component')))
        _stampResults('has_expresses_component')

        return

def _identifyMutationInvolves():
        # builds a has_mutation_involves temp table with (genotype key, allele key) pairs 
        # for those genotypes and alleles with 'mutation involves' relationships.  

        _stamp('1a:_identifyMutationInvolves/has_mutation_involves')

        cmd = '''
                select distinct aa.accid, gag._Genotype_key, gag._Allele_key, mr._Object_key_2 as _Marker_key, m.symbol
                into temp table has_mutation_involves
                from VOC_Annot a, GXD_AlleleGenotype gag, MGI_Relationship mr, ACC_Accession aa, MRK_Marker m
                where a._AnnotType_key = %d
                        and a._Term_key != %d
                        and a._Object_key = gag._Genotype_key
                        and gag._Allele_key = mr._Object_key_1
                        and mr._Category_key = 1003
                        and gag._Genotype_key = aa._Object_key
                        and aa._MGIType_key = 12
                        and aa._logicaldb_key = 1
                        and aa.prefixpart = 'MGI:'
                        and mr._Object_key_2 = m._Marker_key
                        %s
                ''' % (CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS, testSQL)

        db.sql(cmd, None)
        db.sql('create index hmi1 on has_mutation_involves (_Genotype_key)', None)
        db.sql('create index hmi2 on has_mutation_involves (_Allele_key)', None)
        db.sql('create index hmi3 on has_mutation_involves (_Marker_key)', None)
        _stamp('1a:mutation_involves : %d\n' % (_getCount('has_mutation_involves')))
        _stampResults('has_mutation_involves')

        return

def _countAllelePairsPerGenotype(table):
        # collect into a temp table (genotype_pair_counts) the genotypes that
        # have annotations attached, along with the count of allele pairs for each genotype

        _stamp('2:_countAllelePairsPerGenotype')

        db.sql('drop table if exists genotype_pair_counts', None);

        if (table == 'GXD_AllelePair'):
                cmd = '''
                        select gag._Genotype_key, count(1) as pair_count
                        into temp table genotype_pair_counts
                        from GXD_AlleleGenotype gag
                        where exists (select 1 from VOC_Annot v
                                where v._AnnotType_key in (%d)
                                and v._Term_key != %d
                                and v._Object_key = gag._Genotype_key
                                )
                        %s
                        group by gag._Genotype_key
                        ''' % (CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS, testSQL)
        else:
                cmd = '''
                        select gag._Genotype_key, count(1) as pair_count
                        into temp table genotype_pair_counts
                        from scratchpad gag
                        group by gag._Genotype_key
                        '''

        db.sql(cmd, None)
        db.sql('create index tmp_by_count on genotype_pair_counts (pair_count, _Genotype_key)', None)
        _stamp('2:genotype_pair_counts : %d\n' % _getCount('genotype_pair_counts'))

        return

def _buildKeepersTable():
        # build the genotype_keepers table, which will contain the genotype/
        # marker pairs where we can identify a single causative marker for the
        # genotype's annotations.  Note that the table will not be populated
        # by this method.

        _stamp('3:_buildKeepersTable')

        cmd = '''
                create temp table genotype_keepers (
                        _Genotype_key int not null,
                        genotype_type text null,
                        _Marker_key int null,
                        symbol text null,
                        gaccid text not null,
                        maccid text not null
                )
                '''
        db.sql(cmd, None)
        _stamp('\n')

        return

def _keepNaturallySimpleGenotypes():
        # add to the keepers table the naturally simple genotypes:
        #
        # Rule #1
        # One Marker Genotype. Naturally simple genotypes are included in the roll-up, 
        # so we associate these genotypes with their corresponding markers
        # 1. Have only one marker
        # 2. Have no "inserted expressed sequence” attribute
        # 3. Have no “mutation involves” relationships
        # 4. Are not conditional genotypes
        # 5. Ignore wild-type alleles
        #

        _stamp('4:_keepNaturallySimpleGenotypes/genotype_keepers/rule #1 : one marker genotype')
        _stamp('4.1:only one marker')
        _stamp('4.2:has_mutation_involves = false')
        _stamp('4.3:allele attribute "inserted expressed sequence" = false')
        _stamp('4.4:isConditional = false')
        _stamp('4.5:wildtype = false')

        cmd = '''
                insert into genotype_keepers
                select distinct g._Genotype_key,
                        'rule #1 : one marker genotype',
                        p._Marker_key, m.symbol, a1.accid, a2.accid
                from genotype_pair_counts g,
                        GXD_AlleleGenotype p,
                        GXD_Genotype gg,
                        MRK_Marker m,
                        ACC_Accession a1,
                        ACC_Accession a2
                where g.pair_count = 1
                        and g._Genotype_key = p._Genotype_key
                        and p._Genotype_key = gg._Genotype_key
                        -- 4.4:isConditional = false
                        and gg.isConditional = 0        
                        and p._Marker_key = m._Marker_key

                        and g._Genotype_key = a1._Object_key
                        and a1._MGIType_key = 12
                        and a1._logicaldb_key = 1
                        and a1.prefixpart = 'MGI:'
                        and p._Marker_key = a2._Object_key
                        and a2._MGIType_key = 2
                        and a2._logicaldb_key = 1
                        and a2.prefixpart = 'MGI:'
                        and a2.preferred = 1

                        -- 4.2:has_mutation_involves = false
                        and not exists (select 1 from has_mutation_involves mi where g._Genotype_key = mi._Genotype_key)

                        -- 4.3:allele attribute "inserted expressed sequence" = false
                        and not exists (select 1 from VOC_Annot b -- allele subtype annotation
                                where b._AnnotType_key = 1014   
                                and p._Allele_key = b._Object_key -- inserted expressed sequence
                                and b._Term_key = 11025597     
                                )

                        -- 4.5:wildtype = false')
                        and exists (select 1
                                from GXD_AlleleGenotype gag, ALL_Allele a
                                where g._Genotype_key = gag._Genotype_key
                                and gag._Allele_key = a._Allele_key
                                and a.isWildType = 0
                                )
                '''

        db.sql(cmd, None)
        _stamp('4:add naturally simple genotypes: %d\n' % _getCount('genotype_keepers'))
        _addKeeper()

        return

def _indexKeepersTable():
        # add relevant indexes to the genotype_keepers table
        _stamp('5:_indexKeepersTable\n\n')
        db.sql('create index gk_genotype on genotype_keepers (_Genotype_key)', None)
        db.sql('create index gk_marker on genotype_keepers (_Marker_key, _Genotype_key)', None)
        return

def _identifyReporterTransgenes():
        # collect into a temp table (reporter_transgenes) the alleles that
        # are reporter transgenes, defined as alleles with:
        #   1. allele type/generation type "Transgenic" = true
        #   2. allele subtype/attribute "Reporter" = true
        #   3. allele subtype/attribute no other selected

        _stamp('7:_identifyReporterTransgenes/reporter_transgenes')
        _stamp('7.1. allele type/generation type "Transgenic" = true')
        _stamp('7.2. allele subtype/attribute "Reporter" = true')
        _stamp('7.3. allele subtype/attribute no other subtype selected')

        cmd = '''
                select distinct a._Allele_key
                into temp table reporter_transgenes
                from ALL_Allele a
                where a._Allele_Type_key = 847126              
                        and exists (select 1 from VOC_Annot v   -- transgenic
                                where a._Allele_key = v._Object_key -- allele subtype annotation
                                and v._AnnotType_key = 1014         -- reporter
                                and v._Term_key = 11025589      
                                )
                        and not exists (select 1 from VOC_Annot v
                                where a._Allele_key = v._Object_key -- allele subtype annotation
                                and v._AnnotType_key = 1014         -- not reporter
                                and v._Term_key != 11025589     
                                )
                '''

        db.sql(cmd, None)
        db.sql('create unique index tmp_reportertg on reporter_transgenes (_Allele_key)', None)
        _stamp('7:built reporter_transgenes table rows: %d\n' % _getCount('reporter_transgenes'))

        return

def _identifyTransactivators():
        # collect into a temp table (transactivators) the alleles that:
        #   1. are of allele type (generation type) "Transgenic"
        #   2. have an attribute (subtype) of "Transactivator"
        #   3. do NOT have an "inserted expressed sequence" attribute

        # The "distinct" should be superfluous, but we'll include it just in
        # case something flaky comes up.

        _stamp('10:_identifyTransactivators/transactivators')
        _stamp('10.1: allele type/generation type "Transgenic" = true')
        _stamp('10.2: allele subtype/attribute "Transactivator" = true')
        _stamp('10.3: allele subtype/attribute "Inserted_expressed_sequence" = false')

        cmd = '''
                select distinct a._Allele_key
                into temp table transactivators
                from all_allele a, voc_annot t
                where a._Allele_Type_key = 847126         -- transgenic
                        and exists (select 1 from VOC_Annot v
                                where a._Allele_key = v._Object_key
                                and v._AnnotType_key = 1014    -- allele subtype annotation
                                and v._Term_key = 13289567     -- transactivator
                        )
                        and not exists (select 1 from VOC_Annot v
                                where a._Allele_key = v._Object_key
                                and v._AnnotType_key = 1014     -- allele subtype annotation
                                and v._Term_key = 11025597      -- inserted_expressed_sequence
                                )
                '''

        db.sql(cmd, None)
        db.sql('create unique index tmp_transactivators on transactivators (_Allele_key)', None)
        _stamp('10:built transactivators table rows: %d\n' % _getCount('transactivators'))

        return

def _buildScratchPad():
        # For genotypes with multiple allele pairs, we need to try to apply 
        # some rules to nail things down to a single causative marker for each.
        # To begin, we'll collect a table of genotype/marker/allele data to
        # use as a scratch pad for further calculations.

        _stamp('6:_buildScratchPad/scratchpad')
        _stamp('6:exists in genotype_pair_counts')
        _stamp('6:does not exist in genotype_keepers')

        cmd = '''
                select distinct gag._Genotype_key,
                        gag._Marker_key,
                        gag._Allele_key,
                        g.isConditional,
                        m.symbol,
                        a1.accid as gaccid,
                        a2.accid as maccid
                into temp table scratchpad
                from genotype_pair_counts c,
                        GXD_AlleleGenotype gag,
                        GXD_Genotype g,
                        MRK_Marker m,
                        ACC_Accession a1,
                        ACC_Accession a2
                where not exists (select 1 from genotype_keepers k 
                        where c._Genotype_key = k._Genotype_key)
                        and c._Genotype_key = gag._Genotype_key
                        and c._Genotype_key = g._Genotype_key
                        and gag._Marker_key = m._Marker_key
                        and gag._Genotype_key = a1._Object_key
                        and a1._MGIType_key = 12
                        and a1._logicaldb_key = 1
                        and a1.prefixpart = 'MGI:'
                        and gag._Marker_key = a2._Object_key
                        and a2._MGIType_key = 2
                        and a2._logicaldb_key = 1
                        and a2.prefixpart = 'MGI:'
                        and a2.preferred = 1
                '''

        db.sql(cmd, None)
        db.sql('create index scratch_alleles on scratchpad (_Allele_key)', None)
        _stamp('6:built scratchpad table rows: %d\n' % _getCount('scratchpad'))
        _deleteScratchpad()

        return

def _removeConditionalGenotypes():
        # If the genotype is conditional, then exclude recombinase alleles that do not have subtype = "Inserted expressed sequence"

        _stamp('8:_removeConditionalGenotypes/scratchpad')
        _stamp('8.1:isConditional = true')
        _stamp('8.2:allele attribute Recombinase = true')
        _stamp('8.3:allele attribute "inserted expressed sequence" = false')

        before = _getCount('scratchpad')

        cmd = '''
                delete from scratchpad p
                -- 8.1:isConditional = true
                where p.isConditional = 1

                        -- 8.2:allele attribute Recombinase = true
                        and exists (select 1 from VOC_Annot v
                                where p._Allele_key = v._Object_key
                                and v._AnnotType_key = 1014       -- allele subtype annotation
                                and v._Term_key = 11025588        -- recombinase
                                )

                        -- 8.3:allele attribute "inserted expressed sequence" = false
                        and not exists (select 1 from VOC_Annot v
                                where p._Allele_key = v._Object_key
                                and v._AnnotType_key = 1014    -- allele subtype annotation
                                and v._Term_key = 11025597     -- inserted expressed sequence
                                )

               '''

        db.sql(cmd, None)
        _stamp('8:delete recombinase alleles from scratchpad: %d\n' % (before - _getCount('scratchpad'))) 
        _deleteScratchpad()

        return

def _cleanupTempTables():
        # drop any temp tables that we're done with

        tables = [ 
                'has_expresses_component',
                'has_mutation_involves',
                'genotype_pair_counts',
                'genotype_keepers',
                'reporter_transgenes',
                'transactivators',
                'scratchpad',
                'wildtype_alleles',
                'trad',
                'trad_ct',
                'mi',
                'mi_ct',
                'ec',
                'ec_ct',
                'mi1',
                'mi2',
                'mi3'
                ]

        for table in tables:
                db.sql('drop table if exists %s' % table, None)
                _stamp('Drop temp table: %s' % table)
        db.commit()

        return

def _removeReporterTransgenes():
        # exclude reporter transgene alleles (allele type=Transgenic, allele subtype = "reporter" AND NOT "inserted expressed sequence")

        _stamp('9:_removeReporterTransgenes/scratchpad')
        _stamp('9:reporter_transgenes = true')
        before = _getCount('scratchpad')
        #_stampResults('reporter_transgenes')
        db.sql('delete from scratchpad p where exists (select 1 from reporter_transgenes r where p._allele_key = r._allele_key)', None)
        _stamp('9:delete reporter transgenes from scratchpad: %d\n' % (before - _getCount('scratchpad')))
        _deleteScratchpad()

        return 

def _removeTransactivators():
        # remove from the scratch pad any alleles which are transactivators (and are transgenic)

        _stamp('11:_removeTransactivators/scratchpad')
        _stamp('11:transactivators = true')
        #_stampResults('transactivators')
        before = _getCount('scratchpad')
        db.sql('delete from scratchpad p where exists (select 1 from transactivators t where p._allele_key = t._allele_key)', None)
        _stamp('11:delete transactivators from scratchpad: %d\n' % ( before - _getCount('scratchpad')) )
        _deleteScratchpad()

        return 

def _identifyWildTypeAlleles():
        # Identify which remaining alleles are wild-type alleles and include them in a temp table.

        _stamp('12:_identifyWildTypeAlleles/wildtype_alleles')
        _stamp('12:isWildType = true')

        cmd = '''
                select distinct p._Allele_key
                into temp table wildtype_alleles
                from scratchpad p, ALL_Allele a
                where p._Allele_key = a._Allele_key and a.isWildType = 1
               '''

        db.sql(cmd, None)
        db.sql('create unique index wt_alleles on wildtype_alleles (_Allele_key)', None)
        _stamp('12:wildtype_alleles table rows: %d\n' % _getCount('wildtype_alleles'))

        return

def _removeWildTypeAllelesFromScratchPad():
        # Remove any remaining wild-type alleles.

        _stamp('13:_removeWildTypeAllelesFromScratchPad/scratchpad')
        _stamp('13:wildtype_alleles = true')
        before = _getCount('scratchpad')
        _stampResults('wildtype_alleles')
        db.sql('delete from scratchpad p where exists (select 1 from wildtype_alleles w where p._allele_key = w._allele_key)', None)
        _stamp('13:delete wild-type from scratchpad: %d\n' % ( before - _getCount('scratchpad')) )
        _deleteScratchpad()

        return

def _collectMarkerSets():
        # builds several temp tables with markers tied to genotypes by three
        # routes:  traditional marker/allele pairings, expressed component
        # relationships, and mutation involves relationships.  Note that the
        # temp tables will only contain data for genotypes (and only through
        # alleles) which still exist in scratchpad.  Also note that these
        # only include genotypes which have annotations of the current
        # annotation type.  Still another note -- we exclude genotypes where
        # the only annotation is to 'no phenotypic analysis'.

        # The temp tables we will build include:
        # 1. trad - genotype-to-marker via the traditional allele-marker route
        # 2. mi - genotype-to-marker via the mutation involves route
        # 3. ec - genotype-to-marker via the expresses component route
        # 4. trad_ct - genotype to count of its markers in trad
        # 5. mi_ct - genotype to count of its markers in mi
        # 6. ec_ct - genotype to count of its markers in ec

        _stamp('14:_collectMarkerSets')

        # command to build table 1 (distinct genotype/marker pairs)

        tradCmd = '''
                select distinct s._Genotype_key, s._Marker_key
                into temp table trad
                from scratchpad s, ALL_Allele a
                where s._Allele_key = a._Allele_key
                '''

        # commands to build tables 2 & 3 (distinct genotype/marker pairs)

        miCmd = '''
                select distinct s._Genotype_key, mr._Object_key_2 as _Marker_key
                into temp table mi
                from scratchpad s, MGI_Relationship mr
                where s._Allele_key = mr._Object_key_1 
                and mr._Category_key = 1003      -- MI relationship
                '''

        # We also need to track the organism of each expressed marker, as
        # there are special cases down the road which require mouse-only.

        ecCmd = '''
                select distinct s._Genotype_key, mr._Object_key_2 as _Marker_key, m._Organism_key, mr._RelationshipTerm_key
                into temp table ec
                from scratchpad s, MGI_Relationship mr, MRK_Marker m
                where s._Allele_key = mr._Object_key_1
                and mr._Category_key = 1004     -- EC relationship
                and mr._Object_key_2 = m._Marker_key
                '''

        # commands to index tables 1-3

        templateB = 'create index %s on %s (_Genotype_key)'

        tradB = templateB % ('tradIndex1', 'trad')
        ecB = templateB % ('ecIndex1', 'ec')
        miB = templateB % ('miIndex1', 'mi')

        # commands to build tables 4-6 (counts of distinct markers/genotype)

        templateC = '''select _Genotype_key, count(1) as marker_count into temp table %s from %s group by _Genotype_key'''

        tradC = templateC % ('trad_ct', 'trad')
        ecC = templateC % ('ec_ct', 'ec')
        miC = templateC % ('mi_ct', 'mi')

        # commands to index tables 4-6

        tradD = templateB % ('tradIndex2', 'trad_ct')
        ecD = templateB % ('ecIndex2', 'ec_ct')
        miD = templateB % ('miIndex2', 'mi_ct')
        
        templateE = 'create index %s on %s (marker_count)'

        tradE = templateE % ('tradIndex3', 'trad_ct')
        ecE = templateE % ('ecIndex3', 'ec_ct')
        miE = templateE % ('miIndex3', 'mi_ct')

        # build and index each table, and report on completion for profiling

        for (name, tbl1, tbl2, c1, c2, c3, c4, c5) in [
                ('traditional', 'trad', 'trad_ct', tradCmd, tradB, tradC, tradD, tradE),
                ('mutation involves', 'mi', 'mi_ct', miCmd, miB, miC, miD, miE),
                ('expresses component', 'ec', 'ec_ct', ecCmd, ecB, ecC, ecD, ecE)
                ]:

                db.sql(c1, None)
                db.sql(c2, None)
                db.sql(c3, None)
                db.sql(c4, None)
                db.sql(c5, None)

                _stamp('14:built table of %s genotype/marker pairs rows: %d' % (name, _getCount(tbl1)) )
                _stamp('14:built table of %s genotype/marker counts rows: %d' % (name, _getCount(tbl2)) )

                if DEBUG:
                        results = db.sql('''
                                select a.accid, m.symbol, t.* 
                                from %s t, MRK_Marker m, ACC_Accession a
                                where t._Marker_key = m._Marker_key
                                and m._Marker_key = a._Object_key
                                and a._MGIType_key = 2
                                and a._logicaldb_key = 1
                                and a.prefixpart = 'MGI:'
                                and a.preferred = 1
                                order by a.accid
                                ''' % (tbl1), 'auto')
                        for r in results:
                                _stamp(r)
         
                _stamp('\n')

        db.sql('create index ecOrg on ec (_Organism_key)', None)
        db.sql('create index ecTerm on ec (_RelationshipTerm_key)', None)

        return

def _handleMultipleMarkers():
        # for genotypes where multiple markers remain, we only allow one rule
        # by which to roll up annotations:
        #	if there are no markers through 'mutation involves' AND
        #		there are exacty two markers AND
        #		one is a transgene AND
        #		the other is a non-transgene AND
        #		the transgene 'expresses' only the non-transgene
        #			through an 'expresses component' relationship
        #	then we roll up annotations to both markers
        # This function leaves scratchpad with only genotypes having a single
        # marker due to traditional marker-to-allele pairings.

        #EXPRESSES_MOUSE_GENE = 12965808	# term key for 'expresses_mouse_gene'

        _stamp('15:_handleMultipleMarkers/rule #2 : transgene, 1 EC, 0 MI')

        before = _getCount('genotype_keepers')

        template = '''
                insert into genotype_keepers
                select distinct s._Genotype_key, \'rule #2 : transgene, 1 EC, 0 MI\', %s, %s, a1.accid, a2.accid
                from scratchpad s,
                        trad_ct tc, 
                        trad tt, 
                        MRK_Marker mt,
                        trad tn,
                        MRK_Marker nt,
                        ACC_Accession a1,
                        ACC_Accession a2
                where 
                        s._Genotype_key = a1._Object_key
                        and a1._MGIType_key = 12
                        and a1._logicaldb_key = 1
                        and a1.prefixpart = 'MGI:'
                        and %s = a2._Object_key
                        and a2._MGIType_key = 2
                        and a2.prefixpart = 'MGI:'
                        and a2._logicaldb_key = 1
                        and a2.preferred = 1

                        -- no 'mutation involves'
                        and not exists (select 1 from mi_ct mc where s._Genotype_key = mc._Genotype_key)

                        -- exactly two markers
                        and s._Genotype_key = tc._Genotype_key
                        and tc.marker_count = 2

                        -- one transgene
                        and s._Genotype_key = tt._Genotype_key
                        and tt._Marker_key = mt._Marker_key
                        and mt._Marker_Type_key = 12    -- transgene

                        -- one non-transgene
                        and s._Genotype_key = tn._Genotype_key
                        and tn._Marker_key = nt._Marker_key
                        and nt._Marker_Type_key != 12   -- transgene

                        -- transgene expresses mouse non-transgene
                        and exists (select 1
                                from MGI_Relationship r1, ALL_Allele a
                                where r1._Category_key = 1004   -- EC relationship
                                and mt._Marker_key = a._Marker_key
                                and a._Allele_key = r1._Object_key_1
                                and r1._RelationshipTerm_key = 12965808
                                and r1._Object_key_2 = nt._Marker_key
                                )

                        -- transgene does not express any other marker
                        and not exists (select 1
                                from MGI_Relationship r2, ALL_Allele a2
                                where r2._Category_key = 1004   -- EC relationship
                                and mt._Marker_key = a2._Marker_key
                                and a2._Allele_key = r2._Object_key_1
                                and r2._Object_key_2 != nt._Marker_key
                                )

                 ''' % ('%s', '%s', '%s')

        transgeneCmd = template % ('mt._Marker_key', 'mt.symbol', 'mt._Marker_key')
        otherCmd = template % ('nt._Marker_key', 'nt.symbol', 'nt._Marker_key')
        db.sql(transgeneCmd, None)
        db.sql(otherCmd, None)
        _stamp('15:add rows to genotype_keepers for transgene rule A: %d' % (_getCount('genotype_keepers') - before))

        beforeSP = _getCount('scratchpad')
        db.sql('delete from scratchpad where _Genotype_key in (select _Genotype_key from trad_ct where marker_count > 1)', None)
        _stamp('15:delete multi-marker genotypes from scratchpad: %d\n' % ( beforeSP - _getCount('scratchpad')))

        return

def _handleMutationInvolves():
        # see below

        _stamp('16:_handleMutationInvolves/rule #3 : mutation involves')

        results = db.sql('select * from mi_ct', 'auto')
        _stamp('select * from mi_ct: %d\n' % len(results))
        #_stamp(results)

        #
        # important : check scratchpad by genotype AND allele
        #

        #
        # rule #3 Non-Transgene clause
        #
        cmd = '''
                select distinct s.gaccid, s.maccid, s._Genotype_key, s._Marker_key, s.symbol
                into temp table mi1
                from genotype_pair_counts c, scratchpad s, has_mutation_involves mi
                where c.pair_count = 1
                        and c._Genotype_key = s._Genotype_key
                        and s._Genotype_key = mi._Genotype_key
                        and s._Allele_key = mi._Allele_key
                        and (
                             exists (select 1 from VOC_Annot v
                                where s._Marker_key = v._Object_key
                                and v._Annottype_key = 1011
                                and v._Term_key in (6238170,97015607,97015607,103059157,103059158,103059155,15406207)
                                )
                             or exists (select 1 from MRK_Marker m
                                where s._Marker_key = m._Marker_key
                                and m._Marker_Type_key in (3,10)
                                )
                        )
                '''

        db.sql(cmd, None)
        db.sql('create index mi1_1 on mi1 (_Genotype_key)', None)
        db.sql('create index mi1_2 on mi1 (_Marker_key)', None)
        _stamp('16:add rows to mi1 rule#3/Non-Transgene clause')
        _stamp('16:1:the genotype has exactly 1 marker in (M)')
        _stamp('16:2:the genotype has at least 1 mutation involves marker in (I)')
        _stamp('16:3:marker (M) feature type = heritable phenotypic marker, enhancer, silencer, imprinting control region, locus control region, promoter')
        _stamp('16:3:OR marker type = complex/cluster/region) OR cytogenetic marker')

        if DEBUG:
                _stamp('\nmi1')
                results = db.sql('select * from mi1', 'auto')
                for r in results:
                        _stamp(r)

        #
        # rule #3 Transgene clause
        #
        cmd = '''
                select distinct s.gaccid, s.maccid, s._Genotype_key, s._Marker_key, s.symbol
                into temp table mi2
                from genotype_pair_counts c, scratchpad s, has_mutation_involves mi
                where c.pair_count = 1
                        and c._Genotype_key = s._Genotype_key
                        and s._Genotype_key = mi._Genotype_key
                        and s._Allele_key = mi._Allele_key
                        and exists (select 1 from mi_ct where s._Genotype_key = mi_ct._Genotype_key and mi_ct.marker_count = 1)
                        and exists (select 1 from MRK_Marker m
                                where s._Marker_key = m._Marker_key
                                and m._Marker_Type_key in (12)
                        )
                union
                select distinct s.gaccid, s.maccid, s._Genotype_key, mi._Marker_key, mi.symbol
                from genotype_pair_counts c, scratchpad s, has_mutation_involves mi
                where c.pair_count = 1
                        and c._Genotype_key = s._Genotype_key
                        and s._Genotype_key = mi._Genotype_key
                        and s._Allele_key = mi._Allele_key
                        and exists (select 1 from mi_ct where s._Genotype_key = mi_ct._Genotype_key and mi_ct.marker_count = 1)
                        and exists (select 1 from MRK_Marker m
                                where s._Marker_key = m._Marker_key
                                and m._Marker_Type_key in (12)
                        )
                '''
        _stamp('16:add rows to mi2 rule#3/Transgene clause')
        _stamp('16:1:the genotype has exactly 1 marker in (M)')
        _stamp('16:2:the genotype has exactly 1 mutation involves marker in (I)')
        _stamp('16:3:marker type = Transgene (12)')
        _stamp('16:4:rollup to both marker in (M) and marker in (I)')
        db.sql(cmd, None)
        db.sql('create index mi2_1 on mi2 (_Genotype_key)', None)
        db.sql('create index mi2_2 on mi2 (_Marker_key)', None)

        if DEBUG:
                _stamp('\nmi2')
                results = db.sql('select * from mi2', 'auto')
                for r in results:
                        _stamp(r)

        #
        # rule #3 Docking Site clause
        #
        docking_sites = ','.join(map(str, DOCKING_SITES))
        cmd = '''
                select s.gaccid, s.maccid, s._Genotype_key, mi._Marker_key, mi.symbol
                into temp table mi3
                from genotype_pair_counts c, scratchpad s, has_mutation_involves mi
                where c.pair_count = 1
                        and c._Genotype_key = s._Genotype_key
                        and s._Genotype_key = mi._Genotype_key
                        and s._Allele_key = mi._Allele_key
                        and exists (select 1 from mi_ct where s._Genotype_key = mi_ct._Genotype_key and mi_ct.marker_count = 1)
                        and exists (select 1 from MRK_Marker m
                                where s._Marker_key = m._Marker_key
                                and m._Marker_key in (%s)
                        )
                ''' % (docking_sites)

        _stamp('16:add rows to mi3 rule#3/Docking Site clause')
        _stamp('16:1:the genotype has exactly 1 marker')
        _stamp('16:2:the marker is Docking Site (Col1a1, Gt(ROSA)26Sor, Hprt)')
        _stamp('16:3:the genotype has exactly 1 mutation involves marker in (I)')
        db.sql(cmd, None)
        db.sql('create index mi3_1 on mi3 (_Genotype_key)', None)
        db.sql('create index mi3_2 on mi3 (_Marker_key)', None)

        if DEBUG:
                _stamp('\nmi3')
                results = db.sql('select * from mi3', 'auto')
                for r in results:
                        _stamp(r)

        before = _getCount('genotype_keepers')
        cmd = '''
                insert into genotype_keepers
                select s._Genotype_key, \'rule #3 : non-transgene clause\', s._Marker_key, s.symbol, s.gaccid, s.maccid
                from scratchpad s, mi1
                where s._genotype_key = mi1._genotype_key
                '''
        db.sql(cmd, None)
        _stamp('16:added rows to genotype_keepers from mi1: %d\n' % (_getCount('genotype_keepers') - before))
        _addKeeper()

        before = _getCount('genotype_keepers')
        cmd = '''
                insert into genotype_keepers
                select s._Genotype_key, \'rule #3 : transgene clause\', mi2._Marker_key, mi2.symbol, s.gaccid, s.maccid
                from scratchpad s, mi2
                where s._genotype_key = mi2._genotype_key
                '''
        db.sql(cmd, None)
        _stamp('16:added rows to genotype_keepers from mi2: %d\n' % (_getCount('genotype_keepers') - before))
        _addKeeper()

        before = _getCount('genotype_keepers')
        cmd = '''
                insert into genotype_keepers
                select s._Genotype_key, \'rule #3 : docking site clause\', mi3._Marker_key, mi3.symbol, s.gaccid, s.maccid
                from scratchpad s, mi3
                where s._genotype_key = mi3._genotype_key
                '''
        db.sql(cmd, None)
        _stamp('16:added rows to genotype_keepers from mi3: %d\n' % (_getCount('genotype_keepers') - before))
        _addKeeper()

        before2 = _getCount('scratchpad')
        db.sql('delete from scratchpad where _Genotype_key in (select _Genotype_key from mi_ct)', None)
        db.sql('delete from scratchpad where _Genotype_key in (select _Genotype_key from mi1)', None)
        db.sql('delete from scratchpad where _Genotype_key in (select _Genotype_key from mi2)', None)
        db.sql('delete from scratchpad where _Genotype_key in (select _Genotype_key from mi3)', None)
        _stamp('16:delete rows from scratchpad due to mutation involves rule: %d\n' % (before2 - _getCount('scratchpad')) )
        _deleteScratchpad()

        return

def _handleTransgenes():
        # Assumes: all genotypes in scratchpad have a single marker via the
        #	traditional marker-to-allele pairs.  Also assumes that we'll
        #	delete Gt(ROSA) data later.
        # This function handles the case where the associated marker is a
        # transgene, including a special case where the transgene has a single
        # expressed component relationship.

        _stamp('17:_handleTransgenes/rule #4 : transgene')
        _stamp('17:_handleTransgenes/rule #5 : transgene, 1 EC')

        # single marker is a transgene
        cmdTg = '''
                insert into genotype_keepers
                select s._Genotype_key, \'rule #4 : transgene\', s._Marker_key, s.symbol, s.gaccid, s.maccid
                from scratchpad s, MRK_Marker m
                where s._Marker_key = m._Marker_key and m._Marker_Type_key = 12 -- transgene
                '''

        # single marker is a transgene with one expressed component, also
        # include the expressed component marker
        cmdEC = '''
                insert into genotype_keepers
                select s._Genotype_key, \'rule #5 : transgene, 1 EC\', ec._Marker_key, s.symbol, s.gaccid, s.maccid
                from scratchpad s, MRK_Marker m, ec_ct ct, ec ec
                where s._Marker_key = m._Marker_key
                and m._Marker_Type_key = 12     -- transgene
                and s._Genotype_key = ec._Genotype_key
                and ec._RelationshipTerm_key = 12965808 -- expresses_mouse_gene
                and s._Genotype_key = ct._Genotype_key
                and ct.marker_count = 1
                '''

        # delete genotypes with transgene markers from scratchpad
        cmdDel = '''
                delete from scratchpad
                where _Genotype_key in (select s._Genotype_key
                        from scratchpad s, MRK_Marker m
                        where s._Marker_key = m._Marker_key
                        and m._Marker_Type_key = 12     -- transgene
                        )
                '''

        ct1 = _getCount('genotype_keepers')
        db.sql(cmdTg, None)
        ct2 = _getCount('genotype_keepers')
        _stamp('17:add rows to genotype_keepers for transgenes: %d' % ( ct2 - ct1))
        db.sql(cmdEC, None)
        _stamp('17:add rows to genotype_keepers for expressed components: %d' % (_getCount('genotype_keepers') - ct2))
        _addKeeper()

        ct3 = _getCount('scratchpad')
        db.sql(cmdDel, None)
        _stamp('17:delete rows from scratchpad for transgenes: %d\n' % ( ct3 - _getCount('scratchpad')) )
        _deleteScratchpad()

        return

def _handleDockingSites():
        # Assumes: all genotypes in scratchpad have a single marker via the
        #	traditional marker-to-allele pairs.  Also assumes that we'll
        #	delete Gt(ROSA) data later.
        # This function handles the case where the marker associated with a
        # genotype is a docking site.

        _stamp('18:_handleDockingSites/rule #6 : docking site, 1 EC')
        _stamp('18.1: marker count = 1')
        _stamp('18.1: marker in (M) is a docking site')
        _stamp('18.2: EC = expresses_mouse_gene')

        docking_sites = ','.join(map(str, DOCKING_SITES))

        cmd = '''
                insert into genotype_keepers
                select s._Genotype_key, \'rule #6 : docking site, 1 EC\', ec._Marker_key, s.symbol, s.gaccid, s.maccid
                from scratchpad s, ec_ct ct, ec ec
                where s._Marker_key in (%s)
                        and s._Genotype_key = ct._Genotype_key
                        and ct.marker_count = 1
                        and ec._RelationshipTerm_key = 12965808 -- expresses_mouse_gene
                        and s._Genotype_key = ec._Genotype_key
                ''' % (docking_sites)
        ct = _getCount('genotype_keepers')
        db.sql(cmd, None)
        _stamp('18:add rows to genotype_keepers for docking sites rule #6: %d' % (_getCount('genotype_keepers') - ct))
        _addKeeper()

        _stamp('18:_handleDockingSites/rule #9 : docking site, Expresses No Components')
        _stamp('18.1: marker in (M) is a docking site other than Gt(ROSA)26Sor')
        _stamp('18.2: allele subtype/attribute "Inserted_expressed_sequence" = false')

        cmd = '''
                insert into genotype_keepers
                select s._Genotype_key, \'rule #9 : docking site, 0 EC\', s._Marker_key, s.symbol, s.gaccid, s.maccid
                from scratchpad s
                where s._Marker_key in (%d,%d)
                        -- allele attribute "inserted expressed sequence" = false
                        and not exists (select 1 from VOC_Annot b -- allele subtype annotation
                                where b._AnnotType_key = 1014   
                                and s._Allele_key = b._Object_key -- inserted expressed sequence
                                and b._Term_key = 11025597     
                                )

                 ''' % (HPRT, COL1A1)
        ct = _getCount('genotype_keepers')
        db.sql(cmd, None)
        _stamp('18:add rows to genotype_keepers for docking sites rule #9: %d' % (_getCount('genotype_keepers') - ct))
        _addKeeper()

        cmd = '''delete from scratchpad where _Marker_key in (%s)''' % docking_sites
        ct = _getCount('scratchpad')
        db.sql(cmd, None)
        _stamp('18:delete rows from scratchpad for docking sites: %d\n' % ( ct - _getCount('scratchpad')) )
        _deleteScratchpad()

        return

def _handleOtherSingles():
        # Assumes: all genotypes in scratchpad have a single marker via the
        #	traditional marker-to-allele pairs.  Also assumes that we'll
        #	delete Gt(ROSA) data later.  And, assumes that transgenes 
        #	and docking sites have been removed from scratchpad.
        # This function handles genotypes with a single marker that is not
        # a transgene or a docking site, where that single marker may or may
        # not be the sole expressed component of itself.

        _stamp('19:_handleOtherSingles/rule #7 : singles, no EC')
        _stamp('19:_handleOtherSingles/rule #8 : self-expressing single')

        # singles with no expressed components
        cmd1 = '''
                insert into genotype_keepers
                select s._Genotype_key, \'rule #7 : singles, no EC\', s._Marker_key, s.symbol, s.gaccid, s.maccid
                from scratchpad s
                where not exists (select 1 from has_expresses_component ct where s._Genotype_key = ct._Genotype_key)
                '''

        # singles where the marker knows how to express itself
        cmd2 = '''
                insert into genotype_keepers
                select s._Genotype_key, \'rule #8 : self-expressing single\', s._Marker_key, s.symbol, s.gaccid, s.maccid
                from scratchpad s, ec_ct ct, ec ec
                where s._Genotype_key = ct._Genotype_key
                        and ct.marker_count = 1
                        and ec._RelationshipTerm_key = 12965808 -- expresses_mouse_gene
                        and s._Genotype_key = ec._Genotype_key
                        and s._Marker_key = ec._Marker_key
                '''

        ct1 = _getCount('genotype_keepers')
        db.sql(cmd1, None)
        ct2 = _getCount('genotype_keepers')
        _stamp('19:add rows to genotype_keepers for singles with no EC: %d' % (ct2 - ct1))
        db.sql(cmd2, None)
        _stamp('19:add rows to genotype_keepers for self-expressing singles: %d\n' % (_getCount('genotype_keepers') - ct2))
        _addKeeper()

        return

def _removeNullsAndGtRosa():
        # finally, delete all the genotypes associated with Gt(ROSA)26Sor

        _stamp('20:_removeNullsAndGtRosa')
        ct1 = _getCount('genotype_keepers')
        db.sql('delete from genotype_keepers where _Marker_key = %d' % GT_ROSA, None)
        ct2 = _getCount('genotype_keepers')
        _stamp('20:delete genotypes associated with Gt(ROSA)26Sor: %d' % (ct1 - ct2))
        db.sql('delete from genotype_keepers where _Marker_key is null', None)
        _stamp('20:delete genotypes associated with null markers: %d\n' % (ct2 - _getCount('genotype_keepers')))
        _deleteKeeper()

        return

def _getMarkerMetaData():
        # populate global variables with counts of annotations for each marker
        # and an ordered list of marker keys to process.

        global ANNOTATION_COUNTS, MARKER_KEYS, LAST_MARKER_KEY_INDEX

        _stamp('22:_getMarkerMetaData')

        ANNOTATION_COUNTS = {}
        MARKER_KEYS = []
        LAST_MARKER_KEY_INDEX = None

        cmd = '''
                select k._Marker_key, count(1) as annotation_count
                from genotype_keepers k, VOC_Annot a
                where k._Genotype_key = a._Object_key
                and a._AnnotType_key in (%s)
                and a._Term_key != %d
                group by k._Marker_key
                ''' % (CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS)

        results = db.sql(cmd, 'auto')
        for row in results:
                ANNOTATION_COUNTS[row['_Marker_key']] = row['annotation_count']

        if DEBUG:
                for r in results:
                        _stamp(r)
                _stamp('\n')

        MARKER_KEYS = list(ANNOTATION_COUNTS.keys())
        MARKER_KEYS.sort()

        _stamp('22:Retrieved annotation counts for %d markers\n' % len(MARKER_KEYS)) 

        return

def _initializeKeyMaps():
        # initialize the mappings from various database keys to their
        # respective values

        global TERM_MAP, MARKER_MAP, JNUM_MAP, EVIDENCE_MAP
        global QUALIFIER_MAP, USER_MAP, PROPERTY_MAP, CURRENT_ANNOT_TYPE

        _stamp('23:_initializeKeyMaps')

        if not CURRENT_ANNOT_TYPE:
                raise Error('Need to call setAnnotationType()')

        # map from annotated term IDs to their IDs
        if DEBUG:
                accID = "aa.accID || ':' || t.term as accID"
        else:
                accID = "aa.accID as accID"

        termCmd = '''
                select distinct aa._Object_key, %s
                from VOC_Annot va, ACC_Accession aa, VOC_Term t
                where va._AnnotType_key in (%d)
                and va._Term_key = aa._Object_key
                and aa._MGIType_key = 13
                and aa.private = 0
                and aa.preferred = 1
                and va._Term_key = t._Term_key
                ''' % (accID, CURRENT_ANNOT_TYPE)
        TERM_MAP = KeyMap(termCmd, '_Object_key', 'accID')

        # map from annotated markers to their MGI IDs
        if DEBUG:
                accID = "aa.accID || ':' || m.symbol as accID"
        else:
                accID = "aa.accID as accID"

        markerCmd = '''
                select distinct aa._Object_key, %s
                from genotype_keepers k, ACC_Accession aa, MRK_Marker m
                where k._Marker_key = aa._Object_key
                        and aa._MGIType_key = 2
                        and aa.private = 0
                        and aa.preferred = 1
                        and aa._LogicalDB_key = 1
                        and k._Marker_key = m._Marker_key
                ''' % (accID)
        MARKER_MAP = KeyMap(markerCmd, '_Object_key', 'accID')

        # map from reference keys to their Jnum IDs
        jnumCmd = '''
                select distinct r._Refs_key, r.jnumID
                from VOC_Annot va, VOC_Evidence ve, BIB_Citation_Cache r
                where va._AnnotType_key = %d
                and va._Annot_key = ve._Annot_key
                and ve._Refs_key = r._Refs_key
                ''' % CURRENT_ANNOT_TYPE
        JNUM_MAP = KeyMap(jnumCmd, '_Refs_key', 'jnumID')

        # map from evidence term keys to their abbreviations
        evidenceCmd = '''
                select distinct vt._Term_key, vt.abbreviation
                from VOC_AnnotType vat, VOC_Term vt
                where vat._AnnotType_key = %d
                and vat._EvidenceVocab_key = vt._Vocab_key
                ''' % CURRENT_ANNOT_TYPE
        EVIDENCE_MAP = KeyMap(evidenceCmd, '_Term_key', 'abbreviation')

        # map from qualifier term keys to their terms

        qualifierCmd = '''
                select distinct vt._Term_key, vt.term
                from VOC_AnnotType va, VOC_Term vt
                where va._AnnotType_key = %d
                and va._QualifierVocab_key = vt._Vocab_key
                ''' % CURRENT_ANNOT_TYPE
        QUALIFIER_MAP = KeyMap(qualifierCmd, '_Term_key', 'term')

        # map from user key to user login name (small data set - get them all)
        userCmd = '''select u._User_key, u.login from MGI_User u'''
        USER_MAP = KeyMap(userCmd, '_User_key', 'login') 

        # map from property term key to property name
        propertyCmd = '''
                select distinct t._Term_key, t.term
                from VOC_Term t
                where t._Vocab_key = %s
                ''' % os.environ['ANNOTPROPERTY']
        PROPERTY_MAP = KeyMap(propertyCmd, '_Term_key', 'term')

        _stamp('23:Initialized 7 key maps\n')

        return

def _initialize():
        # initialize this module by populating temporary tables; can call this
        # multiple times, as it will be a no-op if we've already initialized
        # the module

        global INITIALIZED

        if INITIALIZED:
                return

        db.useOneConnection(1)
        #_cleanupTempTables()

        _identifyExpressesComponent()
        _identifyMutationInvolves()
        _countAllelePairsPerGenotype('GXD_AllelePair')

        _buildKeepersTable()
        _keepNaturallySimpleGenotypes()
        _indexKeepersTable()

        _buildScratchPad() 

        _identifyReporterTransgenes()
        _removeConditionalGenotypes()
        _removeReporterTransgenes()

        _identifyTransactivators()
        _removeTransactivators()

        _identifyWildTypeAlleles()
        _removeWildTypeAllelesFromScratchPad()

        _collectMarkerSets()

        _handleMultipleMarkers()
        _handleMutationInvolves()
        _handleTransgenes()
        _handleDockingSites()
        _handleOtherSingles()

        _removeNullsAndGtRosa()

        # At this point, genotype_keepers has the genotype/marker pairs where
        # we can definitively identify a causative marker for a genotype's
        # MP/disease annotations.

        _getMarkerMetaData()
        _initializeKeyMaps()

        # And now, we have globals with counts of annotations for each marker
        # and an ordered list of marker keys.  (to use in grouping data when
        # pulling data from the database)  And, we have initialized our key
        # generators for the tables we will be loading data into.

        INITIALIZED = True

        return 

def _getNextMarkerBatch():
        # get the next pair of (startMarker, endMarker) that has a number of
        # annotations that will fit within MAX_ANNOTATIONS.  returns (None,
        # None) if there are no more markers.  In the event that a single
        # marker has more than the allowed number of annotations, then we will
        # return that marker only and let it be processed solo.

        global ANNOTATION_COUNTS, MARKER_KEYS, LAST_MARKER_KEY_INDEX

        maxIndex = len(MARKER_KEYS) - 1

        if LAST_MARKER_KEY_INDEX == None:
                startIndex = 0
        else:
                startIndex = LAST_MARKER_KEY_INDEX + 1

        if startIndex > maxIndex:
                return None, None

        endIndex = startIndex

        total = ANNOTATION_COUNTS[MARKER_KEYS[startIndex]]

        while (total < MAX_ANNOTATIONS) and (endIndex < maxIndex):
                i = endIndex + 1
                total = total + ANNOTATION_COUNTS[MARKER_KEYS[i]]

                if (total < MAX_ANNOTATIONS):
                        endIndex = i 

        LAST_MARKER_KEY_INDEX = endIndex

        return MARKER_KEYS[startIndex], MARKER_KEYS[endIndex]

def _makeDictionary (rows, keyField):
        # take the given list of database 'rows' and turn them into a
        # dictionary using as keys the value in each row corresponding to
        # 'keyField'.
        # Returns: { value : [ row 1, row 2, ... ] }

        out = {}
        for row in rows:
                #print(type(keyField))
                #print(type(row))
                if not row.has_key(keyField):
                        raise Error('Missing key (%s) in row: %s' % (
                                keyField, str(row)))
                key = row[keyField]

                if key not in out:
                        out[key] = []

                out[key].append(row)

        return out

def _getAnnotations (startMarker, endMarker):
        # get all rows from VOC_Annot which can be rolled up to markers between
        # the given 'startMarker' and 'endMarker', inclusive.  
        # Returns: { marker key : [ annotation rows ] }

        _stamp('25:_getAnnotations')

        cmd = '''
                select distinct k._Marker_key, a.*
                from genotype_keepers k, VOC_Annot a
                where k._Genotype_key = a._Object_key
                and a._AnnotType_key in (%d)
                and a._Term_key != %d
                and k._Marker_key >= %d
                and k._Marker_key <= %d
                ''' % ( CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS, startMarker, endMarker)

        return _makeDictionary (db.sql(cmd, 'auto'), '_Marker_key')

def _getEvidence (startMarker, endMarker):
        # get all the rows from VOC_Evidence for annotations which can be
        # rolled up to markers between the given 'startMarker' and 'endMarker',
        # inclusive.
        # Returns: { _Annot_key : [ evidence rows ] }

        _stamp('26:_getEvidence')

        cmd = '''
                select distinct k._Marker_key, e.*
                from genotype_keepers k, VOC_Annot a, VOC_Evidence e
                where k._Genotype_key = a._Object_key
                and a._AnnotType_key in (%d)
                and a._Term_key != %d
                and k._Marker_key >= %d
                and k._Marker_key <= %d
                and a._Annot_key = e._Annot_key
                ''' % (CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS, startMarker, endMarker)

        results = db.sql(cmd, 'auto')

        return _makeDictionary (results, '_Annot_key'), results

def _getEvidenceProperties (startMarker, endMarker, rawEvidence):
        # get all the properties from VOC_Evidence_Property for evidence
        # records, which are for annotations which can be rolled up to markers
        # between the given 'startMarker' and 'endMarker', inclusive.
        # Returns: { _AnnotEvidence_key : [ property rows ] }

        _stamp('27:_getEvidenceProperties')

        cmd = '''
                select distinct k._Marker_key, e._Annot_key, p.*
                from genotype_keepers k, VOC_Annot a, VOC_Evidence e, VOC_Evidence_Property p
                where k._Genotype_key = a._Object_key
                and a._AnnotType_key in (%d)
                and a._Term_key != %d
                and k._Marker_key >= %d
                and k._Marker_key <= %d
                and a._Annot_key = e._Annot_key
                and e._AnnotEvidence_key = p._AnnotEvidence_key
                order by p._AnnotEvidence_key, p.stanza, p.sequenceNum
                ''' % ( CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS, startMarker, endMarker)

        properties = db.sql(cmd, 'auto')
        _stamp('26:Retrieved properties: %d' % len(properties))

        # need to go through the evidence records and add an extra property
        # for each one, to refer back to the _Annot_key of the annotation
        # from which this one is derived.  Note:  These evidence records are
        # all for source annotations (not derived ones), so none of them
        # would already have a _SourceAnnot_key property.

        byEvidenceKey = _makeDictionary (properties, '_AnnotEvidence_key') 

        evidenceKeys = list(byEvidenceKey.keys())
        evidenceKeys.sort()

        added = 0

        # tracks (evidence key, marker key) pairs - indicates when we have a
        # source annotation property associated with that marker through the
        # specified evidence record.
        evidenceMarkerPairs = {}

        for evidenceKey in evidenceKeys:
                rows = byEvidenceKey[evidenceKey]

                # get maximum sequence number for properties tied to this
                # evidenceKey, then increment for the new property

                seqNum = max([x['sequenceNum'] for x in rows])
                if not seqNum:
                        seqNum = 1
                else:
                        seqNum = seqNum + 1

                seqRows = []
                for row in rows:
                        markerKey = row['_Marker_key']

                        # if we don't already have a source row involving this
                        # marker, then copy the current property row as a
                        # starting point, then update the copy with necessary
                        # altered values

                        pair = (evidenceKey, markerKey)
                        if pair in evidenceMarkerPairs:
                                continue

                        evidenceMarkerPairs[pair] = True
                        seqNum = seqNum + 1

                        newProperty = copy.copy(row)

                        newProperty['_PropertyTerm_key'] = SOURCE_ANNOT_KEY
                        newProperty['sequenceNum'] = seqNum
                        newProperty['value'] = newProperty['_Annot_key']

                        seqRows.append(newProperty)

                # and add the new source property
                byEvidenceKey[evidenceKey] = rows + seqRows
                added = added + len(seqRows)

        _stamp('27:add source key properties for records with existing properties: %d' % added)

        # need to handle evidence rows which had no properties previously

        ct = 0
        for row in rawEvidence:
                evidenceKey = row['_AnnotEvidence_key']
                markerKey = row['_Marker_key']

                pair = (evidenceKey, markerKey)

                # Add a source row if the annotation had no evidence at all.

                if pair not in evidenceMarkerPairs:
                        r = {
                                '_Marker_key' : markerKey,
                                '_AnnotEvidence_key' : evidenceKey,
                                '_PropertyTerm_key' : SOURCE_ANNOT_KEY,
                                'stanza' : 1,
                                'sequenceNum' : 1,
                                'value' : row['_Annot_key'],
                                '_CreatedBy_key' : row['_CreatedBy_key'],
                                '_ModifiedBy_key' : row['_ModifiedBy_key'],
                                'creation_date' : row['creation_date'],
                                'modification_date' : row['modification_date'],
                                }

                        if evidenceKey not in byEvidenceKey:
                                byEvidenceKey[evidenceKey] = [ r ]
                        else:
                                byEvidenceKey[evidenceKey].append(r)

                        evidenceMarkerPairs[pair] = True
                        ct = ct + 1

        _stamp('27:add source key properties for records with no existing properties: %d' % ct)
        _stamp('27:add source key properties in all: %d' % len(byEvidenceKey))

        return byEvidenceKey

def _getNotes (startMarker, endMarker):
        # get notes from MGI_Note for evidence records, which
        # are for annotations which can be rolled up to markers between the
        # given 'startMarker' and 'endMarker', inclusive.
        # Returns: { _AnnotEvidence_key : { note key : { record from database } } }
        # handle basic data for each note

        _stamp('28:_getNotes')

        # GENERAL_NOTE = 1008		# note type key for general notes for evidence
        # BACKGROUND_SENSITIVITY_NOTE = 1015	# note type key for background;sensitivity notes for evidence

        cmd = '''
                select distinct k._Marker_key, n.*
                from genotype_keepers k,
                        VOC_Annot a,
                        VOC_Evidence e,
                        MGI_Note n
                where k._Genotype_key = a._Object_key
                        and a._AnnotType_key in (%d)
                        and a._Term_key != %d
                        and k._Marker_key >= %d
                        and k._Marker_key <= %d
                        and a._Annot_key = e._Annot_key
                        and e._AnnotEvidence_key = n._Object_key
                        and n._NoteType_key in (1008, 1015)             -- general note/background;sensitivity note
                order by n._Object_key''' % (
                        CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS,
                        startMarker, endMarker,
                        )

        results = db.sql(cmd, 'auto')

        notes = {}		# evidence key -> notes
        noteToEvidence = {}	# note key -> evidence key

        for row in results:
                evidenceKey = row['_Object_key']
                noteKey = row['_Note_key']

                if evidenceKey in notes:
                        notes[evidenceKey][noteKey] = row
                else:
                        notes[evidenceKey] = { noteKey : row }

                noteToEvidence[noteKey] = evidenceKey

        #print(notes)
        return notes

def _splitByMarker (results):
        # take a dictionary keyed by some other field, with values being a
        # list of data rows (each containing a _Marker_key), and group those
        # rows into a dictionary with the marker keys at the top level, as:
        # 	{ marker key : { original key : [ original rows ] }
        # This is so we can easily determine which rows go with which markers.

        byMarker = {}

        originalKeys = list(results.keys())
        for key in originalKeys:
                for row in results[key]:
                        markerKey = row['_Marker_key']

                        if markerKey not in byMarker:
                                byMarker[markerKey] = { key : [ row ] }

                        elif key not in byMarker[markerKey]:
                                byMarker[markerKey][key] = [ row ]

                        else:
                                byMarker[markerKey][key].append (row)
        return byMarker

def _splitNotesByMarker (notes):
        # takes dictionary of results, as returned by _getNotes() and splits
        # them up by marker, returning a dictionary where a marker key
        # goes to a dictionary of those evidence records for that marker.

        byMarker = {}

        for evidenceKey in list(notes.keys()):
                noteKeys = list(notes[evidenceKey].keys())

                # just need to get the marker from the first note

                if noteKeys:
                        noteKey = noteKeys[0]
                        markerKey = notes[evidenceKey][noteKey]['_Marker_key']

                # Each evidence key is only for one marker, so we don't need
                # to worry about overwriting an existing evidenceKey for a
                # marker.

                if markerKey not in byMarker:
                        byMarker[markerKey] = { evidenceKey : notes[evidenceKey] }
                else:
                        byMarker[markerKey][evidenceKey] = notes[evidenceKey]
        return byMarker

def _getMarkers (startMarker, endMarker):
        # get a list of Marker objects (including annotations, evidence, and
        # properties) for all markers between (and including) the two given
        # marker keys

        _stamp('24:_getMarkers')

        # marker key -> list of annotation rows
        annotations = _getAnnotations(startMarker, endMarker)

        # marker key -> annotation key -> evidence rows 
        evidenceResults, rawEvidence =_getEvidence(startMarker, endMarker)
        evidence = _splitByMarker(evidenceResults)
        _stamp('24:returned rawEvidence rows: %d' % len(rawEvidence))

        # marker key -> evidence key -> property rows
        properties = _splitByMarker(_getEvidenceProperties(startMarker, endMarker, rawEvidence))
        _stamp('24:received properties for markers: %d' % len(properties))

        # marker key -> evidence key -> note rows
        notes = _splitNotesByMarker(_getNotes(startMarker, endMarker))
        _stamp('24:received notes for markers: %d' % len(notes))

        # Returns: { _AnnotEvidence_key : { note key : { record from database } } }

        markerKeys = list(annotations.keys())
        markerKeys.sort()

        markers = []		# list of Marker object to return

        for markerKey in markerKeys:
                marker = Marker(markerKey)
                marker.setAnnotations(annotations[markerKey])

                if markerKey in evidence:
                        marker.setEvidence(evidence[markerKey])
                if markerKey in properties:
                        marker.setProperties(properties[markerKey])
                if markerKey in notes:
                        marker.setNotes(notes[markerKey])

                markers.append(marker)

        gc.collect()
        return markers

###--- public functions ---###

def setAnnotationType (annotType):
        # set up this module to use correct annotType

        global CURRENT_ANNOT_TYPE, SOURCE_ANNOT_KEY

        if annotType in (DO_GENOTYPE, MP_GENOTYPE):
                CURRENT_ANNOT_TYPE = annotType
                if annotType == DO_GENOTYPE:
                        SOURCE_ANNOT_KEY = 13611348
                else:
                        SOURCE_ANNOT_KEY = 13576001 
        else:
                raise Error('Unknown MGI Type: %d' % annotType)
        return

def getNextMarker():
        # get the next Marker object which has rolled-up annotations.
        # Returns: a Marker object, or None if there are no more Markers

        global MARKERS_TO_DO

        if not INITIALIZED:
                _initialize()

        if not MARKERS_TO_DO:
                startMarker, endMarker = _getNextMarkerBatch()

                if startMarker == None:
                        return None

                MARKERS_TO_DO = _getMarkers (startMarker, endMarker)
        
        if not MARKERS_TO_DO:
                return None

        # pop the first marker off the list and return it

        marker = MARKERS_TO_DO[0]
        MARKERS_TO_DO = MARKERS_TO_DO[1:]

        return marker

def addTiming(s):
        # add a timing point to the profiler, identified by item 's'
        _cleanupTempTables()
        _stamp(s)
        return

