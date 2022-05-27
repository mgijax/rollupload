# Purpose: to encapsulate the rollup rules for rolling up genotype-level
#	MP and disease annotations to alleles, eliminating paths where we
#	cannot determine a single causative allele.
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
Error = 'rollupallelelib.Error'
DEBUG = True

MAX_ANNOTATIONS = 5000		# maximum number of annotations to cache in

CURRENT_ANNOT_TYPE = None	# either DO_GENOTYPE or MP_GENOTYPE
DO_GENOTYPE = 1020		# DO/Genotype annotation type
MP_GENOTYPE = 1002		# Mammalian Phenotype/Genotype annotation type

NO_PHENOTYPIC_ANALYSIS = 293594	# term key for 'no phenotypic analysis' term
SOURCE_ANNOT_KEY = None		# term key for _SourceAnnot_key property

INITIALIZED = False		# have we finished initializing this module?

ANNOTATION_COUNTS = {}		# allele key -> count of rolled-up annotations

ALLELE_KEYS = []		# ordered list of allele keys

LAST_ALLELE_KEY_INDEX = None	# index into ALLELE_KEYS of last allele key which had its details loaded

ALLELES_TO_DO = []		# list of alleles loaded and waiting to be processed

TERM_MAP = None			# KeyMap for term key -> term ID
ALLELE_MAP = None		# KeyMap for allele key -> allele ID
JNUM_MAP = None			# KeyMap for refs key -> J: number
EVIDENCE_MAP = None		# KeyMap for evidence key -> evidence abbrev.
QUALIFIER_MAP = None		# KeyMap for qualifier key -> qualifier term
USER_MAP = None			# KeyMap for user key -> user
PROPERTY_MAP = None		# KeyMap for property key -> property name

# no testing
testSQL = ""

testSQL = '''
and exists (select 1 from ACC_Accession testg 
where gag._genotype_key = testg._object_key and testg._mgitype_key = 12 
and testg.accid in (
'MGI:3776488',
'MGI:2667804',
'MGI:2173405',
'MGI:3799369',
'MGI:3723214',
'MGI:2665793',
'MGI:3814544',
'MGI:2678410',
'MGI:5439284',
'MGI:3836962',
'MGI:3510920',
'MGI:5897214',
'MGI:5558945',
'MGI:5882410',
'MGI:3613531',
'MGI:3837691',
'MGI:3776520',
'MGI:3698007',
'MGI:3814907',
'MGI:5690044',
'MGI:6450805'
)
)
'''

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

class Allele:
        # Is: a single allele from the database
        # Has: a allele key and sets of annotations, evidence, evidence
        #	properties, and evidence notes that can be rolled up to that
        #	allele.
        # Does: converts primary keys and annotation types so the data can
        #	be extracted from the object for loading into the database
        #	with no primary key conflicts
        # Notes: Once you call a get*() method to extract data from this
        #	object, you can no longer call any set*() methods.  This is
        #	due to the need to prepare the data (as noted above).

        def __init__ (self, alleleKey):
                # constructor; initializes object for allele with given key

                self.finalized = False

                # input data to be populated from outside this object:

                self.alleleKey = alleleKey

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
                # The finalize() method takes this allele object from its
                # original state (old primary keys, old annotation types) and
                # converts them to rows appropriate to be loaded as new
                # records by the annotation loader.  This can be called
                # multiple times, as it will be a no-op if the allele is
                # already finalized.

                if self.finalized:
                        return

                # The input file for the annotation loader allows up to eleven
                # fields per line.  We will generate rows in that format for
                # self.finalAnnotations, specifically for our data set:
                #    1. vocab term ID
                #    2. allele ID
                #    3. J: num
                #    4. evidence code abbreviation
                #    5. inferred from
                #    6. qualifier
                #    7. username
                #    8. empty -- use the current date by default
                #    9. notes -- can be empty
                #   10. empty -- defaults to MGI IDs for alleles
                #   11. properties -- can be empty

                self.finalAnnotations = []

                for annotRow in self.annotations:
                        annotKey = annotRow['_Annot_key']

                        if annotRow['_Term_key'] == NO_PHENOTYPIC_ANALYSIS:
                                # skip annotations to this term
                                continue

                        termID = TERM_MAP.get(annotRow['_Term_key'])
                        alleleID = ALLELE_MAP.get(annotRow['_Allele_key'])
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
                                                alleleID
                                                ]
                                else:
                                        row = [
                                                termID,
                                                alleleID,
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
                        s = 'Cannot set %s on finalized allele: %d' % (setWhat, self.alleleKey)
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
                _stampResults('genotype_keepers', 'order by gaccid, aaccid')

def _deleteKeeper():
        if DEBUG:
                _stamp('after delete from genotype_keepers')
                _stampResults('genotype_keepers', 'order by gaccid, aaccid')

def _deleteScratchpad():
        if DEBUG:
                _stamp('after delete from scratchpad')
                _stampResults('scratchpad', 'order by gaccid, aaccid')

        # re-create genotype_pair_counts from scratchpad
        _countAllelePairsPerGenotype('scratchpad')

def _getCount(table):
        results = db.sql('select count(1) as get_count from %s' % table, 'auto')
        if not results:
                return 0
        return results[0]['get_count']

def _countAllelePairsPerGenotype(table):
        # collect into a temp table (genotype_pair_counts) the genotypes that
        # have annotations attached, along with the count of allele pairs for each genotype

        _stamp('1:_countAllelePairsPerGenotype')

        db.sql('drop table if exists genotype_pair_counts', None);

        if (table == 'GXD_AlleleGenotype'):
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
        _stamp('1:genotype_pair_counts : %d\n' % _getCount('genotype_pair_counts'))

        return

def _buildKeepersTable():
        # build the genotype_keepers table, which will contain the genotype/
        # allele pairs where we can identify a single causative allele for the
        # genotype's annotations.  Note that the table will not be populated
        # by this method.

        _stamp('3:_buildKeepersTable')

        cmd = '''
                create temp table genotype_keepers (
                        _Genotype_key int not null,
                        genotype_type text null,
                        _Allele_key int null,
                        symbol text null,
                        gaccid text not null,
                        aaccid text not null
                )
                '''
        db.sql(cmd, None)
        _stamp('\n')

        return

def _keepNaturallySimpleGenotypes():
        # add to the keepers table the naturally simple genotypes:
        #
        # Rule #1
        # One Allele Genotype. Naturally simple genotypes are included in the roll-up, 
        # so we associate these genotypes with their corresponding alleles
        # 1. Have only one allele
        # 4. Are not conditional genotypes
        # 5. Ignore wild-type alleles
        #

        _stamp('4:_keepNaturallySimpleGenotypes/genotype_keepers/rule #1 : one allele genotype')
        _stamp('4.1:only one allele')
        _stamp('4.2:isConditional = false')
        _stamp('4.3:wildtype = false')

        cmd = '''
                insert into genotype_keepers
                select distinct g._Genotype_key,
                        'rule #1 : one allele genotype',
                        p._Allele_key, a.symbol, a1.accid, a2.accid
                from genotype_pair_counts g,
                        GXD_AlleleGenotype p,
                        GXD_Genotype gg,
                        ALL_Allele a,
                        ACC_Accession a1,
                        ACC_Accession a2
                where g.pair_count = 1
                        and g._Genotype_key = p._Genotype_key
                        and p._Genotype_key = gg._Genotype_key
                        -- 4.2:isConditional = false
                        and gg.isConditional = 0        
                        and p._Allele_key = a._Allele_key

                        and g._Genotype_key = a1._Object_key
                        and a1._MGIType_key = 12
                        and a1._logicaldb_key = 1
                        and a1.prefixpart = 'MGI:'
                        and p._Allele_key = a2._Object_key
                        and a2._MGIType_key = 11
                        and a2._logicaldb_key = 1
                        and a2.prefixpart = 'MGI:'
                        and a2.preferred = 1

                        -- 4.3:wildtype = false')
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
        db.sql('create index gk_allele on genotype_keepers (_Allele_key, _Genotype_key)', None)
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
        # some rules to nail things down to a single causative allele for each.
        # To begin, we'll collect a table of genotype/allele/allele data to
        # use as a scratch pad for further calculations.

        _stamp('6:_buildScratchPad/scratchpad')
        _stamp('6:exists in genotype_pair_counts')
        _stamp('6:does not exist in genotype_keepers')

        cmd = '''
                select distinct gag._Genotype_key,
                        gag._Marker_key,
                        gag._Allele_key,
                        g.isConditional,
                        a.symbol,
                        a1.accid as gaccid,
                        a2.accid as aaccid
                into temp table scratchpad
                from genotype_pair_counts c,
                        GXD_AlleleGenotype gag,
                        GXD_Genotype g,
                        ALL_Allele a,
                        ACC_Accession a1,
                        ACC_Accession a2
                where not exists (select 1 from genotype_keepers k 
                        where c._Genotype_key = k._Genotype_key)
                        and c._Genotype_key = gag._Genotype_key
                        and c._Genotype_key = g._Genotype_key
                        and gag._Allele_key = a._Allele_key
                        and gag._Genotype_key = a1._Object_key
                        and a1._MGIType_key = 12
                        and a1._logicaldb_key = 1
                        and a1.prefixpart = 'MGI:'
                        and gag._Allele_key = a2._Object_key
                        and a2._MGIType_key = 11
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
                'genotype_pair_counts',
                'genotype_keepers',
                'reporter_transgenes',
                'transactivators',
                'scratchpad',
                'wildtype_alleles',
                'trad',
                'trad_ct'
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

def _collectAlleleSets():
        # builds several temp tables with alleles tied to genotypes by three
        # routes:  traditional allele/allele pairings, expressed component
        # relationships, and mutation involves relationships.  Note that the
        # temp tables will only contain data for genotypes (and only through
        # alleles) which still exist in scratchpad.  Also note that these
        # only include genotypes which have annotations of the current
        # annotation type.  Still another note -- we exclude genotypes where
        # the only annotation is to 'no phenotypic analysis'.

        # The temp tables we will build include:
        # 1. trad - genotype-to-allele via the traditional allele-allele route
        # 2. trad_ct - genotype to count of its alleles in trad

        _stamp('14:_collectAlleleSets')

        # command to build table 1 (distinct genotype/allele pairs)

        tradCmd = '''
                select distinct s._Genotype_key, s._Allele_key
                into temp table trad
                from scratchpad s, ALL_Allele a
                where s._Allele_key = a._Allele_key
                '''

        # commands to index tables 1-3

        templateB = 'create index %s on %s (_Genotype_key)'

        tradB = templateB % ('tradIndex1', 'trad')

        # commands to build tables 4-6 (counts of distinct alleles/genotype)

        templateC = '''select _Genotype_key, count(1) as allele_count into temp table %s from %s group by _Genotype_key'''

        tradC = templateC % ('trad_ct', 'trad')

        # commands to index tables 4-6

        tradD = templateB % ('tradIndex2', 'trad_ct')
        
        templateE = 'create index %s on %s (allele_count)'

        tradE = templateE % ('tradIndex3', 'trad_ct')

        # build and index each table, and report on completion for profiling

        for (name, tbl1, tbl2, c1, c2, c3, c4, c5) in [
                ('traditional', 'trad', 'trad_ct', tradCmd, tradB, tradC, tradD, tradE)
                ]:

                db.sql(c1, None)
                db.sql(c2, None)
                db.sql(c3, None)
                db.sql(c4, None)
                db.sql(c5, None)

                _stamp('14:built table of %s genotype/allele pairs rows: %d' % (name, _getCount(tbl1)) )
                _stamp('14:built table of %s genotype/allele counts rows: %d' % (name, _getCount(tbl2)) )

                if DEBUG:
                        results = db.sql('''
                                select ac.accid, a.symbol, t.* 
                                from %s t, ALL_Allele a, ACC_Accession ac
                                where t._Allele_key = a._Allele_key
                                and a._Allele_key = ac._Object_key
                                and ac._MGIType_key = 11
                                and ac._logicaldb_key = 1
                                and ac.prefixpart = 'MGI:'
                                and ac.preferred = 1
                                order by ac.accid
                                ''' % (tbl1), 'auto')
                        for r in results:
                                _stamp(r)
         
                _stamp('\n')

        return

def _getAlleleMetaData():
        # populate global variables with counts of annotations for each allele
        # and an ordered list of allele keys to process.

        global ANNOTATION_COUNTS, ALLELE_KEYS, LAST_ALLELE_KEY_INDEX

        _stamp('22:_getAlleleMetaData')

        ANNOTATION_COUNTS = {}
        ALLELE_KEYS = []
        LAST_ALLELE_KEY_INDEX = None

        cmd = '''
                select k._Allele_key, count(1) as annotation_count
                from genotype_keepers k, VOC_Annot a
                where k._Genotype_key = a._Object_key
                and a._AnnotType_key in (%s)
                and a._Term_key != %d
                group by k._Allele_key
                ''' % (CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS)

        results = db.sql(cmd, 'auto')
        for row in results:
                ANNOTATION_COUNTS[row['_Allele_key']] = row['annotation_count']

        if DEBUG:
                for r in results:
                        _stamp(r)
                _stamp('\n')

        ALLELE_KEYS = list(ANNOTATION_COUNTS.keys())
        ALLELE_KEYS.sort()

        _stamp('22:Retrieved annotation counts for %d alleles\n' % len(ALLELE_KEYS)) 

        return

def _initializeKeyMaps():
        # initialize the mappings from various database keys to their
        # respective values

        global TERM_MAP, ALLELE_MAP, JNUM_MAP, EVIDENCE_MAP
        global QUALIFIER_MAP, USER_MAP, PROPERTY_MAP, CURRENT_ANNOT_TYPE

        _stamp('23:_initializeKeyMaps')

        if not CURRENT_ANNOT_TYPE:
                raise Error('Need to call setAnnotationType()')

        # map from annotated term IDs to their IDs
        if DEBUG:
                accID = "ac.accID || ':' || t.term as accID"
        else:
                accID = "ac.accID as accID"

        termCmd = '''
                select distinct ac._Object_key, %s
                from VOC_Annot va, ACC_Accession ac, VOC_Term t
                where va._AnnotType_key in (%d)
                and va._Term_key = ac._Object_key
                and ac._MGIType_key = 13
                and ac.private = 0
                and ac.preferred = 1
                and va._Term_key = t._Term_key
                ''' % (accID, CURRENT_ANNOT_TYPE)
        TERM_MAP = KeyMap(termCmd, '_Object_key', 'accID')

        # map from annotated alleles to their MGI IDs
        if DEBUG:
                accID = "ac.accID || ':' || a.symbol as accID"
        else:
                accID = "ac.accID as accID"

        alleleCmd = '''
                select distinct ac._Object_key, %s
                from genotype_keepers k, ACC_Accession ac, ALL_Allele a
                where k._Allele_key = ac._Object_key
                        and ac._MGIType_key = 11
                        and ac.private = 0
                        and ac.preferred = 1
                        and ac._LogicalDB_key = 1
                        and k._Allele_key = a._Allele_key
                ''' % (accID)
        ALLELE_MAP = KeyMap(alleleCmd, '_Object_key', 'accID')

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

        _countAllelePairsPerGenotype('GXD_AlleleGenotype')

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

        _collectAlleleSets()

        # At this point, genotype_keepers has the genotype/allele pairs where
        # we can definitively identify a causative allele for a genotype's
        # MP/disease annotations.

        _getAlleleMetaData()
        _initializeKeyMaps()

        # And now, we have globals with counts of annotations for each allele
        # and an ordered list of allele keys.  (to use in grouping data when
        # pulling data from the database)  And, we have initialized our key
        # generators for the tables we will be loading data into.

        INITIALIZED = True

        return 

def _getNextAlleleBatch():
        # get the next pair of (startAllele, endAllele) that has a number of
        # annotations that will fit within MAX_ANNOTATIONS.  returns (None,
        # None) if there are no more alleles.  In the event that a single
        # allele has more than the allowed number of annotations, then we will
        # return that allele only and let it be processed solo.

        global ANNOTATION_COUNTS, ALLELE_KEYS, LAST_ALLELE_KEY_INDEX

        maxIndex = len(ALLELE_KEYS) - 1

        if LAST_ALLELE_KEY_INDEX == None:
                startIndex = 0
        else:
                startIndex = LAST_ALLELE_KEY_INDEX + 1

        if startIndex > maxIndex:
                return None, None

        endIndex = startIndex

        total = ANNOTATION_COUNTS[ALLELE_KEYS[startIndex]]

        while (total < MAX_ANNOTATIONS) and (endIndex < maxIndex):
                i = endIndex + 1
                total = total + ANNOTATION_COUNTS[ALLELE_KEYS[i]]

                if (total < MAX_ANNOTATIONS):
                        endIndex = i 

        LAST_ALLELE_KEY_INDEX = endIndex

        return ALLELE_KEYS[startIndex], ALLELE_KEYS[endIndex]

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

def _getAnnotations (startAllele, endAllele):
        # get all rows from VOC_Annot which can be rolled up to alleles between
        # the given 'startAllele' and 'endAllele', inclusive.  
        # Returns: { allele key : [ annotation rows ] }

        _stamp('25:_getAnnotations')

        cmd = '''
                select distinct k._Allele_key, a.*
                from genotype_keepers k, VOC_Annot a
                where k._Genotype_key = a._Object_key
                and a._AnnotType_key in (%d)
                and a._Term_key != %d
                and k._Allele_key >= %d
                and k._Allele_key <= %d
                ''' % ( CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS, startAllele, endAllele)

        return _makeDictionary (db.sql(cmd, 'auto'), '_Allele_key')

def _getEvidence (startAllele, endAllele):
        # get all the rows from VOC_Evidence for annotations which can be
        # rolled up to alleles between the given 'startAllele' and 'endAllele',
        # inclusive.
        # Returns: { _Annot_key : [ evidence rows ] }

        _stamp('26:_getEvidence')

        cmd = '''
                select distinct k._Allele_key, e.*
                from genotype_keepers k, VOC_Annot a, VOC_Evidence e
                where k._Genotype_key = a._Object_key
                and a._AnnotType_key in (%d)
                and a._Term_key != %d
                and k._Allele_key >= %d
                and k._Allele_key <= %d
                and a._Annot_key = e._Annot_key
                ''' % (CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS, startAllele, endAllele)

        results = db.sql(cmd, 'auto')

        return _makeDictionary (results, '_Annot_key'), results

def _getEvidenceProperties (startAllele, endAllele, rawEvidence):
        # get all the properties from VOC_Evidence_Property for evidence
        # records, which are for annotations which can be rolled up to alleles
        # between the given 'startAllele' and 'endAllele', inclusive.
        # Returns: { _AnnotEvidence_key : [ property rows ] }

        _stamp('27:_getEvidenceProperties')

        cmd = '''
                select distinct k._Allele_key, e._Annot_key, p.*
                from genotype_keepers k, VOC_Annot a, VOC_Evidence e, VOC_Evidence_Property p
                where k._Genotype_key = a._Object_key
                and a._AnnotType_key in (%d)
                and a._Term_key != %d
                and k._Allele_key >= %d
                and k._Allele_key <= %d
                and a._Annot_key = e._Annot_key
                and e._AnnotEvidence_key = p._AnnotEvidence_key
                order by p._AnnotEvidence_key, p.stanza, p.sequenceNum
                ''' % ( CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS, startAllele, endAllele)

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

        # tracks (evidence key, allele key) pairs - indicates when we have a
        # source annotation property associated with that allele through the
        # specified evidence record.
        evidenceAllelePairs = {}

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
                        alleleKey = row['_Allele_key']

                        # if we don't already have a source row involving this
                        # allele, then copy the current property row as a
                        # starting point, then update the copy with necessary
                        # altered values

                        pair = (evidenceKey, alleleKey)
                        if pair in evidenceAllelePairs:
                                continue

                        evidenceAllelePairs[pair] = True
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
                alleleKey = row['_Allele_key']

                pair = (evidenceKey, alleleKey)

                # Add a source row if the annotation had no evidence at all.

                if pair not in evidenceAllelePairs:
                        r = {
                                '_Allele_key' : alleleKey,
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

                        evidenceAllelePairs[pair] = True
                        ct = ct + 1

        _stamp('27:add source key properties for records with no existing properties: %d' % ct)
        _stamp('27:add source key properties in all: %d' % len(byEvidenceKey))

        return byEvidenceKey

def _getNotes (startAllele, endAllele):
        # get notes from MGI_Note for evidence records, which
        # are for annotations which can be rolled up to alleles between the
        # given 'startAllele' and 'endAllele', inclusive.
        # Returns: { _AnnotEvidence_key : { note key : { record from database } } }
        # handle basic data for each note

        _stamp('28:_getNotes')

        # GENERAL_NOTE = 1008		# note type key for general notes for evidence
        # BACKGROUND_SENSITIVITY_NOTE = 1015	# note type key for background;sensitivity notes for evidence

        cmd = '''
                select distinct k._Allele_key, n.*
                from genotype_keepers k,
                        VOC_Annot a,
                        VOC_Evidence e,
                        MGI_Note n
                where k._Genotype_key = a._Object_key
                        and a._AnnotType_key in (%d)
                        and a._Term_key != %d
                        and k._Allele_key >= %d
                        and k._Allele_key <= %d
                        and a._Annot_key = e._Annot_key
                        and e._AnnotEvidence_key = n._Object_key
                        and n._NoteType_key in (1008, 1015)             -- general note/background;sensitivity note
                order by n._Object_key''' % (
                        CURRENT_ANNOT_TYPE, NO_PHENOTYPIC_ANALYSIS,
                        startAllele, endAllele,
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

def _splitByAllele (results):
        # take a dictionary keyed by some other field, with values being a
        # list of data rows (each containing a _Allele_key), and group those
        # rows into a dictionary with the allele keys at the top level, as:
        # 	{ allele key : { original key : [ original rows ] }
        # This is so we can easily determine which rows go with which alleles.

        byAllele = {}

        originalKeys = list(results.keys())
        for key in originalKeys:
                for row in results[key]:
                        alleleKey = row['_Allele_key']

                        if alleleKey not in byAllele:
                                byAllele[alleleKey] = { key : [ row ] }

                        elif key not in byAllele[alleleKey]:
                                byAllele[alleleKey][key] = [ row ]

                        else:
                                byAllele[alleleKey][key].append (row)
        return byAllele

def _splitNotesByAllele (notes):
        # takes dictionary of results, as returned by _getNotes() and splits
        # them up by allele, returning a dictionary where a allele key
        # goes to a dictionary of those evidence records for that allele.

        byAllele = {}

        for evidenceKey in list(notes.keys()):
                noteKeys = list(notes[evidenceKey].keys())

                # just need to get the allele from the first note

                if noteKeys:
                        noteKey = noteKeys[0]
                        alleleKey = notes[evidenceKey][noteKey]['_Allele_key']

                # Each evidence key is only for one allele, so we don't need
                # to worry about overwriting an existing evidenceKey for a
                # allele.

                if alleleKey not in byAllele:
                        byAllele[alleleKey] = { evidenceKey : notes[evidenceKey] }
                else:
                        byAllele[alleleKey][evidenceKey] = notes[evidenceKey]
        return byAllele

def _getAlleles (startAllele, endAllele):
        # get a list of Allele objects (including annotations, evidence, and
        # properties) for all alleles between (and including) the two given
        # allele keys

        _stamp('24:_getAlleles')

        # allele key -> list of annotation rows
        annotations = _getAnnotations(startAllele, endAllele)

        # allele key -> annotation key -> evidence rows 
        evidenceResults, rawEvidence =_getEvidence(startAllele, endAllele)
        evidence = _splitByAllele(evidenceResults)
        _stamp('24:returned rawEvidence rows: %d' % len(rawEvidence))

        # allele key -> evidence key -> property rows
        properties = _splitByAllele(_getEvidenceProperties(startAllele, endAllele, rawEvidence))
        _stamp('24:received properties for alleles: %d' % len(properties))

        # allele key -> evidence key -> note rows
        notes = _splitNotesByAllele(_getNotes(startAllele, endAllele))
        _stamp('24:received notes for alleles: %d' % len(notes))

        # Returns: { _AnnotEvidence_key : { note key : { record from database } } }

        alleleKeys = list(annotations.keys())
        alleleKeys.sort()

        alleles = []		# list of Allele object to return

        for alleleKey in alleleKeys:
                allele = Allele(alleleKey)
                allele.setAnnotations(annotations[alleleKey])

                if alleleKey in evidence:
                        allele.setEvidence(evidence[alleleKey])
                if alleleKey in properties:
                        allele.setProperties(properties[alleleKey])
                if alleleKey in notes:
                        allele.setNotes(notes[alleleKey])

                alleles.append(allele)

        gc.collect()
        return alleles

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

def getNextAllele():
        # get the next Allele object which has rolled-up annotations.
        # Returns: a Allele object, or None if there are no more Alleles

        global ALLELES_TO_DO

        if not INITIALIZED:
                _initialize()

        if not ALLELES_TO_DO:
                startAllele, endAllele = _getNextAlleleBatch()

                if startAllele == None:
                        return None

                ALLELES_TO_DO = _getAlleles (startAllele, endAllele)
        
        if not ALLELES_TO_DO:
                return None

        # pop the first allele off the list and return it

        allele = ALLELES_TO_DO[0]
        ALLELES_TO_DO = ALLELES_TO_DO[1:]

        return allele

def addTiming(s):
        # add a timing point to the profiler, identified by item 's'
        _cleanupTempTables()
        _stamp(s)
        return

