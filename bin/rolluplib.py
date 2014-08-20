# Name: rolluplib.py
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
import db
import Profiler

###--- globals ---###

Error = 'rolluplib.Error'

MAX_ANNOTATIONS = 10000		# maximum number of annotations to cache in
				# ...memory at once

PROFILING_ON = True		# is profiling turned on currently?

MGI = 1				# logical DB key for MGI

ANNOT_EVIDENCE = 25		# MGI Type key for annotation evidence record
VOCAB_TERM = 13			# MGI Type key for vocabulary term
MARKER = 2			# MGI Type key for markers
REFERENCE = 1			# MGI Type key for references

GENERAL_NOTE = 1008		# note type key for general notes for evidence
NORMAL_NOTE = 1031		# note type key for normal notes for evidence

BACKGROUND_SENSITIVITY_NOTE = 1015	# note type key for background
					# ...sensitivity notes for evidence

CURRENT_ANNOT_TYPE = None	# either OMIM_GENOTYPE or MP_GENOTYPE
OMIM_GENOTYPE = 1005		# OMIM/Genotype annotation type
MP_GENOTYPE = 1002		# Mammalian Phenotype/Genotype annotation type

ALLELE_SUBTYPE = 1014		# Allele/Subtype annotation type
OMIM_MARKER = 1016		# new OMIM/Marker annotation type
MP_MARKER = 1015		# new Mammalian Phenotype/Marker annotation type

REPORTER = 11025589		# term key for Reporter allele attribute
TRANSGENIC = 847126		# term key for Transgenic allele type
RECOMBINASE = 11025588		# term key for Recombinase allele attribute

GT_ROSA = 37270			# marker key for Gt(ROSA)26Sor marker

INITIALIZED = False		# have we finished initializing this module?

PROFILER = Profiler.Profiler()	# used for timing of code, to aid optimization

ANNOTATION_COUNTS = {}		# marker key -> count of rolled-up annotations

MARKER_KEYS = []		# ordered list of marker keys

LAST_MARKER_KEY_INDEX = None	# index into MARKER_KEYS of last marker key
				# ...which had its details loaded

MARKERS_TO_DO = []		# list of markers loaded and waiting to be
				# ...processed

TERM_MAP = None			# KeyMap for term key -> term ID
MARKER_MAP = None		# KeyMap for marker key -> marker ID
JNUM_MAP = None			# KeyMap for refs key -> J: number
EVIDENCE_MAP = None		# KeyMap for evidence key -> evidence abbrev.
QUALIFIER_MAP = None		# KeyMap for qualifier key -> qualifier term
USER_MAP = None			# KeyMap for user key -> user
PROPERTY_MAP = None		# KeyMap for property key -> property name

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

		if self.mapping.has_key(key):
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

		# evidence key -> { note key : { note record } +
		#	{ chunks : [ note chunk rows ] } }
		self.notes = {}

		# output data to be computed by finalize() method, in columns
		# expected by annotload product 
		self.finalAnnotations = []
		return

	def _concatenateNotes (self, evidenceKey):
		# concatenate the various note chunks together from noteDict
		# for the given 'evidenceKey'

		notes = ''
		if self.notes.has_key(evidenceKey):
			for (noteKey, noteRow) in \
			    self.notes[evidenceKey].items():

				# if this row has note chunks (it should), 
				# add them to 'notes'

				if noteRow.has_key('chunks'):

					# make sure we have a space between
					# this note and the previous one (and
					# that we eliminate extra whitespace)

					if notes:
						notes = notes.strip() + ' '

					for row in noteRow['chunks']:
						notes = notes + row['note']
		return notes.replace('\n', ' ').replace('\t', ' ').strip()

	def _buildPropertiesValue (self, evidenceKey):
		# build a string to encapsulate the various properties in the
		# manner expected by the annotload

		stanzas = []

		if self.evidenceProperties.has_key(evidenceKey):
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

				clause = '%s&=&%s' % (
				    PROPERTY_MAP.get(row['_PropertyTerm_key']),
				    row['value'] )

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

			termID = TERM_MAP.get(annotRow['_Term_key'])
			markerID = MARKER_MAP.get(annotRow['_Marker_key'])
			qualifier = QUALIFIER_MAP.get(
				annotRow['_Qualifier_key'])

			if qualifier == None:
				qualifier = ''

			if not self.evidence.has_key(annotKey):
				continue

			for evidRow in self.evidence[annotKey]:
				evidKey = evidRow['_AnnotEvidence_key']

				inferredFrom = evidRow['inferredFrom']
				evidenceCode = EVIDENCE_MAP.get(
					evidRow['_EvidenceTerm_key'])
				jnumID = JNUM_MAP.get(
					evidRow['_Refs_key'])
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

				# build properties string to include any
				# properties

				properties = self._buildPropertiesValue(evidKey)

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
			s = 'Cannot set %s on finalized marker: %d' % (
				setWhat, self.markerKey)
			raise Error, s
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

	global PROFILER, PROFILING_ON

	if PROFILING_ON:
		PROFILER.stamp (s)
	return

def _countAllelePairsPerGenotype():
	# collect into a temp table (#genotype_pair_counts) the genotypes that
	# have annotations attached, along with the count of allele pairs for
	# each genotype

	cmd0 = '''select p._Genotype_key, count(1) as pair_count
		into #genotype_pair_counts
		from GXD_AllelePair p
		where exists (select 1 from VOC_Annot v
			where v._AnnotType_key in (%d)
			and v._Object_key = p._Genotype_key)
		group by p._Genotype_key''' % CURRENT_ANNOT_TYPE
	
	db.sql(cmd0, 'auto')
	_stamp('Built & filled #genotype_pair_counts table')

	# index by pair_count, as we will be using that for processing

	cmd1 = '''create index #tmp_by_count
		on #genotype_pair_counts (pair_count, _Genotype_key)'''
	db.sql(cmd1, 'auto')
	_stamp('Indexed #genotype_pair_counts')
	return

def _buildKeepersTable():
	# build the #genotype_keepers table, which will contain the genotype/
	# marker pairs where we can identify a single causative marker for the
	# genotype's annotations.  Note that the table will not be populated
	# by this method.

	cmd2 = '''create table #genotype_keepers (
		_Genotype_key int not null,
		genotype_type varchar(40) null,
		_Marker_key int null)'''
	db.sql(cmd2, 'auto')
	_stamp('Created #genotype_keepers table')
	return

def _keepSuperSimpleGenotypes():
	# add to the keepers table the super-simple genotypes, which are those
	# with:
	#   1. only one allele pair
	#   2. at least one mutant (non-wildtype) allele

	cmd3 = '''insert into #genotype_keepers
		select g._Genotype_key, "super-simple", p._Marker_key
		from #genotype_pair_counts g,
			GXD_AllelePair p
		where g.pair_count = 1
			and g._Genotype_key = p._Genotype_key
			and exists (select 1
				from GXD_AlleleGenotype gag, ALL_Allele a
				where g._Genotype_key = gag._Genotype_key
				and gag._Allele_key = a._Allele_key
				and a.isWildType = 0)'''
	db.sql(cmd3, 'auto')
	_stamp('Added super-simple genotypes to #genotype_keepers')
	return

def _indexKeepersTable():
	# add relevant indexes to the #genotype_keepers table

	cmdA = '''create index #gk_genotype
			on #genotype_keepers (_Genotype_key)'''
	cmdB = '''create index #gk_marker
			on #genotype_keepers (_Marker_key, _Genotype_key)'''

	db.sql(cmdA, 'auto')
	db.sql(cmdB, 'auto')
	_stamp('Indexed #genotype_keepers')
	return

def _identifyReporterTransgenes():
	# collect into a temp table (#reporter_transgenes) the alleles that
	# are reporter transgenes, defined as alleles with:
	#   1. an allele type (generation type) of "Transgenic"
	#   2. an allele subtype (attribute) of "Reporter"
	#   3. no other subtypes selected

	# The "distinct" should be superfluous, but we'll include it just in
	# case something flaky comes up.

	cmd4 = '''select distinct a._Object_key as _Allele_key
		into #reporter_transgenes
		from VOC_Annot a,
			ALL_Allele aa
		where a._AnnotType_key = %d
			and a._Term_key = %d
			and a._Object_key = aa._Allele_key
			and aa._Allele_Type_key = %d
			and not exists (select 1 from VOC_Annot b
				where b._AnnotType_key = %d
				and a._Object_key = b._Object_key
				and b._Term_key != %d)''' % (
		ALLELE_SUBTYPE, REPORTER, TRANSGENIC, ALLELE_SUBTYPE, REPORTER)
	db.sql(cmd4, 'auto')
	_stamp('Built & filled #reporter_transgenes table')

	# build a unique index on the allele key for the reporter transgenes

	cmd5 = '''create unique index #tmp_reportertg
		on #reporter_transgenes (_Allele_key)'''
	db.sql(cmd5, 'auto')
	_stamp('Indexed #reporter_transgenes')
	return

def _buildScratchPad():
	# For genotypes with multiple allele pairs, we need to try to apply 
	# some rules to nail things down to a single causative marker for each.
	# To begin, we'll collect a table of genotype/marker/allele data to
	# use as a scratch pad for further calculations.

	cmd6 = '''select gag._Genotype_key,
			gag._Marker_key,
			gag._Allele_key,
			g.isConditional
		into #scratchpad
		from #genotype_pair_counts c,
			GXD_AlleleGenotype gag,
			GXD_Genotype g
		where c.pair_count != 1
			and c._Genotype_key = gag._Genotype_key
			and c._Genotype_key = g._Genotype_key
			and not exists (select 1 from #genotype_keepers k
				where c._Genotype_key = k._Genotype_key)'''
	db.sql(cmd6, 'auto')
	_stamp('Built & filled #scratchpad table')

	# build an index on allele key for performance

	cmd7 = '''create index #scratch_alleles on #scratchpad (_Allele_key)'''
	db.sql(cmd7, 'auto')
	_stamp('Indexed alleles in #scratchpad')
	return

def _removeReporterTransgenes():
	# remove from the scratch pad any alleles which are reporter transgenes

	# Remove alleles which are reporter transgenes.  These appear to be
	# either hemizygous (only one allele in the pair) or homozygous (with
	# both alleles matching), so the entirety of each pair should be
	# eliminated here.
	cmd8 = '''delete from #scratchpad
		where _Allele_key in (select _Allele_key
			from #reporter_transgenes)'''
	db.sql(cmd8, 'auto')
	_stamp('Removed reporter transgenes from #scratchpad')
	return 

def _identifyWildTypeAlleles():
	# Identify which remaining alleles are wild-type alleles and include
	# them in a temp table.

	cmd9 = '''select distinct s._Allele_key
		into #wildtype_alleles
		from #scratchpad s,
			ALL_Allele a
		where s._Allele_key = a._Allele_key
			and a.isWildType = 1'''
	db.sql(cmd9, 'auto')
	_stamp('Built & filled #wildtype_alleles table')

	# build an index on allele key for performance

	cmd10 = '''create unique index #wt_alleles
		on #wildtype_alleles (_Allele_key)'''
	db.sql(cmd10, 'auto')
	_stamp('Indexed alleles in #wildtype_alleles')
	return

def _cleanupConditionalGenotypes():
	# For conditional genotypes, we need to remove recombinase alleles
	# from consideration.

	cmd11 = '''delete from #scratchpad 
		where isConditional = 1
			and _Allele_key in (select _Object_key
				from VOC_Annot 
				where _AnnotType_key = %d
				and _Term_key = %d)''' % (
					ALLELE_SUBTYPE, RECOMBINASE)
	db.sql(cmd11, 'auto')
	_stamp('Removed recombinase alleles from #scratchpad') 

	# For conditional genotypes, we also need to remove wild-type alleles.

	cmd12 = '''delete from #scratchpad
		where isConditional = 1
			and _Allele_key in (select _Allele_key
			from #wildtype_alleles)'''
	db.sql(cmd12, 'auto')
	_stamp('Removed wild-type alleles from #scratchpad')
	return

def _buildSingleMarkerTable():
	# We need to now identify those genotypes which have only one marker
	# remaining.  We'll dump these into a temp #genotype_marker_counts
	# table.

	cmd13 = '''select _Genotype_key,
			count(distinct _Marker_key) as marker_count
		into #genotype_marker_counts
		from #scratchpad
		group by _Genotype_key
		having count(distinct _Marker_key) = 1'''
	db.sql(cmd13, 'auto')
	_stamp('Built & filled #genotype_marker_counts')

	# build an index on genotype key for performance

	cmd14 = '''create unique index #gmc_genotype
		on #genotype_marker_counts (_Genotype_key)'''
	db.sql(cmd14, 'auto')
	_stamp('Indexed genotypes in #genotype_marker_counts')
	return

def _keepSimpleGenotypes():
	# Add the genotypes which now have a single marker to the set of
	# keepers.

	cmd15 = '''insert into #genotype_keepers
		select distinct s._Genotype_key,
			"simple",
			s._Marker_key
		from #genotype_marker_counts gmc,
			#scratchpad s
		where gmc._Genotype_key = s._Genotype_key'''
	db.sql(cmd15, 'auto')
	_stamp('Added pairs to #genotype_keepers')
	return

def _removeSimpleGenotypesFromScratchPad():
	# Remove the genotypes we just handled from further consideration.

	cmd16 = '''delete from #scratchpad 
		where _Genotype_key in (
			select _Genotype_key from #genotype_marker_counts)'''
	db.sql(cmd16, 'auto')
	_stamp('Removed those genotypes from #scratchpad')
	return

def _removeWildTypeAllelesFromScratchPad():
	# Remove any remaining wild-type alleles.  (We had removed them only
	# for conditional genotypes previously.)

	cmd17 = '''delete from #scratchpad 
		where _Allele_key in (select _Allele_key
			from #wildtype_alleles)'''
	db.sql(cmd17, 'auto')
	_stamp('Removed remaining wild-type alleles')
	return

def _refreshSingleMarkerTable():
	# clean out the #genotype_marker_counts table and repopulate it with
	# new data

	cmd18 = '''delete from #genotype_marker_counts'''
	db.sql(cmd18, 'auto')
	_stamp('Truncated #genotype_marker_counts')

	# re-populate #genotype_marker_counts with the genotypes that have now
	# been reduced to a single causative marker.

	cmd19 = '''insert into #genotype_marker_counts
		select _Genotype_key,
			count(distinct _Marker_key) as marker_count
		from #scratchpad
		group by _Genotype_key
		having count(distinct _Marker_key) = 1'''
	db.sql(cmd19, 'auto')
	_stamp('Refilled #genotype_marker_counts')
	return

def _removeNullsAndGtRosa():
	# finally, delete all the genotypes associated with Gt(ROSA)26Sor

	cmd20 = '''delete from #genotype_keepers
		where _Marker_key = %d''' % GT_ROSA
	db.sql(cmd20, 'auto')
	cmd21 = '''delete from #genotype_keepers
		where _Marker_key is null'''
	db.sql(cmd21, 'auto')
	_stamp('Removed genotypes for null marker and Gt(ROSA)')
	return

def _cleanupTempTables():
	# drop any temp tables that we're done with

	tables = [ '#genotype_pair_counts',
		'#reporter_transgenes',
		'#scratchpad',
		'#wildtype_alleles',
		'#genotype_marker_counts'
		]

	for table in tables:
		db.sql('drop table %s' % table, 'auto')

	_stamp('Removed unneeded temp tables')
	return

def _getMarkerMetaData():
	# populate global variables with counts of annotations for each marker
	# and an ordered list of marker keys to process.

	global ANNOTATION_COUNTS, MARKER_KEYS, LAST_MARKER_KEY_INDEX

	ANNOTATION_COUNTS = {}
	MARKER_KEYS = []
	LAST_MARKER_KEY_INDEX = None

	cmd22 = '''select k._Marker_key, count(1) as annotation_count
		from #genotype_keepers k, voc_annot a
		where k._Genotype_key = a._Object_key
		and a._AnnotType_key in (%s)
		group by k._Marker_key''' % CURRENT_ANNOT_TYPE

	results = db.sql(cmd22, 'auto')

	for row in results:
		ANNOTATION_COUNTS[row['_Marker_key']] = row['annotation_count']

	MARKER_KEYS = ANNOTATION_COUNTS.keys()
	MARKER_KEYS.sort()

	_stamp('Retrieved annotation counts for %d markers' % len(MARKER_KEYS)) 
	return

def _initializeKeyMaps():
	# initialize the mappings from various database keys to their
	# respective values

	global TERM_MAP, MARKER_MAP, JNUM_MAP, EVIDENCE_MAP
	global QUALIFIER_MAP, USER_MAP, PROPERTY_MAP, CURRENT_ANNOT_TYPE

	if not CURRENT_ANNOT_TYPE:
		raise Error, 'Need to call setAnnotationType()'

	# map from annotated term IDs to their IDs

	termCmd = '''select distinct aa._Object_key, aa.accID
		from VOC_Annot va, ACC_Accession aa
		where va._AnnotType_key in (%d)
			and va._Term_key = aa._Object_key
			and aa._MGIType_key = %d
			and aa.private = 0
			and aa.preferred = 1''' % (CURRENT_ANNOT_TYPE,
				VOCAB_TERM)
	TERM_MAP = KeyMap(termCmd, '_Object_key', 'accID')

	# map from annotated markers to their MGI IDs

	markerCmd = '''select distinct aa._Object_key, aa.accID
		from #genotype_keepers k, ACC_Accession aa
		where k._Marker_key = aa._Object_key
			and aa._MGIType_key = %d
			and aa.private = 0
			and aa.preferred = 1
			and aa._LogicalDB_key = %d''' % (MARKER, MGI)
	MARKER_MAP = KeyMap(markerCmd, '_Object_key', 'accID')

	# map from reference keys to their Jnum IDs

	jnumCmd = '''select distinct r._Refs_key, r.jnumID
		from VOC_Annot va, VOC_Evidence ve, BIB_Citation_Cache r
		where va._AnnotType_key = %d
			and va._Annot_key = ve._Annot_key
			and ve._Refs_key = r._Refs_key''' % CURRENT_ANNOT_TYPE
	JNUM_MAP = KeyMap(jnumCmd, '_Refs_key', 'jnumID')

	# map from evidence term keys to their abbreviations

	evidenceCmd = '''select distinct vt._Term_key, vt.abbreviation
		from VOC_AnnotType vat, VOC_Term vt
		where vat._AnnotType_key = %d
		and vat._EvidenceVocab_key = vt._Vocab_key''' % \
			CURRENT_ANNOT_TYPE
	EVIDENCE_MAP = KeyMap(evidenceCmd, '_Term_key', 'abbreviation')

	# map from qualifier term keys to their terms

	qualifierCmd = '''select distinct vt._Term_key, vt.term
		from VOC_AnnotType va, VOC_Term vt
		where va._AnnotType_key = %d
			and va._QualifierVocab_key = vt._Vocab_key''' % \
				CURRENT_ANNOT_TYPE
	QUALIFIER_MAP = KeyMap(qualifierCmd, '_Term_key', 'term')

	# map from user key to user login name (small data set - get them all)

	userCmd = '''select u._User_key, u.login from MGI_User u'''
	USER_MAP = KeyMap(userCmd, '_User_key', 'login') 

	# map from property term key to property name

	propertyCmd = '''select distinct t._Term_key, t.term
		from VOC_Annot va, VOC_Evidence ve, VOC_Evidence_Property p,
			VOC_Term t
		where va._AnnotType_key in (%d)
			and va._Annot_key = ve._Annot_key
			and ve._AnnotEvidence_key = p._AnnotEvidence_key
			and p._PropertyTerm_key = t._Term_key''' % \
				CURRENT_ANNOT_TYPE
	PROPERTY_MAP = KeyMap(propertyCmd, '_Term_key', 'term')

	_stamp('Initialized 7 key maps')
	return

def _initialize():
	# initialize this module by populating temporary tables; can call this
	# multiple times, as it will be a no-op if we've already initialized
	# the module

	global INITIALIZED

	if INITIALIZED:
		return

	db.useOneConnection(1)

	_countAllelePairsPerGenotype()

	_buildKeepersTable()
	_keepSuperSimpleGenotypes()
	_indexKeepersTable()

	_identifyReporterTransgenes()
	_buildScratchPad()
	_removeReporterTransgenes()

	_identifyWildTypeAlleles()
	_cleanupConditionalGenotypes()

	_buildSingleMarkerTable()
	_keepSimpleGenotypes()

	_removeSimpleGenotypesFromScratchPad()
	_removeWildTypeAllelesFromScratchPad()

	_refreshSingleMarkerTable()
	_keepSimpleGenotypes()

	_removeNullsAndGtRosa()
	_cleanupTempTables()

	# At this point, #genotype_keepers has the genotype/marker pairs where
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
		if not row.has_key(keyField):
			raise Error, 'Missing key (%s) in row: %s' % (
				keyField, str(row))
		key = row[keyField]

		if not out.has_key(key):
			out[key] = []
		out[key].append(row)

	return out

def _getAnnotations (startMarker, endMarker):
	# get all rows from VOC_Annot which can be rolled up to markers between
	# the given 'startMarker' and 'endMarker', inclusive.  
	# Returns: { marker key : [ annotation rows ] }

	cmd23 = '''select k._Marker_key, a.*
		from #genotype_keepers k, VOC_Annot a
		where k._Genotype_key = a._Object_key
			and a._AnnotType_key in (%d)
			and k._Marker_key >= %d
			and k._Marker_key <= %d''' % (
				CURRENT_ANNOT_TYPE, startMarker, endMarker)

	return _makeDictionary (db.sql(cmd23, 'auto'), '_Marker_key')

def _getEvidence (startMarker, endMarker):
	# get all the rows from VOC_Evidence for annotations which can be
	# rolled up to markers between the given 'startMarker' and 'endMarker',
	# inclusive.
	# Returns: { _Annot_key : [ evidence rows ] }

	cmd24 = '''select k._Marker_key, e.*
		from #genotype_keepers k,
			VOC_Annot a,
			VOC_Evidence e
		where k._Genotype_key = a._Object_key
			and a._AnnotType_key in (%d)
			and k._Marker_key >= %d
			and k._Marker_key <= %d
			and a._Annot_key = e._Annot_key''' % (
				CURRENT_ANNOT_TYPE, startMarker, endMarker)

	return _makeDictionary (db.sql(cmd24, 'auto'), '_Annot_key')

def _getEvidenceProperties (startMarker, endMarker):
	# get all the properties from VOC_Evidence_Property for evidence
	# records, which are for annotations which can be rolled up to markers
	# between the given 'startMarker' and 'endMarker', inclusive.
	# Returns: { _AnnotEvidence_key : [ property rows ] }

	cmd25 = '''select k._Marker_key, p.*
		from #genotype_keepers k, VOC_Annot a, VOC_Evidence e,
			VOC_Evidence_Property p
		where k._Genotype_key = a._Object_key
			and a._AnnotType_key in (%d)
			and k._Marker_key >= %d
			and k._Marker_key <= %d
			and a._Annot_key = e._Annot_key
			and e._AnnotEvidence_key = p._AnnotEvidence_key
		order by p._AnnotEvidence_key, p.stanza, p.sequenceNum''' % (
			CURRENT_ANNOT_TYPE, startMarker, endMarker)

	return _makeDictionary (db.sql(cmd25, 'auto'), '_AnnotEvidence_key') 

def _getNotes (startMarker, endMarker):
	# get notes from MGI_Note & MGI_NoteChunk for evidence records, which
	# are for annotations which can be rolled up to markers between the
	# given 'startMarker' and 'endMarker', inclusive.
	# Returns: { _AnnotEvidence_key : { note key : { record from database
	#	+ chunks : [ rows from note chunk table ] } } }

	# handle basic data for each note

	cmd26 = '''select k._Marker_key, n.*
		from #genotype_keepers k,
			VOC_Annot a,
			VOC_Evidence e,
			MGI_Note n
		where k._Genotype_key = a._Object_key
			and a._AnnotType_key in (%d)
			and k._Marker_key >= %d
			and k._Marker_key <= %d
			and a._Annot_key = e._Annot_key
			and e._AnnotEvidence_key = n._Object_key
			and n._NoteType_key in (%d, %d)
		order by e._AnnotEvidence_key''' % (
			CURRENT_ANNOT_TYPE, startMarker, endMarker,
			GENERAL_NOTE, BACKGROUND_SENSITIVITY_NOTE)

	results = db.sql(cmd26, 'auto')

	notes = {}		# evidence key -> notes
	noteToEvidence = {}	# note key -> evidence key

	for row in results:
		evidenceKey = row['_Object_key']
		noteKey = row['_Note_key']

		if notes.has_key(evidenceKey):
			notes[evidenceKey][noteKey] = row
		else:
			notes[evidenceKey] = { noteKey : row }

		noteToEvidence[noteKey] = evidenceKey

	# now get the actual note chunks and associate them with their notes

	cmd27 = '''select c.*
		from #genotype_keepers k,
			VOC_Annot a,
			VOC_Evidence e,
			MGI_Note n,
			MGI_NoteChunk c
		where k._Genotype_key = a._Object_key
			and a._AnnotType_key in (%d)
			and k._Marker_key >= %d
			and k._Marker_key <= %d
			and a._Annot_key = e._Annot_key
			and e._AnnotEvidence_key = n._Object_key
			and n._NoteType_key in (%d, %d)
			and n._Note_key = c._Note_key
		order by c._Note_key, c.sequenceNum''' % (
			CURRENT_ANNOT_TYPE, startMarker, endMarker,
			GENERAL_NOTE, BACKGROUND_SENSITIVITY_NOTE)

	results = db.sql(cmd27, 'auto')

	for row in results:
		noteKey = row['_Note_key']
		evidenceKey = noteToEvidence[noteKey]

		note = notes[evidenceKey][noteKey]
		if note.has_key('chunks'):
			note['chunks'].append(row)
		else:
			note['chunks'] = [ row ]
	return notes

def _splitByMarker (results):
	# take a dictionary keyed by some other field, with values being a
	# list of data rows (each containing a _Marker_key), and group those
	# rows into a dictionary with the marker keys at the top level, as:
	# 	{ marker key : { original key : [ original rows ] }
	# This is so we can easily determine which rows go with which markers.

	byMarker = {}

	originalKeys = results.keys()
	for key in originalKeys:
		for row in results[key]:
			markerKey = row['_Marker_key']

			if not byMarker.has_key(markerKey):
				byMarker[markerKey] = { key : [ row ] }

			elif not byMarker[markerKey].has_key(key):
				byMarker[markerKey][key] = [ row ]

			else:
				byMarker[markerKey][key].append (row)
	return byMarker

def _splitNotesByMarker (notes):
	# takes dictionary of results, as returned by _getNotes() and splits
	# them up by marker, returning a dictionary where a marker key
	# goes to a dictionary of those evidence records for that marker.

	byMarker = {}

	for evidenceKey in notes.keys():
		noteKeys = notes[evidenceKey].keys()

		# just need to get the marker from the first note

		if noteKeys:
			noteKey = noteKeys[0]
			markerKey = notes[evidenceKey][noteKey]['_Marker_key']

		# Each evidence key is only for one marker, so we don't need
		# to worry about overwriting an existing evidenceKey for a
		# marker.

		if not byMarker.has_key(markerKey):
			byMarker[markerKey] = {
				evidenceKey : notes[evidenceKey] }
		else:
			byMarker[markerKey][evidenceKey] = notes[evidenceKey]
	return byMarker

def _getMarkers (startMarker, endMarker):
	# get a list of Marker objects (including annotations, evidence, and
	# properties) for all markers between (and including) the two given
	# marker keys

	# marker key -> list of annotation rows
	annotations = _getAnnotations(startMarker, endMarker)

	# marker key -> annotation key -> evidence rows 
	evidence = _splitByMarker(_getEvidence(startMarker, endMarker))

	# marker key -> evidence key -> property rows
	properties = _splitByMarker(_getEvidenceProperties(startMarker,
		endMarker))

	# marker key -> evidence key -> note rows -> note chunk rows
	notes = _splitNotesByMarker(_getNotes(startMarker, endMarker))

	# Returns: { _AnnotEvidence_key : { note key : { record from database
	#	+ chunks : [ rows from note chunk table ] } } }

	markerKeys = annotations.keys()
	markerKeys.sort()

	markers = []		# list of Marker object to return

	for markerKey in markerKeys:
		marker = Marker(markerKey)
		marker.setAnnotations(annotations[markerKey])

		if evidence.has_key(markerKey):
			marker.setEvidence(evidence[markerKey])
		if properties.has_key(markerKey):
			marker.setProperties(properties[markerKey])
		if notes.has_key(markerKey):
			marker.setNotes(notes[markerKey])

		markers.append(marker)

	gc.collect()
	return markers

###--- public functions ---###

def setAnnotationType (annotType):
	# set up this module to use either OMIM_GENOTYPE or MP_GENOTYPE

	global CURRENT_ANNOT_TYPE

	if annotType in (OMIM_GENOTYPE, MP_GENOTYPE):
		CURRENT_ANNOT_TYPE = annotType
	else:
		raise Error, 'Unknown MGI Type: %d' % annotType
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

	_stamp(s)
	return

def dumpTimings():
	# if we collected profiling data, write it out

	if PROFILING_ON:
		print 'Profiling data:'
		PROFILER.write()
	return
