#!/usr/local/bin/python

# Name: rollup_check.py
# Purpose: to verify the rolled-up annotations in a particular database,
#	ensuring that their 'source annotation' properties are correct

import sys
import db
import time
import gc

USAGE = '''Usage: %s <server> <database>
''' % sys.argv[0]

DEBUG = False

# annotation types

OMIM_MARKER = 1016
MP_MARKER = 1015
OMIM_GENOTYPE = 1005
MP_GENOTYPE = 1002

DERIVED_TYPE = None
SOURCE_TYPE = None

CHUNK_SIZE = 200

START = time.time()

ANNOTATION = 'annotation'
EVIDENCE = 'evidence'
PROPERTY = 'property'
NOTE = 'note'

def bailout (s, showUsage = False):
	if showUsage:
		print USAGE
	print 'Error: %s' % s
	sys.exit(1)

def debug (s):
	if DEBUG:
		sys.stderr.write ('%8.2f : %s\n' % (time.time() - START, s))
	return
	
def processCommandLine():
	if len(sys.argv) != 3:
		bailout('Incorrect command-line; need two parameters.')

	db.set_sqlLogin('mgd_public', 'mgdpub', sys.argv[1], sys.argv[2])
	db.useOneConnection(1)

	try:
		db.sql('select count(1) from MGI_dbInfo', 'auto')
	except:
		bailout('Cannot query database %s..%s' % (sys.argv[1],
			sys.argv[2]))
	return

def getDerivedAnnotationKeys():
	cmd = '''select _Annot_key
		from VOC_Annot
		where _AnnotType_key = %d
		order by _Annot_key''' % DERIVED_TYPE

	results = db.sql(cmd, 'auto')

	keys = []
	for row in results:
		keys.append(row['_Annot_key'])
	debug('Got %d derived annotation keys'% len(keys))
	return keys 

def chunk(keys):
	# split keys up into sublists of 100 keys
	
	sublists = []

	lenKeys = len(keys)

	start = 0
	end = CHUNK_SIZE

	while start <= lenKeys:
		sublists.append(keys[start:end])
		start = end
		end = end + CHUNK_SIZE

	return sublists

def flatten (rows, fields):
	# flatten the rows into a single string, including the given 'fields'
	# from each

	items = []
	for row in rows:
		s = ''
		if row.has_key('property') and \
			row['property'] == '_SourceAnnot_key':
			continue

		for field in fields:
			s = s + '||' + str(row[field])
		items.append(s)

	items.sort()
	
	return '@@'.join(items)

def superset(a,b):
	for s in b.split('@@'):
		if a.find(s) < 0:
			return 0
	return 1

class Annotation:
	def __init__ (self, annotKey):
		self.annotKey = annotKey
		self.annotationRow = None
		self.evidenceRows = []
		self.propertyRows = []
		self.noteRows = []
		return

	def addRow(self, rowType, row):
		if rowType == ANNOTATION:
			if self.annotationRow:
				bailout('Two annotation rows for annotKey %d' \
					% self.annotKey)
			self.annotationRow = row
		elif rowType == EVIDENCE:
			self.evidenceRows.append(row)
		elif rowType == PROPERTY:
			self.propertyRows.append(row)
		elif rowType == NOTE:
			self.noteRows.append(row)
		return

	def getSourceKey(self):
		for row in self.propertyRows:
			if row['property'] == '_SourceAnnot_key':
				return int(row['value'])
		return None

	def matches(self, a):
		# compare data in (derived annot) self with (source annot) a

		if not self.annotationRow:
			bailout('No annotation row for %d' % self.annotKey)

		fields = [ 'term', 'qualifier' ]

		mykey = self.annotKey
		akey = a.annotKey

		if flatten([self.annotationRow], fields) \
			!= flatten([a.annotationRow], fields):
				debug('A (derived %d, source %d)' % (mykey, akey))
				debug('self: ' + flatten([self.annotationRow], fields))
				debug('a:    ' + flatten([a.annotationRow], fields))

				return 0

		fields = [ 'evidence', 'jnumID', 'inferredFrom', ]

		if not superset(flatten(self.evidenceRows, fields),
			flatten(a.evidenceRows, fields)):
				debug('B (derived %d, source %d)' % (mykey, akey))
				debug('self: ' + flatten(self.evidenceRows, fields))
				debug('a:    ' + flatten(a.evidenceRows, fields))

				return 0

		fields = [ 'property', 'stanza', 'value', ]

		if not superset(flatten(self.propertyRows, fields),
			flatten(a.propertyRows, fields)):
				debug('C (derived %d, source %d)' % (mykey, akey))
				debug("self: " + flatten(self.propertyRows, fields))
				debug("a:    " + flatten(a.propertyRows, fields))
		
				return 0

#		fields = [ 'noteType', 'private', 'note', 'sequenceNum' ]
#
#		if flatten(self.noteRows, fields).find(
#			flatten(a.noteRows, fields)) < 0:
#				debug('D (derived %d, source %d)' % (mykey, akey))
#				debug('self: ' + flatten(self.noteRows, fields))
#				debug('a:    ' + flatten(a.noteRows, fields))
#
#				return 0

		return 1

def loadAnnotations(annotKeys):
	keyString = ','.join(map(str,annotKeys))

	# basic annotation data
	cmd1 = '''select va._Annot_key, vat.name as annotation_type,
			t.term, q.term as qualifier, va.creation_date,
			va.modification_date
		from voc_annot va, voc_annottype vat, voc_term t,
			voc_term q
		where va._Annot_key in (%s)
			and va._AnnotType_key = vat._AnnotType_key
			and va._Term_key = t._Term_key
			and va._Qualifier_key = q._Term_key
		order by va._Annot_key, vat.name, t.term''' % keyString

	# evidence rows
	cmd2 = '''select ve._Annot_key, ve._AnnotEvidence_key,
			e.abbreviation as evidence, c.jnumID,
			ve.inferredFrom, cu.login as created_by,
			mu.login as modified_by, c.numericPart
		from voc_evidence ve, voc_term e,
			bib_citation_cache c, mgi_user cu, mgi_user mu
		where ve._Annot_key in (%s)
			and ve._EvidenceTerm_key = e._Term_key
			and ve._CreatedBy_key = cu._User_key
			and ve._ModifiedBy_key = mu._User_key
			and ve._Refs_key = c._Refs_key
		order by ve._Annot_key, ve._AnnotEvidence_key, e.abbreviation,
			c.numericPart''' % keyString

	# property rows
	cmd3 = '''select ve._Annot_key, vep._AnnotEvidence_key,
			p.term as property, 
			vep.stanza, vep.sequenceNum, vep.value,
			cu.login as created_by, mu.login as modified_by
		from voc_evidence ve,
			voc_evidence_property vep, voc_term p,
			mgi_user cu, mgi_user mu
		where ve._Annot_key in (%s)
			and ve._AnnotEvidence_key = vep._AnnotEvidence_key
			and vep._PropertyTerm_key = p._Term_key
			and vep._CreatedBy_key = cu._User_key
			and vep._ModifiedBy_key = mu._User_key
		order by vep._AnnotEvidence_key, vep.sequenceNum''' % keyString

	# note rows
	cmd4 = '''select ve._Annot_key, t.noteType, t.private, c.note,
			c.sequenceNum, ve._AnnotEvidence_key
		from voc_evidence ve, mgi_note n, mgi_notetype t,
			mgi_notechunk c
		where ve._Annot_key in (%s)
			and ve._AnnotEvidence_key = n._Object_key
			and t._NoteType_key = n._NoteType_key
			and t._NoteType_key in (1008, 1015)
			and t._MGIType_key = 25
			and n._Note_key = c._Note_key
		order by ve._Annot_key, ve._AnnotEvidence_key,
			t.noteType, c.sequenceNum''' % keyString

	annotations = {}

	for row in db.sql(cmd1, 'auto'):
		annotKey = row['_Annot_key']

		if annotations.has_key(annotKey):
			annot = annotations[annotKey]
		else:
			annot = Annotation(annotKey)
			annotations[annotKey] = annot

		annot.addRow(ANNOTATION, row)
	
	debug('Got %d annotations from %d-%d' % (len(annotations), annotKeys[0],
		annotKeys[-1]))

	to_do = [ (cmd2, EVIDENCE), (cmd3, PROPERTY), (cmd4, NOTE) ]

	for (cmd, dataType) in to_do:
		rows = db.sql(cmd, 'auto')
		for row in rows:
			annotKey = row['_Annot_key']

			if not annotations.has_key(annotKey):
				bailout('Unknown annot key %d (in %s)' % (
					annotKey, dataType))

			annot = annotations[annotKey]
			annot.addRow(dataType, row)

		debug('Got %d %s rows for annotations %d-%d' % (
			len(rows), dataType, annotKeys[0], annotKeys[-1]))

	return annotations

def extractSourceKeys(annotations):
	sourceKeys = {}
	for annotKey in annotations.keys():
		sourceKey = annotations[annotKey].getSourceKey()
		if sourceKey:
			sourceKeys[sourceKey] = 1
	keys = sourceKeys.keys()
	keys.sort()
	debug('Found %d source annnotation keys' % len(keys))
	return keys

def report (s):
	print 'Mismatch: %s' % s
	return

def compare (derivedAnnotations, sourceAnnotations):
	mismatches = 0

	keys = derivedAnnotations.keys()
	keys.sort()

	for derivedKey in keys:
		derivedAnnot = derivedAnnotations[derivedKey]
		sourceKey = derivedAnnot.getSourceKey()

		if not sourceKey:
			report('No source annotation for %d' % derivedKey)
			continue

		if not sourceAnnotations.has_key(sourceKey):
			report('Missing source annotation for %d' % derivedKey)

		sourceAnnot = sourceAnnotations[sourceKey]

		if not derivedAnnot.matches(sourceAnnot):
			mismatches = mismatches + 1

	return mismatches

def main():
	global DERIVED_TYPE, SOURCE_TYPE

	processCommandLine()

	pairs = [ (OMIM_MARKER, OMIM_GENOTYPE), (MP_MARKER, MP_GENOTYPE) ]

	mismatches = 0

	for (derivedAnnotType, sourceAnnotType) in pairs:
		debug('Processing %s' % derivedAnnotType)

		DERIVED_TYPE = derivedAnnotType
		SOURCE_TYPE = sourceAnnotType

		allKeys = getDerivedAnnotationKeys()
		chunks = chunk(allKeys)

		del allKeys
		gc.collect()

		for derivedKeys in chunks:
			derivedAnnotations = loadAnnotations(derivedKeys)
			sourceKeys = extractSourceKeys(derivedAnnotations)
			sourceAnnotations = loadAnnotations(sourceKeys)

			found = compare(derivedAnnotations,
				sourceAnnotations)

			if found:
				debug('--> %d mismatches' % found)

			mismatches = mismatches + found

	if mismatches:
		bailout('Failed with %d mismatches' % len(mismatches))

	print 'All records matched'
	return

if __name__ == '__main__':
	main()
