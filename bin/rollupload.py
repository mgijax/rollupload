#
# Name: rollupload.py
#
# see the README file for testing instructions
#
# History:
#
# jsb	08/14/2014
#	- initial implementation

import os
import sys
import rollupmarkerlib
import rollupallelelib

###--- globals ---###

DEBUG = False

# annotation formatted file
annotFileName = None
annotFileName2 = None

# annotation file pointer
annotFile = None
annotFile2 = None

annotType = ""

###--- functions ---###

def initialize():
        # Purpose: do any necessary initialization for this load
        # Returns: 0 if initialized okay, 1 if errors occured
        # Assumes: nothing
        # Effects: opens output file
        # Throws: nothing

        global annotFileName, annotFile, annotType
        global annotFileName2, annotFile2

        try:
                annotFileName = os.environ['INFILE_NAME']
        except:
                print('INFILE_NAME not defined in environment')
                return 1

        try:
                annotFile = open(annotFileName, 'w')
        except:
                print('Cannot open output file: %s' % annotFileName)
                return 1

        # used for mpMarkerNonMouse
        if "INFILE_NAME2" in os.environ:
                try:
                        annotFileName2 = os.environ['INFILE_NAME2']
                except:
                        print('INFILE_NAME2 not defined in environment')
                        return 1
        
                try:
                        annotFile2 = open(annotFileName2, 'w')
                except:
                        print('Cannot open output file: %s' % annotFileName2)
                        return 1

        if len(sys.argv) < 2:
                print('Need to specify annotation type on command-line')
                return 1

        annotType = sys.argv[1]
        if annotType == 'diseaseMarker':
                rollupmarkerlib.setAnnotationType(rollupmarkerlib.DO_GENOTYPE)
        elif annotType == 'mpMarker':
                rollupmarkerlib.setAnnotationType(rollupmarkerlib.MP_GENOTYPE)
        elif annotType == 'diseaseAllele':
                rollupallelelib.setAnnotationType(rollupallelelib.DO_GENOTYPE)
        elif annotType == 'mpAllele':
                rollupallelelib.setAnnotationType(rollupallelelib.MP_GENOTYPE)
        else:
                print('Unknown annotation type: %s' % annotType)
                return 1

        return 0

def finalize():
        # Purpose: do any necessary finalization for this load
        # Returns: 0 if finalized okay, 1 if errors occurred
        # Assumes: nothing
        # Effects: closes output file
        # Throws: nothing

        global annotFile, annotFile2

        try:
                annotFile.close()
        except:
                print('Failed to close output file properly: %s' % annotFileName)
                return 1

        if "INFILE_NAME2" in os.environ:
                try:
                        annotFile2.close()
                except:
                        print('Failed to close output file properly: %s' % annotFileName2)
                        return 1

        return 0

###--- main program ---###

if initialize():
        sys.exit(1)

if annotType in ('diseaseMarker', 'mpMarker'):
        rollupmarkerlib.addTiming('Finished initialization')
elif annotType in ('diseaseAllele', 'mpAllele'):
        rollupallelelib.addTiming('Finished initialization')

# template for one line for the annnotation loader:
#   1. vocab term ID
#   2. marker ID
#   3. J: number
#   4. evidence code abbreviation
#   5. inferred from
#   6. qualifier
#   7. username
#   8. empty (use the current date by default)
#   9. notes - optional
#  10. 'MGI' (mouse) or 'Entrez Gene' (non-mouse)
#  11. properties - optional

if DEBUG:
        annotLine = '%s\t%s\n'
else:
        annotLine = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'

# split diseaseMarker -> diseaseMarker (INFILE_NAME), diseaseMarkerNonMouse (INFILE_NAME2)
# split mpMarker -> mpMarker (INFILE_NAME), mpMarkerNonMouse (INFILE_NAME2)
if annotType in ('diseaseMarker', 'mpMarker'):
        marker = rollupmarkerlib.getNextMarker()
        while marker:
                for annot in marker.getAnnotations():
                        # skip; no entrez gene id
                        if annot[1] == None:
                                continue
                        elif DEBUG == True:
                                annotFile.write(annotLine % tuple(annot))
                        # mouse
                        elif annot[9] == 'MGI':
                                annotFile.write(annotLine % tuple(annot))
                        # non-mouse
                        else:
                                annotFile2.write(annotLine % tuple(annot))
                marker = rollupmarkerlib.getNextMarker()

elif annotType in ('diseaseAllele', 'mpAllele'):
        allele = rollupallelelib.getNextAllele()
        while allele:
                for annot in allele.getAnnotations():
                        annotFile.write(annotLine % tuple(annot))
                allele = rollupallelelib.getNextAllele()
        
if finalize():
        sys.exit(1)

print('Outcome:  Generated %s successfully' % annotFileName)

if annotType in ('diseaseMarker', 'mpMarker'):
        rollupmarkerlib.addTiming('Finished writing output file')
elif annotType in ('diseaseAllele', 'mpAllele'):
        rollupallelelib.addTiming('Finished writing output file')

