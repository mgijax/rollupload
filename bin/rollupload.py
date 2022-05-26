
# Name: rollupload.py
#
# Inputs:
#
# Outputs:
#
# Report:
# 
# Usage:
#
# History:
#
# jsb	08/14/2014
#	- initial implementation

import os
import sys
import rollupmarkerlib

###--- globals ---###

DEBUG = False

# annotation formatted file
annotFileName = None

# annotation file pointer
annotFile = None

###--- functions ---###

def initialize():
        # Purpose: do any necessary initialization for this load
        # Returns: 0 if initialized okay, 1 if errors occured
        # Assumes: nothing
        # Effects: opens output file
        # Throws: nothing

        global annotFileName, annotFile

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

        if len(sys.argv) < 2:
                print('Need to specify annotation type on command-line')
                return 1

        annotType = sys.argv[1]
        if annotType == 'diseaseMarker':
                rollupmarkerlib.setAnnotationType(rollupmarkerlib.DO_GENOTYPE)
        elif annotType == 'mpMarker':
                rollupmarkerlib.setAnnotationType(rollupmarkerlib.MP_GENOTYPE)
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

        global annotFile

        try:
                annotFile.close()
        except:
                print('Failed to close output file properly: %s' % annotFileName)
                return 1

        return 0

###--- main program ---###

if initialize():
        sys.exit(1)

rollupmarkerlib.addTiming('Finished initialization')

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
#  10. empty (default column 2 to MGI marker IDs)
#  11. properties - optional

if DEBUG:
        annotLine = '%s\t%s\n'
else:
        annotLine = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'

marker = rollupmarkerlib.getNextMarker()
while marker:
        for annot in marker.getAnnotations():
                annotFile.write(annotLine % tuple(annot))
        marker = rollupmarkerlib.getNextMarker()

if finalize():
        sys.exit(1)

print('Outcome:  Generated %s successfully' % annotFileName)

rollupmarkerlib.addTiming('Finished writing output file')
