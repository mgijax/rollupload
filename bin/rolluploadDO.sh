#!/bin/sh
#
#  rolluploadDO.sh
###########################################################################
#
#  Purpose:
# 	This script interrogates the database to find MP and disease
#	annotations which can be rolled up from genotypes to the marker level
#	and generates an output file of annotations for markers.  These new
#	annotations are then loaded into MGD using the annotation loader
#	(annotload).
#
#  Env Vars:
#
#      See the configuration file
#
#  Inputs:
#
#      - Common configuration file -
#               /usr/local/mgi/live/mgiconfig/master.config.sh
#      - rollupload load configuration file - rollupload.config
#
#  Outputs:
#
#      - An archive file
#      - Log files defined by the environment variables ${LOG_PROC},
#        ${LOG_DIAG}, ${LOG_CUR} and ${LOG_VAL}
#      - Output file: annotload
#      - see annotload outputs
#      - Records written to the database tables
#      - Exceptions written to standard error
#      - Configuration and initialization errors are written to a log file
#        for the shell script
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Fatal error occurred
#      2:  Non-fatal error occurred
#
#  Assumes:  Nothing
#
# History:
#
# 11/02/2016	lec
#	- TR12427/Disease Ontology (DO)
#

cd `dirname $0`

COMMON_CONFIG=rollupload.config

USAGE="Usage: rolluploadDO.sh"

#
# Make sure the common configuration file exists and source it.
#
if [ -f ../${COMMON_CONFIG} ]
then
    . ../${COMMON_CONFIG}
else
    echo "Missing configuration file: ${COMMON_CONFIG}"
    exit 1
fi

#
# Initialize the log file.
#
LOG=${LOG_FILE}
rm -rf ${LOG}
touch ${LOG}

#
#  Source the DLA library functions.
#

if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" | tee -a ${LOG}
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." | tee -a ${LOG}
    exit 1
fi

#####################################
#
# Main
#
#####################################

#
# createArchive including OUTPUTDIR, startLog, getConfigEnv
# sets "JOBKEY"
preload ${OUTPUTDIR}

#
# create disease/marker input file for annotation load
#
echo 'Running rollupload.py diseaseMarkerDO' >> ${LOG_DIAG}

ANNOTPROPERTY=`grep '^setenv ANNOTPROPERTY' ${ROLLUPLOAD}/diseaseMarkerDO.csh.config | awk '{print $3}'`
export ANNOTPROPERTY

INFILE_NAME=${DISEASEDO_INFILE_NAME}
export INFILE_NAME

${ROLLUPLOAD}/bin/rollupload.py diseaseMarkerDO >> ${LOG_DIAG}
STAT=$?
checkStatus ${STAT} "${ROLLUPLOAD}/bin/rollupload.py diseaseMarkerDO"

#
# run annotation load for DO-disease/marker
#

COMMON_CONFIG_CSH=${ROLLUPLOAD}/diseaseMarkerDO.csh.config
echo "Running disease rollupload annotation load" >> ${LOG_DIAG}
cd ${OUTPUTDIR}
${ANNOTLOADER_CSH} ${COMMON_CONFIG_CSH} diseaseMarker >> ${LOG_DIAG}
STAT=$?
checkStatus ${STAT} "${ANNOTLOADER_CSH} ${COMMON_CONFIG_CSH} diseaseMarkerDO"

#
# run postload cleanup and email logs
#
shutDown

