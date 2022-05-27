#!/bin/sh
#
#  rollupload.sh
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
# 8/14/14 jsb
#	- initial addition

cd `dirname $0`

COMMON_CONFIG=rollupload.config

USAGE="Usage: rollupload.sh"

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
#echo 'Running rollupload.py diseaseMarker' >> ${LOG_DIAG}
#ANNOTPROPERTY=`grep '^setenv ANNOTPROPERTY' ${ROLLUPLOAD}/diseaseMarker.csh.config | awk '{print $3}'`
#export ANNOTPROPERTY
#INFILE_NAME=${DISEASEMARKER_INFILE_NAME}
#export INFILE_NAME
#${PYTHON} ${ROLLUPLOAD}/bin/rollupload.py diseaseMarker >> ${LOG_DIAG} 2>&1
#STAT=$?
#checkStatus ${STAT} "${ROLLUPLOAD}/bin/rollupload.py diseaseMarker"

#
# create MP/marker input file for annotation load
#
#echo 'Running rollupload.py mpMarker' >> ${LOG_DIAG}
#ANNOTPROPERTY=`grep '^setenv ANNOTPROPERTY' ${ROLLUPLOAD}/mpMarker.csh.config | awk '{print $3}'`
#export ANNOTPROPERTY
#INFILE_NAME=${MPMARKER_INFILE_NAME}
#export INFILE_NAME
#${PYTHON} ${ROLLUPLOAD}/bin/rollupload.py mpMarker >> ${LOG_DIAG} 2>&1
#STAT=$?
#checkStatus ${STAT} "${ROLLUPLOAD}/bin/rollupload.py mpMarker"

#
# create disease/allele input file for annotation load
#
echo 'Running rollupload.py diseaseAllele' >> ${LOG_DIAG}
ANNOTPROPERTY=`grep '^setenv ANNOTPROPERTY' ${ROLLUPLOAD}/diseaseAllele.csh.config | awk '{print $3}'`
export ANNOTPROPERTY
INFILE_NAME=${DISEASEALLELE_INFILE_NAME}
export INFILE_NAME
${PYTHON} ${ROLLUPLOAD}/bin/rollupload.py diseaseAllele >> ${LOG_DIAG} 2>&1
STAT=$?
checkStatus ${STAT} "${ROLLUPLOAD}/bin/rollupload.py diseaseAllele"

#
# create MP/allele input file for annotation load
#
echo 'Running rollupload.py mpAllele' >> ${LOG_DIAG}
ANNOTPROPERTY=`grep '^setenv ANNOTPROPERTY' ${ROLLUPLOAD}/mpAllele.csh.config | awk '{print $3}'`
export ANNOTPROPERTY
INFILE_NAME=${MPALLELE_INFILE_NAME}
export INFILE_NAME
${PYTHON} ${ROLLUPLOAD}/bin/rollupload.py mpAllele >> ${LOG_DIAG} 2>&1
STAT=$?
checkStatus ${STAT} "${ROLLUPLOAD}/bin/rollupload.py mpAllele"

#
# run annotation load for disease/marker
#
#COMMON_CONFIG_CSH=${ROLLUPLOAD}/diseaseMarker.csh.config
#echo "Running disease rollupload annotation load" >> ${LOG_DIAG}
#cd ${OUTPUTDIR}
#${ANNOTLOADER_CSH} ${COMMON_CONFIG_CSH} diseaseMarker >> ${LOG_DIAG}
#STAT=$?
#checkStatus ${STAT} "${ANNOTLOADER_CSH} ${COMMON_CONFIG_CSH} diseaseMarker"

#
# run annotation load for MP/marker
#
#COMMON_CONFIG_CSH=${ROLLUPLOAD}/mpMarker.csh.config
#echo "Running MP rollupload annotation load" >> ${LOG_DIAG}
#cd ${OUTPUTDIR}
#${ANNOTLOADER_CSH} ${COMMON_CONFIG_CSH} mpMarker >> ${LOG_DIAG}
#STAT=$?
#checkStatus ${STAT} "${ANNOTLOADER_CSH} ${COMMON_CONFIG_CSH} mpMarker"

#
# run annotation load for disease/allele
#
#COMMON_CONFIG_CSH=${ROLLUPLOAD}/diseaseAllele.csh.config
#echo "Running disease rollupload annotation load" >> ${LOG_DIAG}
#cd ${OUTPUTDIR}
#${ANNOTLOADER_CSH} ${COMMON_CONFIG_CSH} diseaseAllele >> ${LOG_DIAG}
#STAT=$?
#checkStatus ${STAT} "${ANNOTLOADER_CSH} ${COMMON_CONFIG_CSH} diseaseAllele"

#
# run annotation load for MP/allele
#
#COMMON_CONFIG_CSH=${ROLLUPLOAD}/mpAllele.csh.config
#echo "Running MP rollupload annotation load" >> ${LOG_DIAG}
#cd ${OUTPUTDIR}
#${ANNOTLOADER_CSH} ${COMMON_CONFIG_CSH} mpAllele >> ${LOG_DIAG}
#STAT=$?
#checkStatus ${STAT} "${ANNOTLOADER_CSH} ${COMMON_CONFIG_CSH} mpAllele"

#
# run postload cleanup and email logs
#
shutDown

