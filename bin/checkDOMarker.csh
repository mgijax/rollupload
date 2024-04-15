#!/bin/csh -f

#
# Template
#


if ( ${?MGICONFIG} == 0 ) then
        setenv MGICONFIG /usr/local/mgi/live/mgiconfig
endif

source ${MGICONFIG}/master.config.csh

cd `dirname $0`

setenv LOG $0.log
rm -rf $LOG
touch $LOG
 
date | tee -a $LOG
 
cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 | tee -a $LOG

select distinct m.symbol, aa.accid, o.commonname, ga.accid
from VOC_Annot a, MRK_Marker m, ACC_Accession aa, MGI_Organism o
, VOC_Evidence e, VOC_Evidence_Property p
, VOC_Annot mp, ACC_Accession ga

where a._annottype_key = 1015
and a._object_key = m._marker_key
and m._organism_key = o._organism_key
and a._term_key = aa._object_key
and aa._logicaldb_key = 34
and aa.preferred = 1

-- 1015 annotations
and a._annot_key = e._annot_key
and e._annotevidence_key = p._annotevidence_key
and p._propertyterm_key = 13576001

-- match 1015 annotations to 1002 annotations
and p.value::int = mp._annot_key
and mp._annottype_key = 1002
and mp._object_key = ga._object_key
and ga._mgitype_key = 12
and ga._logicaldb_key = 1
and ga.preferred = 1

order by m.symbol, aa.accid
;

EOSQL

date |tee -a $LOG

