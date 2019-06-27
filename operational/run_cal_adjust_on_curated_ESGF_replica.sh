#!/bin/bash

# This script will work it's way through the directories and convert all the files
ESGF_DIR="/data/CMIP/curated_ESGF_replica"
THIS_DIR=`pwd`
echo $PWD
info_file="cal_adj_info.csv"

pmip3_gcms=""
pmip3_expts="midHolocene"

#PMIP4
pmip4_gcms="IPSL-CM6A-LR"
pmip4_expts="midHolocene lig127k"

for gcm in $pmip4_gcms
do
  for expt in $pmip4_expts
  do
    echo "activity,variable,time_freq,model,experiment,ensemble,grid_label,begdate,enddate,suffix,adj_name,calendar_type,begageBP,endageBP,agestep,begyrCE,nsimyrs,source_path,adjusted_path" > $info_file    
    cd $ESGF_DIR/$gcm/$expt
    ncfiles=`ls -d *.nc`
    echo "$ncfiles"
    cd $THIS_DIR
    for ncfile in $ncfiles
    do
      if [ $ncfile != "*.nc" ]; then
        #manipulate string
        no_nc=`echo ${ncfile%.nc}`
        yr_str=${no_nc##*_}
        prior_str=${no_nc%_*}
        start_yr=`echo $yr_str | cut -c-4`
        end_yr=`echo ${yr_str##*-} | cut -c-4`
        let length=$((10#$end_yr))-$((10#$start_yr))+1
        #write names into csv file
        echo `pwd`
        echo 'PMIP4,'${prior_str//_/,},$start_yr'01',$end_yr'12,,_cal_adj,noleap,-6000,-6000,1,1000,'$length',"'$ESGF_DIR/$gcm/$expt'/","'$ESGF_DIR/$gcm/$expt'_cal_adj/"' >> $info_file 
      fi
   done
  done
done