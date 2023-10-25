#!/bin/bash

#usage: INFO Beta XtX
#given an INFO file, merge it with betai output from baypass

#info=$1
#beta=(tr ' ' '\t' < ${2})
#tr ' ' '\t' < $2 > ${2}.tabbed
#sed -i "s|' '|'\t'|g" $2 > ${2}.tabbed
#echo "$2"
#echo "$3"
if [ ! -f "${2}.tabbed" ]; then 
    python3 /home/cole/scratch2/scripts/to_tabs.py $2 $2.tabbed
fi

if [ ! -f "${3}.tabbed" ]; then 
    #echo "does not exist"
    python3 /home/cole/scratch2/scripts/to_tabs.py $3 $3.tabbed
fi

#make this multiline...
#printf "CHROMOSOME\tPOSITION\tREF\tALT\tCOVARIABLE\tMRK\tM_Pearson\tSD_Pearson\tM_Spearman\tSD_Spearman\tBF(dB)\tBeta_is\tSD_Beta_is\teBPis\tMRK\tM_P\tSD_P\tM_XtX\tSD_XtX\tXtXst\tlog10(1/pval)\n"
#in place
#add chromosome/position/ref
#we may want to fix this upstream in the prep_chunks.smk file, especially if we feel the need to add more info. 
#echo $(printf "CHROMOSOME\tPOSITON\tREF\tALT\n"; cat $1) > ${1}.header
#info should just have header now 
paste -d '\t' ${1} ${2}.tabbed ${3}.tabbed



