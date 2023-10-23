#!/bin/bash

#usage: INFO betai.out
#given an INFO file, merge it with betai output from baypass



#info=$1
#beta=(tr ' ' '\t' < ${2})
#tr ' ' '\t' < $2 > ${2}.tabbed
#sed -i "s|' '|'\t'|g" $2 > ${2}.tabbed
if [ ! -f "${2}.tabbed" ]; then 
    python3 /home/cole/scratch2/scripts/to_tabs.py $2 $2.tabbed
fi

#printf "CHROMOSOME\tPOSITION\tREF\tALT\tCOVARIABLE\tMRK\tM_Pearson\tSD_Pearson\tM_Spearman\tSD_Spearman\tBF(dB)\tBeta_is\tSD_Beta_is\teBPis\n"
paste -d '\t' $1 ${2}.tabbed


#next filter on size
