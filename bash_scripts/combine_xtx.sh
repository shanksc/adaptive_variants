#!/bin/bash

#make to to include last forward slash
#ex dataset/snp_files/
#usage <snp_files_dir_path> <out>

#merge all xtx files
echo $1
#declare -a chunks
xtxs=(${1}*summary_pi_xtx.out)
echo ${#xtxs[@]}
first=${xtxs[0]}

head -n 1 $first > $2
for xtx in "${xtxs[@]}"
do
    sed 1d $xtx >> $2
done


