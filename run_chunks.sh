#!/bin/bash

#for each dataset, create chunks with prep_chunks.smk, and then for each chunk, run run_chunk.smk
#pass by using snakemake's config ie snakemake --config muscle-params="-msf
#tempting to use GNU parallel
for dataset in "46-Haliotis"
do
    echo $dataset
    #make sure that snp_files contains only the chunks
    for n in ${dataset}/snp_files/*.SNPS
    #for n in ${dataset}/snp_files/chunk_00.SNPS
    do  
        echo $n
        #declare $file_name=$(basename $n)
        eval file_name=$(basename $n)
        eval file_name_without_ext=${file_name%.*}
        snakemake -s run_chunk.smk --cores 1 --quiet all --rerun-incomplete --config chunk=$file_name_without_ext dataset=$dataset >> ${dataset}.out.txt
        #/usr/bin/time -o ${dataset}/snp_files/${file_name_without_ext}.time ./g_baypass -gfile $n -efile ${dataset}/46-Haliotis.ENVS -nthreads 64 -outprefix ${dataset}/snp_files/${file_name_without_ext} > ${dataset}/${file_name_without_ext}.out
        #snakemake -s run_chunk.smk --quiet "all" --cores 64 --scheduler "greedy" --rerun-incomplete --config chunk=$file_name_without_ext dataset=$dataset > ${dataset}/${file_name_without_ext}.out
    done 
    #combine the bfs into 1 
    #cat ${dataset}/snp_files/*.bfs > ${dataset}/snps.bfs
    summary_counts=$(ls -1 ${dataset}/snp_files/*summary_betai_reg.out | wc -l)
    snps_counts=$(ls -1 ${dataset}/snp_files/*SNPS | wc -l)
    if [[$summary_counts -eq $snp_counts]]
    then 
        #now we combine summary files into one to map back to .INFO files
        #remove header and then add back on later
        chunks=(${dataset}/snp_files/*.summary_betai_reg.out)
        first=$(chunks[0]) #this includes path 
        head -n 1 $first > ${dataset}/summary_betai_reg.out
        sed 1d ${dataset}/snp_files/*summary_betai_reg.out >> ${dataset}/summary_betai_reg.out
    else
        echo "Not all chunks ran. Could not concatenate output."
    fi
done
