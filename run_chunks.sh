#!/bin/bash

#for each dataset, create chunks with prep_chunks.smk, and then for each chunk, run run_chunk.smk
#pass by using snakemake's config ie snakemake --config muscle-params="-msf"
#in the future if we're on a cluster we could submit to a cluster N jobs at a time and wait, where each job is a 50k or so chunk.
for dataset in "73-Haliotis"
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

    if [ $(ls -1 ${dataset}/snp_files/*summary_betai_reg.out | wc -l) == $(ls -1 ${dataset}/snp_files/*SNPS | wc -l) ]; then
        #now we combine summary files into one to map back to .INFO files
        #remove header and then add back on later
        #maybe also put this in separate script. 
        chunks=(${dataset}/snp_files/*summary_betai_reg.out)
        first=${chunks[0]} #this includes path 
        head -n 1 $first > ${dataset}/summary_betai_reg.out
        sed 1d ${dataset}/snp_files/*summary_betai_reg.out >> ${dataset}/summary_betai_reg.out
        #discard other environmental covariates for now CHANGE THIS IN FUTURE
        head -n 1 $first > ${dataset}/summary_betai_reg_latitude.out
        #here we select covariable 1, which corresponds to latitude in ENV file. 
        cat ${dataset}/summary_betai_reg.out | awk -v dataset="$dataset" '$1 == 1{print >> (dataset"/summary_betai_reg_latitude.out")}'
        #add other scripts here to create final TSV containing INFO, XtX, and betai_reg output
        bash bash_scripts/combine_xtx.sh ${dataset}/snp_files/ ${dataset}/xtxs.out
        #add allele freqs for each population
        #this will fail if that rule isn't ran previously. Currently prep_chunks doesn't have that specificed in it's all rule. 
        bash bash_scripts/add_info_to_betai.sh ${dataset}/${dataset}.INFO ${dataset}/summary_betai_reg_latitude.out ${dataset}/xtxs.out  ${dataset}/${dataset}.AF.tsv > ${dataset}/latitude.tsv
    else
        echo "Not all chunks ran. Could not concatenate output."
    fi
done
