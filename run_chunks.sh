#for each dataset, create chunks with prep_chunks.smk, and then for each chunk, run run_chunk.smk
#pass by using snakemake's config ie snakemake --config muscle-params="-msf
#tempting to use GNU parallel
for dataset in "46-Haliotis"
do
    echo $dataset
    #make sure that snp_files contains only the chunks
    for n in $dataset/snp_files/*.SNPS
    #for n in $dataset/snp_files/chunk_00.SNPS
    do  
        echo $n
        #declare $file_name=$(basename $n)
        eval file_name=$(basename $n)
        eval file_name_without_ext=${file_name%.*}
        snakemake -s run_chunk.smk --quiet="all" --cores 64 --scheduler "greedy" --rerun-incomplete --config chunk=$file_name_without_ext dataset=$dataset > ${dataset}/${file_name_without_ext}.out
    done 
    #combine the bfs into 1 
    cat ${dataset}/snp_files/*.bfs > ${dataset}/snps.bfs

done