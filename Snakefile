
#fix formatting 
#need other ~20 marine species.
DATASETS=['46-Haliotis']
#number of lines per chunk 
#increase this for full dataset, it's better to have larger files if possible.
LINES_PER_CHUNK=2000

rule all:
    input:
        #expand("{dataset}/snp_files/", dataset=DATASETS),
        expand("{dataset}/all_snps.bfs", dataset=DATASETS)

#also add a rule for creating populations s.t. they're just individuals
rule create_populations_by_latitude:
    input:
        "{dataset}/{dataset}.coords.txt"
    output:
        "{dataset}/{dataset}.done.txt"
    shell:
        "python3 scripts/create_pops_by_latitude.py"
        " --coords {wildcards.dataset}/{wildcards.dataset}.coords.txt" 
        " --o {wildcards.dataset}/populations/{wildcards.dataset}"
        " && echo 'done' > {wildcards.dataset}/{wildcards.dataset}.done.txt" 

#this might be entirely redundant if we're just using done.txt
#wait for populations to be done
checkpoint wait_for_pops:
    input:
        "{dataset}/{dataset}.done.txt"
    output:
        "{dataset}/{dataset}.num_pops.txt"
    shell:
        #wouldn't be a bad idea to remove the done.txt
        "ls -1 {wildcards.dataset}/populations/ | wc -l > {wildcards.dataset}/{wildcards.dataset}.num_pops.txt"

#we need a more stringent pruning of the vcf 
#for now just use pruned LD vcf 
rule prune_vcf:
    input:
        #prune again
        "{dataset}/{dataset}_annotated_pruned_0.6.vcf.gz"
    output:
        #"46-Haliotis/46-Haliotis_pruned.vcf.gz"
        #"{dataset}/{dataset}_pruned.vcf.gz"
        #"{dataset}/{dataset}.pruned.vcf.gz"
        "{dataset}/{dataset}_annotated_pruned_0.6.vcf.gz"
    #shell:
        #the F_MISSING filter removed few
        #"bcftools filter --threads 16 -O z --exclude 'F_MISSING>0.25' {wildcards.dataset}/{wildcards.dataset}_annotated_pruned_0.6.vcf.gz > {wildcards.dataset}/{wildcards.dataset}.pruned.vcf.gz"
        #"bcftools query --output-type z -f'%AF\n' --exclude 'INFO/AF < 0.1' 46-Haliotis/46-Haliotis_clean_snps.vcf.gz > 46-Haliotis/46-Haliotis_pruned.vcf.gz"

#may need a rule to further separate vcfs by chromosom
rule vcf_to_bayenv:
    input:
        #expand("{{dataset}}/populations/{{dataset}}.{lat}.txt", lat=LATS, allow_missing=True),
        "{dataset}/{dataset}.num_pops.txt",
        #"{dataset}/{dataset}.pruned.vcf.gz"
        #"{dataset}/{dataset}.pruned.vcf.gz"
        "{dataset}/{dataset}_annotated_pruned_0.6.vcf.gz"
    output:
        "{dataset}/{dataset}.SNPS",
        "{dataset}/{dataset}.INFO"
    shell:
        "bash bash_scripts/create_bayenv_format.sh {wildcards.dataset}/populations/ {wildcards.dataset}/{wildcards.dataset}_annotated_pruned_0.6.vcf.gz {wildcards.dataset}/{wildcards.dataset}"

#these are still present
rule rm_fixed_alleles:
    input:
        "{dataset}/{dataset}.SNPS",
        "{dataset}/{dataset}.INFO"
    output:
        "{dataset}/{dataset}.clean.SNPS",
        "{dataset}/{dataset}.clean.INFO"
    shell:
        #handle case where there are no fixed alleles, script fails if length of fixed alleles is zero otherwise
        """
        cat {wildcards.dataset}/{wildcards.dataset}.SNPS | bash bash_scripts/find_fixed_alleles.sh > {wildcards.dataset}/{wildcards.dataset}.fixed_alleles.txt
        if [ -s {wildcards.dataset}/{wildcards.dataset}.fixed_alleles.txt ]
        then
            bash bash_scripts/clean_bayenv_format.sh {wildcards.dataset}/{wildcards.dataset}.fixed_alleles.txt {wildcards.dataset}/{wildcards.dataset}.SNPS > {wildcards.dataset}/{wildcards.dataset}.clean.SNPS
            bash bash_scripts/clean_bayenv_format.sh {wildcards.dataset}/{wildcards.dataset}.fixed_alleles.txt {wildcards.dataset}/{wildcards.dataset}.INFO > {wildcards.dataset}/{wildcards.dataset}.clean.INFO
        else
            cp {wildcards.dataset}/{wildcards.dataset}.SNPS {wildcards.dataset}/{wildcards.dataset}.clean.SNPS
            cp {wildcards.dataset}/{wildcards.dataset}.INFO {wildcards.dataset}/{wildcards.dataset}.clean.INFO
        fi
        """    

#just lats for now, to be extended with other environmental variables
#adding sea surface temperatures - make sure to change -n parameter for bayenv command
rule generate_env:
    input:
        #expand("{{dataset}}/populations/{{dataset}}.{lat}.txt", lat=LATS, allow_missing=True),
        "{dataset}/{dataset}.num_pops.txt"
    output:
        "{dataset}/{dataset}.ENVS"
    shell:
        #python3 scripts/format_coords.py --coords {wildcards.dataset}/{wildcards.dataset}.coords.txt --samples $(ls {wildcards.dataset}/populations/ | grep .txt | sed 's_.*_{wildcards.dataset}/populations/&_' | tr '\n' ' ') --o {wildcards.dataset}/{wildcards.dataset}"
        "python3 scripts/create_env.py --coords {wildcards.dataset}/{wildcards.dataset}.coords.txt"
        " --samples $(ls {wildcards.dataset}/populations/ | grep .txt | sed 's_.*_{wildcards.dataset}/populations/&_' | tr '\n' ' ')"
        " --sst sst/sst.mon.ltm.1991-2020.nc"
        " --o {wildcards.dataset}/{wildcards.dataset}"

rule create_random_sample:
    input:
        "{dataset}/{dataset}.clean.SNPS"
    output:
        "{dataset}/{dataset}.clean.random.SNPS"
    shell:
        #REMEMBER TO CHANGE K
        "python3 scripts/random_sample.py --SNPS {wildcards.dataset}/{wildcards.dataset}.clean.SNPS --k 10000 --o {wildcards.dataset}/{wildcards.dataset}.clean.random"

''''
rule build_covar:
    input:
        "{dataset}/{dataset}.clean.random.SNPS"
    output:
        "{dataset}/{dataset}.MATRIX.OUT"
    #may need to adjust 
    shell:
        "./bayenv -i {wildcards.dataset}/{wildcards.dataset}.clean.random.SNPS -p $(ls -1 {wildcards.dataset}/populations | wc -l)"
        " -k 10000 > {wildcards.dataset}/{wildcards.dataset}.MATRIX.OUT"

#get last n lines to define covariance matrix 
rule final_matrix:
    input:
        "{dataset}/{dataset}.MATRIX.OUT"
    output:
        "{dataset}/{dataset}.MATRIX"
    shell:
        "sed '$d' {wildcards.dataset}/{wildcards.dataset}.MATRIX.OUT | tail -n $(ls -1 {wildcards.dataset}/populations | wc -l) > {wildcards.dataset}/{wildcards.dataset}.MATRIX"
'''

checkpoint prep_bayenv:
    input:
        "{dataset}/{dataset}.MATRIX",
        "{dataset}/{dataset}.ENVS",
        "{dataset}/{dataset}.clean.random.SNPS",
        #"{dataset}/{dataset}.clean.SNPS"
    output:
        "{dataset}/bayenv_ready.txt"
    shell:
        "echo 'ready' > {wildcards.dataset}/bayenv_ready.txt"


#split into chunks of N lines
#split snps file into chunks
#make this that this maintains order
#test with .random for now
#before we create chunks lets require that bayenv is ready
checkpoint create_SNPS_chunks:
    input:
        "{dataset}/{dataset}.clean.random.SNPS",
        "{dataset}/bayenv_ready.txt"
    output:
        directory("{dataset}/snp_files/")
    shell:
        """
        rm {wildcards.dataset}/bayenv_ready.txt
        mkdir {wildcards.dataset}/snp_files/
        split --additional-suffix=.SNPS -d --lines={LINES_PER_CHUNK} {wildcards.dataset}/{wildcards.dataset}.clean.random.SNPS {wildcards.dataset}/snp_files/chunk_
        """  
#this maybe could just use a glob? see if it works
from pathlib import Path
def get_snps(wc):
    chk_output= checkpoints.create_SNPS_chunks.get(**wc).output[0]
    out = Path(chk_output).glob("*.SNPS")
    ns = []
    for f in out:
        s = f.name.replace(".SNPS","").split("_")[1]
        ns.append(s)
    fs = expand("{{dataset}}/snp_files/chunk_{n}.SNPS", n=ns)

    return ns
#split snps file into individual snp files
checkpoint split_SNPS_chunk:
    input:
        #for now just use a smaller example
        #"{dataset}/{dataset}.clean.random.SNPS"
        "{dataset}/snp_files/chunk_{n}.SNPS"
        
    output: 
        directory("{dataset}/snp_files/chunk_{n}_snps/")
    shell: 
        """
        mkdir {wildcards.dataset}/snp_files/chunk_{wildcards.n}_snps/
        split --additional-suffix=.SNP -d -l 2 {wildcards.dataset}/snp_files/chunk_{wildcards.n}.SNPS {wildcards.dataset}/snp_files/chunk_{wildcards.n}_snps/
        #now we can remove chunk since we have directory of individual snps
        rm {wildcards.dataset}/snp_files/chunk_{wildcards.n}.SNPS
        """


rule run_bayenv_on_snp:
    input:
        "{dataset}/snp_files/chunk_{n}_snps/{snp}.SNP"
    output:
        "{dataset}/snp_files/chunk_{n}_snps/{snp}.bf"

    #compute each snp file and then remove them 
    #we would add the memoization here in some form - we could just have some hashmap locally
    #would need to handle access from each snakemake rule, could maybe be slow
    #we also need a way to increase iterations if a file contains NaNs 
    
    #CHANGE -n PARAMETER WHEN WE ADD MORE ENV VARS, should be -n 3 if we're using coldest and warmest month
    shell:
        "./bayenv -i {wildcards.dataset}/snp_files/chunk_{wildcards.n}_snps/{wildcards.snp}.SNP -m {wildcards.dataset}/{wildcards.dataset}.MATRIX"
        " -e {wildcards.dataset}/{wildcards.dataset}.ENVS"
        " -p $(ls -1 {wildcards.dataset}/populations | wc -l)"
        #setting an explicit seed
        " -k 10000 -n 3 -r 4536 -t"
        #have to use -o otherwise it segfaults, likely some environmental variable causing this in bayenv
        #.bf extension is added by bayenv
        " -o {wildcards.dataset}/snp_files/chunk_{wildcards.n}_snps/{wildcards.snp}"
        #error handling
        """
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then   
            #this printf would need to be changed for different number of environmental variables
            #create a placeholder bf - this problem might be random seed related, in which case we might just want to re-run
            #we can't produce bfs with 0s, so this is a fine placeholder value for bayenv failing. 
            printf '{wildcards.snp}\t0\t0\t0\n' > {wildcards.dataset}/snp_files/chunk_{wildcards.n}_snps/{wildcards.snp}.bf
            exit 0
        else
            exit 0
        fi
        """

def get_chunk_bfs(wildcards):
    checkpoint_output = checkpoints.split_SNPS_chunk.get(**wildcards).output[0]
    snps = glob_wildcards(f"{wildcards.dataset}/snp_files/chunk_{wildcards.n}_snps/{{snp}}.SNP").snp
    #return expand(f"")
    return expand(f"{wildcards.dataset}/snp_files/chunk_{wildcards.n}_snps/{{snp}}.bf", snp=snps)


#input will call apply_bayenv_to_chunk
rule run_by_chunk:
    input:
        #note that this expansion is limited to all the snp files within the given chunk/dataset
        get_chunk_bfs
    output:
        "{dataset}/snp_files/chunk_{n}.bfs"
    shell:
        """
        #combine bfs files into one and delete rest
        cat {wildcards.dataset}/snp_files/chunk_{wildcards.n}_snps/*.bf > {wildcards.dataset}/snp_files/chunk_{wildcards.n}.bfs 
        #now delete all the individual files
        rm -r {wildcards.dataset}/snp_files/chunk_{wildcards.n}_snps/
        """

#ensure that chunk_n_snps directories have been made before run_all_chunks
def get_all_chunks_snps(wildcards):
    ns = glob_wildcards(f"{wildcards.dataset}/snp_files/chunk_{{n}}").n
    #checkpoint_output = checkpoints.chunks
    #return expand(f"{wildcards.dataset}/snp_files/chunk_{{n}}_snps/", n=ns)
    return expand(f"{wildcards.dataset}/snp_files/chunk_{{n}}_snps", n=ns)


#don't need this anymore
'''
checkpoint prep_chunks:
    input:
        get_all_chunks_snps,
        "{dataset}/chunks_ready.txt",
    output:
        "{dataset}/chunks_ready.txt"
    shell:
        "echo 'ready' > {wildcards.dataset}/chunks_ready.txt"
'''

def get_all_chunks_bfs(wildcards):
    # checkpoint_output = checkpoints.split_SNPS_chunk.get(**wildcards)
    ns = get_snps(wildcards)
    #print(ns)    # ns = glob_wildcards(f"{wildcards.dataset}/snp_files/chunk_{{n}}_snps/").n
    # ns = Path(f"{widlcards.dataset}/snp_files/chunk_{n}_snps/")
    #checkpoint_output = checkpoints.
    #with open('chunks_seen.txt', 'w') as f: 
    #    for n in ns:
    #        f.write(f'{n}\n')
    # ns=range(10)
    return expand("{dataset}/snp_files/chunk_{n}.bfs", dataset=wildcards.dataset, n=ns)

    
#this will require everything we want as input
rule run_all_chunks:
    input:
        #this needs to be an expansion
        #"{dataset}/snp_files/chunk_{n}_snps/all.bfs"
        # "{dataset}/chunks_ready.txt",
        get_all_chunks_bfs
    output:
        #combine all chunks into single file 
        "{dataset}/all_snps.bfs"
    shell:
        #this final all_snps.bfs should be same order and length as clean.INFO and clean.SNPS
        """
        cat {wildcards.dataset}/snp_files/*.bfs > {wildcards.dataset}/all_snps.bfs
        #remove everything that we don't need 
        rm -r {wildcards.dataset}/snp_files/
        #rm {wildcards.dataset}/chunks_ready.txt
        """
