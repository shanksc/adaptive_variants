
#fix formatting 
#need other ~20 marine species.


DATASETS=['46-Haliotis']
#number of lines per chunk 
#increase this for full dataset, it's better to have larger files if possible.
LINES_PER_CHUNK=2000

rule all:
    input:
        expand("{dataset}/snp_files/", dataset=DATASETS)

#also add a rule for creating populations s.t. they're just individuals
rule create_populations_by_latitude:
    input:
        "{dataset}/{dataset}.coords.txt"
    output:
        "{dataset}/{dataset}.done.txt"
    shell:
        """
        mkdir {wildcards.dataset}/populations
        python3 scripts/create_pops_by_latitude.py \
         --coords {wildcards.dataset}/{wildcards.dataset}.coords.txt \
         --o {wildcards.dataset}/populations/{wildcards.dataset} \
         && echo 'done' > {wildcards.dataset}/{wildcards.dataset}.done.txt 
        """
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
        "python3 scripts/random_sample.py --SNPS {wildcards.dataset}/{wildcards.dataset}.clean.SNPS --k 100000 --o {wildcards.dataset}/{wildcards.dataset}.clean.random"
'''
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
        #mkdir {wildcards.dataset}/snp_bfs/
        split --additional-suffix=.SNPS -d --lines={LINES_PER_CHUNK} {wildcards.dataset}/{wildcards.dataset}.clean.random.SNPS {wildcards.dataset}/snp_files/chunk_
        """












