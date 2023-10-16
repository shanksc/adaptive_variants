
#fix formatting 
#need other ~20 marine species.
DATASETS=['46-Haliotis']

rule all:
    #temporary
    input:
        #expand("{dataset}/{dataset}.clean.random.SNPS", dataset=DATASETS),
        #expand("{dataset}/{dataset}.MATRIX.OUT", dataset=DATASETS),
        #expand("{dataset}/{dataset}.ENVS", dataset=DATASETS)
        #expand("{dataset}/{dataset}.MATRIX", dataset=DATASETS)
        #expand("{dataset}/snp_files", dataset=DATASETS)
        expand("{dataset}/snp_bfs.all", dataset=DATASETS)
        #expand("{dataset}/snp_files", dataset=DATASETS, glob)

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

#CHANGE EVERY RULE BELOW THIS TO USE ANNOTATED PRUNED SNPS
#split snps file into individual snp files
rule split_SNPS:
    input:
        #for now just use a smaller example
        #"{dataset}/{dataset}.clean.random.SNPS"
        "{dataset}/{dataset}_annotated_pruned_0.6.vcf.gz"
    output:
        directory("{dataset}/snp_files/")
    shell:
        """
        mkdir {wildcards.dataset}/snp_files/
        split -d -l 2 {wildcards.dataset}/{wildcards.dataset}_annotated_pruned_0.6.vcf.gz {wildcards.dataset}/snp_files/
        """

#memoize allele counts for each SNP s.t. that we're not just recomputing what we've already computed 
#we need to split up SNPS into individal SNP files, then run separately. Before running we check if we've already computed
#the same SNP in terms of the counts. We need this feature eventually anyways, even if its not very effective for LD
#filtered vcfs
rule prep_bayenv:
    input:
        "{dataset}/snp_files/",
        "{dataset}/{dataset}.MATRIX",
        "{dataset}/{dataset}.ENVS"
    output:
        directory("{dataset}/snp_bfs/")
    shell:
        "mkdir {wildcards.dataset}/snp_bfs"


#rule for individual SNP files
rule run_bayenv:
    input:
        "{dataset}/snp_files/{n}",
        "{dataset}/snp_bfs/",
        "{dataset}/{dataset}.MATRIX",
        "{dataset}/{dataset}.ENVS"
    output:
        "{dataset}/snp_bfs/{n}.bf"
    
    #compute each snp file and then remove them 
    #we would add the memoization here in some form - we could just have some hashmap locally
    #would need to handle access from each snakemake rule, could maybe be slow
    #we also need a way to increase iterations if a file contains NaNs 
    
    #CHANGE -n PARAMETER WHEN WE ADD MORE ENV VARS, should be -n 3 if we're using coldest and warmest month
    shell:
        "./bayenv -i {wildcards.dataset}/snp_files/{wildcards.n} -m {wildcards.dataset}/{wildcards.dataset}.MATRIX"
        " -e {wildcards.dataset}/{wildcards.dataset}.ENVS"
        " -p $(ls -1 {wildcards.dataset}/populations | wc -l)"
        " -k 10000 -n 3 -t"
        #have to use -o otherwise it segfaults, likely some environmental variable
        " -o {wildcards.dataset}/snp_bfs/{wildcards.n}"

def find_snp_files(wildcards):
    #ns = glob_wildcards("{wildcards.dataset}/snp_files/{n}")
    ns = glob_wildcards(f"{wildcards.dataset}/snp_files/{{n}}").n
    #with open(f'{wildcards.dataset}/snp_files_names.txt', 'w') as f:
    #    for n in ns:
    #        f.write(n+'\n')
    return expand(f"{wildcards.dataset}/snp_bfs/{{n}}.bf", n=ns)

#calculate bayes factors using bayenv for all snps 
#we can remove temporary snp_files dir here
rule calc_bfs:
    input:
        "{dataset}/snp_bfs/",
        "{dataset}/snp_files/",
        #this returns an expansion
        find_snp_files
    output:
        "{dataset}/snp_bfs.all"
    #remove this probably
    benchmark:
        "{dataset}/calc_bfs.benchmark.txt"
    shell:
        #rm -r {wildcards.dataset}/snp_files
        #"echo 'snps_bfs done' > {wildcards.dataset}/snp_bfs.done.txt"
        "cat {wildcards.dataset}/snp_bfs/* > {wildcards.dataset}/snp_bfs.all"
        "rm -r {wildcards.dataset}/snp_files/"
        #"rm -r {wildcards.dataset}"


