


#in the future we probably want to utilize the config to extend our options
#DATASETS=['46-Haliotis']
DATASETS=['73-Haliotis']
#how to many lines per chunk 
LINES_PER_CHUNK=50000

rule all:
    input:
        expand("{dataset}/{dataset}.SNPS", dataset=DATASETS),
        expand("{dataset}/{dataset}.ENVS", dataset=DATASETS),
        expand("{dataset}/snp_files", dataset=DATASETS)

#add creating a txt containing samples per population ie just save the .intersect files we produce temporarly 
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

#calculate the allele freqs for each population we create. 
rule get_AF_by_population:
    input:
        "{dataset}/{dataset}.done.txt"
    output:
        "{dataset}/{dataset}.AF.tsv"
    shell:
        """
        bash bash_scripts/get_AF_for_population.sh {wildcards.dataset}/populations/ {wildcards.dataset}/{wildcards.dataset}_annotated_pruned_0.6.vcf.gz {wildcards.dataset}/{wildcards.dataset}.AF.tsv
        """

#placeholder for further filtering of vcf in future 
rule prune_vcf:
    input:
        "{dataset}/{dataset}_annotated_pruned_0.6.vcf.gz"
    output:
        "{dataset}/{dataset}_annotated_pruned_0.6.vcf.gz"
    #shell:

#convert vcf to baypass format 
rule vcf_to_baypass:
    input:
        "{dataset}/{dataset}.done.txt",
        "{dataset}/{dataset}_annotated_pruned_0.6.vcf.gz"
    output:
        "{dataset}/{dataset}.SNPS",
        "{dataset}/{dataset}.INFO"
    shell:
        #"ls {wildcards.dataset}/populations/ | grep temp| sed 's_.*_{wildcards.dataset}/populations/&_' > out.txt"
        "bash bash_scripts/create_baypass_format.sh {wildcards.dataset}/populations/ {wildcards.dataset}/{wildcards.dataset}_annotated_pruned_0.6.vcf.gz {wildcards.dataset}/{wildcards.dataset}"

#Don't think we need this for baypass, easy to re-implement though in the case that we do need it
#this just calls a script to remove fixed alleles from out SNPs file, could or maybe should be done upstream in VCFs?
'''
#these are still present
rule rm_fixed_alleles:
    input:
        "{dataset}/{dataset}.JAJLRC010000004.1.SNPS",
        "{dataset}/{dataset}.JAJLRC010000004.1.INFO"
    output:
        "{dataset}/{dataset}.JAJLRC010000004.1.clean.SNPS",
        "{dataset}/{dataset}.JAJLRC010000004.1.clean.INFO"
    shell:
        #handle case where there are no fixed alleles, script fails if length of fixed alleles is zero otherwise
        """
        cat {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.SNPS | bash bash_scripts/find_fixed_alleles.sh > {wildcards.dataset}/{wildcards.dataset}.fixed_alleles.txt
        if [ -s {wildcards.dataset}/{wildcards.dataset}.fixed_alleles.txt ]
        then
            bash bash_scripts/clean_bayenv_format.sh {wildcards.dataset}/{wildcards.dataset}.fixed_alleles.txt {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.SNPS > {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.clean.SNPS
            bash bash_scripts/clean_bayenv_format.sh {wildcards.dataset}/{wildcards.dataset}.fixed_alleles.txt {wildcards.dataset}/{wildcards.dataset}.INFO > {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.clean.INFO
        else
            cp {wildcards.dataset}/{wildcards.dataset}.SNPS {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.clean.SNPS
            cp {wildcards.dataset}/{wildcards.dataset}.INFO {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.clean.INFO
        fi
'''  

#currently this is just latitude and SST. Although the .ENVs file contains the min monthly max monthly temps in addition to latitude. 
#this would need expanded to 20+ variables in future, with n PCA components so they don't correlate. 
rule generate_env:
    input:
        #expand("{{dataset}}/populations/{{dataset}}.{lat}.txt", lat=LATS, allow_missing=True),
        "{dataset}/{dataset}.done.txt"
    output:
        "{dataset}/{dataset}.ENVS"
    shell:
        #python3 scripts/format_coords.py --coords {wildcards.dataset}/{wildcards.dataset}.coords.txt --samples $(ls {wildcards.dataset}/populations/ | grep .txt | sed 's_.*_{wildcards.dataset}/populations/&_' | tr '\n' ' ') --o {wildcards.dataset}/{wildcards.dataset}"
        "python3 scripts/create_env.py --coords {wildcards.dataset}/{wildcards.dataset}.coords.txt"
        " --samples $(find {wildcards.dataset}/populations/ | sort | grep .txt | tr '\n' ' ')"
        " --sst sst/sst.mon.ltm.1991-2020.nc"
        " --o {wildcards.dataset}/{wildcards.dataset}"

#generate a random sample of k SNPs to build covariance matrix
rule create_random_sample:
    input:
        "{dataset}/{dataset}.SNPS"
    output:
        "{dataset}/{dataset}.random.SNPS"
    shell:
        #REMEMBER TO CHANGE K
        "python3 scripts/random_sample.py --SNPS {wildcards.dataset}/{wildcards.dataset}.SNPS --k 100000 --o {wildcards.dataset}/{wildcards.dataset}.random"

#to build the covariance matrix, we can just run baypass on the random sample
rule build_covar:
    input:
        "{dataset}/{dataset}.random.SNPS",
        "{dataset}/{dataset}.ENVS"
    output:
        "{dataset}/{dataset}.random_mat_omega.out"
    #may need to adjust
    #note the number of threads as well, should probably be passed in config 
    shell: 
        "./g_baypass" 
        " -gfile {wildcards.dataset}/{wildcards.dataset}.random.SNPS"
        " -efile {wildcards.dataset}/{wildcards.dataset}.ENVS -nthreads 64"
        " -outprefix {wildcards.dataset}/{wildcards.dataset}.random"

#make sure everything is ready before running baypass 
#may be rudundant. 
checkpoint prep_baypass:
    input:
        "{dataset}/{dataset}.random_mat_omega.out",
        "{dataset}/{dataset}.ENVS",
        "{dataset}/{dataset}.SNPS",
        #"{dataset}/{dataset}.clean.SNPS"
    output:
        "{dataset}/baypass_ready.txt"
    shell:
        "echo 'ready' > {wildcards.dataset}/baypass_ready.txt"

#before we create chunks lets require that baypass is ready
checkpoint create_SNPS_chunks:
    input:
        "{dataset}/{dataset}.SNPS",
        "{dataset}/baypass_ready.txt"
    output:
        directory("{dataset}/snp_files/")
    shell:
        """
        rm {wildcards.dataset}/baypass_ready.txt
        mkdir {wildcards.dataset}/snp_files/
        split --additional-suffix=.SNPS -d --lines={LINES_PER_CHUNK} {wildcards.dataset}/{wildcards.dataset}.SNPS {wildcards.dataset}/snp_files/chunk_
        """











