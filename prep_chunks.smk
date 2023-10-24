
#fix formatting 
#need other ~20 marine species.

#in the future we probably want to utilize the config to extend our options
DATASETS=['46-Haliotis']
#with baypass we are building a covariance matrix for each chunk, which is fine as long as they are similar?
LINES_PER_CHUNK=50000

rule all:
    input:
        expand("{dataset}/{dataset}.SNPS", dataset=DATASETS),
        expand("{dataset}/{dataset}.ENVS", dataset=DATASETS),
        expand("{dataset}/snp_files", dataset=DATASETS)

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


#filtering on chrom for now
#we need a more stringent pruning of the vcf 
#for now just use pruned LD vcf 
rule prune_vcf:
    input:
        #prune again
        "{dataset}/{dataset}_annotated_pruned_0.6.vcf.gz"
    output:
        #"46-Haliotis/46-Haliotis_pruned.vcf.gz"
        #"{dataset}/{dataset}_pruned.vcf.gz"
        #"{dataset}/{dataset}.JAJLRC010000004.1.vcf.gz"
        "{dataset}/{dataset}_annotated_pruned_0.6.vcf.gz"
    #shell:
        #the F_MISSING filter removed few
        #"""
        #bcftools index {wildcards.dataset}/{wildcards.dataset}_annotated_pruned_0.6.vcf.gz -t -o {wildcards.dataset}/{wildcards.dataset}_annotated_pruned_0.6.vcf.gz.tbi --threads 32
        #bcftools view {wildcards.dataset}/{wildcards.dataset}_annotated_pruned_0.6.vcf.gz --regions JAJLRC010000004.1 > {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.vcf.gz --threads 32
        #"""
        #"bcftools filter --threads 32 -O z --exclude '%CHROM!=JAJLRC010000004.1' {wildcards.dataset}/{wildcards.dataset}_annotated_pruned_0.6.vcf.gz > {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.vcf.gz"
        #"bcftools query --output-type z -f'%AF\n' --exclude 'INFO/AF < 0.1' 46-Haliotis/46-Haliotis_clean_snps.vcf.gz > 46-Haliotis/46-Haliotis_pruned.vcf.gz"

#may need a rule to further separate vcfs by chromosom
rule vcf_to_baypass:
    input:
        "{dataset}/{dataset}.num_pops.txt",
        "{dataset}/{dataset}_annotated_pruned_0.6.vcf.gz"
    output:
        "{dataset}/{dataset}.SNPS",
        "{dataset}/{dataset}.INFO"
    shell:
        #"ls {wildcards.dataset}/populations/ | grep temp| sed 's_.*_{wildcards.dataset}/populations/&_' > out.txt"
        "bash bash_scripts/create_baypass_format.sh {wildcards.dataset}/populations/ {wildcards.dataset}/{wildcards.dataset}_annotated_pruned_0.6.vcf.gz {wildcards.dataset}/{wildcards.dataset}"

#Don't think we need this for baypass, easy to re-implement though in the case that we do need it 
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

#this is in an identical format as before
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
        " --samples $(find {wildcards.dataset}/populations/ | sort | grep .txt | tr '\n' ' ')"
        " --sst sst/sst.mon.ltm.1991-2020.nc"
        " --o {wildcards.dataset}/{wildcards.dataset}"

#we don't actually need to explicity construct the covariance matrix separately unless we want to do AUX model 
rule create_random_sample:
    input:
        "{dataset}/{dataset}.SNPS"
    output:
        "{dataset}/{dataset}.random.SNPS"
    shell:
        #REMEMBER TO CHANGE K
        "python3 scripts/random_sample.py --SNPS {wildcards.dataset}/{wildcards.dataset}.SNPS --k 100000 --o {wildcards.dataset}/{wildcards.dataset}.random"

#to build the covariance matrix, we can just run 
rule build_covar:
    input:
        "{dataset}/{dataset}.random.SNPS",
        "{dataset}/{dataset}.ENVS"
    output:
        "{dataset}/{dataset}.random_mat_omega.out"
    #may need to adjust 
    shell: 
        "./g_baypass" 
        " -gfile {wildcards.dataset}/{wildcards.dataset}.random.SNPS"
        " -efile {wildcards.dataset}/{wildcards.dataset}.ENVS -nthreads 64"
        " -outprefix {wildcards.dataset}/{wildcards.dataset}.random"
        
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

#test with .random for now
#before we create chunks lets require that bayenv is ready
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











