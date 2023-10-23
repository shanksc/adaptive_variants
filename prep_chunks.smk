
#fix formatting 
#need other ~20 marine species.


DATASETS=['46-Haliotis']
#with baypass we are building a covariance matrix for each chunk, which is fine as long as they are similar?
LINES_PER_CHUNK=50000

rule all:
    input:
        expand("{dataset}_baypass/{dataset}.SNPS", dataset=DATASETS),
        expand("{dataset}_baypass/{dataset}.ENVS", dataset=DATASETS),
        expand("{dataset}_baypass/snp_files", dataset=DATASETS)

#also add a rule for creating populations s.t. they're just individuals
rule create_populations_by_latitude:
    input:
        "{dataset}_baypass/{dataset}.coords.txt"
    output:
        "{dataset}_baypass/{dataset}.done.txt"
    shell:
        """
        mkdir {wildcards.dataset}_baypass/populations
        python3 scripts/create_pops_by_latitude.py \
         --coords {wildcards.dataset}_baypass/{wildcards.dataset}.coords.txt \
         --o {wildcards.dataset}_baypass/populations/{wildcards.dataset} \
         && echo 'done' > {wildcards.dataset}_baypass/{wildcards.dataset}.done.txt 
        """
#this might be entirely redundant if we're just using done.txt
#wait for populations to be done
checkpoint wait_for_pops:
    input:
        "{dataset}_baypass/{dataset}.done.txt"
    output:
        "{dataset}_baypass/{dataset}.num_pops.txt"
    shell:
        #wouldn't be a bad idea to remove the done.txt
        "ls -1 {wildcards.dataset}_baypass/populations/ | wc -l > {wildcards.dataset}_baypass/{wildcards.dataset}.num_pops.txt"


#filtering on chrom for now
#we need a more stringent pruning of the vcf 
#for now just use pruned LD vcf 
rule prune_vcf:
    input:
        #prune again
        "{dataset}_baypass/{dataset}_annotated_pruned_0.6.vcf.gz"
    output:
        #"46-Haliotis/46-Haliotis_pruned.vcf.gz"
        #"{dataset}/{dataset}_pruned.vcf.gz"
        #"{dataset}/{dataset}.JAJLRC010000004.1.vcf.gz"
        "{dataset}_baypass/{dataset}_annotated_pruned_0.6.vcf.gz"
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
        "{dataset}_baypass/{dataset}.num_pops.txt",
        "{dataset}_baypass/{dataset}_annotated_pruned_0.6.vcf.gz"
    output:
        "{dataset}_baypass/{dataset}.SNPS",
        "{dataset}_baypass/{dataset}.INFO"
    shell:
        #"ls {wildcards.dataset}_baypass/populations/ | grep temp| sed 's_.*_{wildcards.dataset}/populations/&_' > out.txt"
        "bash bash_scripts/create_baypass_format.sh {wildcards.dataset}_baypass/populations/ {wildcards.dataset}_baypass/{wildcards.dataset}_annotated_pruned_0.6.vcf.gz {wildcards.dataset}_baypass/{wildcards.dataset}"

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
        "{dataset}_baypass/{dataset}.num_pops.txt"
    output:
        "{dataset}_baypass/{dataset}.ENVS"
    shell:
        #python3 scripts/format_coords.py --coords {wildcards.dataset}/{wildcards.dataset}.coords.txt --samples $(ls {wildcards.dataset}/populations/ | grep .txt | sed 's_.*_{wildcards.dataset}/populations/&_' | tr '\n' ' ') --o {wildcards.dataset}/{wildcards.dataset}"
        "python3 scripts/create_env.py --coords {wildcards.dataset}_baypass/{wildcards.dataset}.coords.txt"
        " --samples $(find {wildcards.dataset}_baypass/populations/ | sort | grep .txt | tr '\n' ' ')"
        " --sst sst/sst.mon.ltm.1991-2020.nc"
        " --o {wildcards.dataset}_baypass/{wildcards.dataset}"

#we don't actually need to explicity construct the covariance matrix separately unless we want to do AUX model 
'''
rule create_random_sample:
    input:
        "{dataset}/{dataset}.JAJLRC010000004.1.clean.SNPS"
    output:
        "{dataset}/{dataset}.JAJLRC010000004.1.clean.random.SNPS"
    shell:
        #REMEMBER TO CHANGE K
        "python3 scripts/random_sample.py --SNPS {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.clean.SNPS --k  --o {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.clean.random"

rule build_covar:
    input:
        "{dataset}/{dataset}.JAJLRC010000004.1.clean.random.SNPS"
    output:
        "{dataset}/{dataset}.JAJLRC010000004.1.MATRIX.OUT"
    #may need to adjust 
    shell:
        "./bayenv -i {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.clean.random.SNPS -p $(ls -1 {wildcards.dataset}/populations | wc -l)"
        " -k 10000 > {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.MATRIX.OUT"

#get last n lines to define covariance matrix 
rule final_matrix:
    input:
        "{dataset}/{dataset}.JAJLRC010000004.1.MATRIX.OUT"
    output:
        "{dataset}/{dataset}.JAJLRC010000004.1.MATRIX"
    shell:
        "sed '$d' {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.MATRIX.OUT | tail -n $(ls -1 {wildcards.dataset}/populations | wc -l) > {wildcards.dataset}/{wildcards.dataset}.JAJLRC010000004.1.MATRIX"
'''

checkpoint prep_baypass:
    input:
        #"{dataset}/{dataset}.JAJLRC010000004.1.MATRIX",
        "{dataset}_baypass/{dataset}.ENVS",
        "{dataset}_baypass/{dataset}.SNPS",
        #"{dataset}/{dataset}.clean.SNPS"
    output:
        "{dataset}_baypass/bayenv_ready.txt"
    shell:
        "echo 'ready' > {wildcards.dataset}_baypass/bayenv_ready.txt"

#test with .random for now
#before we create chunks lets require that bayenv is ready
checkpoint create_SNPS_chunks:
    input:
        "{dataset}_baypass/{dataset}.SNPS",
        "{dataset}_baypass/bayenv_ready.txt"
    output:
        directory("{dataset}_baypass/snp_files/")
    shell:
        """
        rm {wildcards.dataset}_baypass/bayenv_ready.txt
        mkdir {wildcards.dataset}_baypass/snp_files/
        split --additional-suffix=.SNPS -d --lines={LINES_PER_CHUNK} {wildcards.dataset}_baypass/{wildcards.dataset}.SNPS {wildcards.dataset}_baypass/snp_files/chunk_
        """











