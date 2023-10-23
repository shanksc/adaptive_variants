#this isn't ideal but preserves compatability with rest of rules
DATASET=config['dataset']
CHUNK=config['chunk']
print(DATASET)
print(CHUNK)
rule all:
    input:
        expand("{dataset}_baypass/snp_files/{chunk}_summary_betai_reg.out", dataset=DATASET, chunk=CHUNK)

rule run_baypass:
    input:
        "{dataset}_baypass/snp_files/{chunk}.SNPS"
    output:
        "{dataset}_baypass/snp_files/{chunk}_summary_betai_reg.out"
    shell:
        "./g_baypass" 
        " -gfile {wildcards.dataset}_baypass/snp_files/{wildcards.chunk}.SNPS"
        " -efile {wildcards.dataset}_baypass/{wildcards.dataset}.ENVS -nthreads 64"
        " -outprefix {wildcards.dataset}/snp_files/{wildcards.chunk}"