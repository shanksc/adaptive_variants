#this isn't ideal but preserves compatability with rest of rules
DATASET=config['dataset']
CHUNK=config['chunk']
print(DATASET)
print(CHUNK)
rule all:
    input:
        expand("{dataset}/snp_files/{chunk}_summary_betai_reg.out", dataset=DATASET, chunk=CHUNK)

rule run_baypass:
    input:
        "{dataset}/snp_files/{chunk}.SNPS",
        "{dataset}/{dataset}.random_mat_omega.out"
    output:
        "{dataset}/snp_files/{chunk}_summary_betai_reg.out"
    shell:
        "./g_baypass" 
        " -gfile {wildcards.dataset}/snp_files/{wildcards.chunk}.SNPS"
        " -efile {wildcards.dataset}/{wildcards.dataset}.ENVS -nthreads 64"
        " -omegafile {wildcards.dataset}/{wildcards.dataset}.random_mat_omega.out"
        " -outprefix {wildcards.dataset}/snp_files/{wildcards.chunk}"