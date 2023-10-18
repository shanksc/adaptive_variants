#this isn't ideal but preserves compatability with rest of rules
DATASET=config['dataset']
CHUNK=config['chunk']
print(DATASET)
print(CHUNK)
rule all:
    input:
        #expand("{dataset}/{dataset}_{chunk}.info.txt", dataset=[DATASETS], chunk=[CHUNKS])
        expand("{dataset}/snp_files/{chunk}.bfs", dataset=[DATASET], chunk=[CHUNK])
        #"{DATASET}/snp_bfs/{CHUNK}.bfs"
'''
rule check_config:
    output:
        "{dataset}/{dataset}_{chunk}.info.txt"
    shell:
        "echo 'read' > {wildcards.dataset}/{wildcards.dataset}_{wildcards.chunk}.info.txt"
'''
#this maybe could just use a glob? see if it works
'''
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
'''

checkpoint split_SNPS_chunk:
    input:
        #for now just use a smaller example
        #"{dataset}/{dataset}.clean.random.SNPS"
        #"{dataset}/{dataset}_annotated_pruned_0.6.vcf.gz"
        "{dataset}/snp_files/{chunk}.SNPS"
    output:
        directory("{dataset}/snp_files/{chunk}_snps/"),
        directory("{dataset}/snp_files/{chunk}_bfs/")
    shell:
        """
        mkdir {wildcards.dataset}/snp_files/{wildcards.chunk}_bfs/
        mkdir {wildcards.dataset}/snp_files/{wildcards.chunk}_snps/
        split --additional-suffix=.SNP -d -l 2 {wildcards.dataset}/snp_files/{wildcards.chunk}.SNPS {wildcards.dataset}/snp_files/{wildcards.chunk}_snps/
        #now we can remove chunk since we have directory of individual snps
        #rm {wildcards.dataset}/snp_files/{wildcards.chunk}.SNPS
        """

rule run_bayenv_on_snp:
    input:
        "{dataset}/snp_files/{chunk}_snps/{snp}.SNP"
    output:
        "{dataset}/snp_files/{chunk}_bfs/{snp}.bf"

    #compute each snp file and then remove them 
    #we would add the memoization here in some form - we could just have some hashmap locally
    #would need to handle access from each snakemake rule, could maybe be slow
    #we also need a way to increase iterations if a file contains NaNs 
    
    #CHANGE -n PARAMETER WHEN WE ADD MORE ENV VARS, should be -n 3 if we're using coldest and warmest month
    shell:
        "./bayenv -i {wildcards.dataset}/snp_files/{wildcards.chunk}_snps/{wildcards.snp}.SNP -m {wildcards.dataset}/{wildcards.dataset}.MATRIX"
        " -e {wildcards.dataset}/{wildcards.dataset}.ENVS"
        " -p $(ls -1 {wildcards.dataset}/populations | wc -l)"
        #setting an explicit seed
        " -k 10000 -n 3 -r 4536 -t"
        #have to use -o otherwise it segfaults, likely some environmental variable causing this in bayenv
        #.bf extension is added by bayenv
        " -o {wildcards.dataset}/snp_files/{wildcards.chunk}_bfs/{wildcards.snp}"
        #error handling
        """
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then   
            #this printf would need to be changed for different number of environmental variables
            #create a placeholder bf - this problem might be random seed related, in which case we might just want to re-run
            #we can't produce bfs with 0s, so this is a fine placeholder value for bayenv failing. 
            printf '{wildcards.snp}\t0\t0\t0\n' > {wildcards.dataset}/snp_files/{wildcards.chunk}_bfs/{wildcards.snp}.bf
            exit 0
        else
            exit 0
        fi
        """

def get_chunk_bfs(wildcards):
    #this checkpoint makes sure that the correct directories were made and that snps were created
    checkpoint_output = checkpoints.split_SNPS_chunk.get(**wildcards).output[0]
    snps = glob_wildcards(f"{wildcards.dataset}/snp_files/{wildcards.chunk}_snps/{{snp}}.SNP").snp
    #return expand(f"")
    return expand(f"{wildcards.dataset}/snp_files/{wildcards.chunk}_bfs/{{snp}}.bf", snp=snps)

#calculate bayes factors using bayenv for all snps 
#we can remove temporary snp_files dir here
rule run_chunk:
    input:
        #"{dataset}/snp_bfs/",
        #"{dataset}/snp_files/",
        #this returns an expansion
        #find_snp_files
        get_chunk_bfs
    output:
        "{dataset}/snp_files/{chunk}.bfs"
    shell:
        """
        cat {wildcards.dataset}/snp_files/{wildcards.chunk}_bfs/* > {wildcards.dataset}/snp_files/{wildcards.chunk}.bfs
        rm -r {wildcards.dataset}/snp_files/{wildcards.chunk}_snps/
        #delete everything but bfs - test this later
        rm -r {wildcards.dataset}/snp_files/{wildcards.chunk}_bfs/
        """