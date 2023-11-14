# finding adaptive variants with BayPass

Snakemake workflows and scripts for running BayPass. See [here](https://forgemia.inra.fr/mathieu.gautier/baypass_public) for the BayPass repo. The manual is the best resource for understanding how to best use BayPass, and its file formats. Reading the [publication](https://academic.oup.com/genetics/article/201/4/1555/5930067) is required for understanding the different models. 

The workflow begins with `prep_chunks.smk`. Note that this requires a VCF file and the coordinates for all the samples in a directory. This creates all the needed files before running BayPass on all the SNPs.

Once this is done, running `run_chunks.sh` iterates through all the chunks and applies `run_chunk.smk` to each one. Once BayPass has been run on each chunk, `run_chunks.sh` will also assemble a large TSV of many variables for each SNP.

These can be passed to `scripts/plot_manhattan.py` to plot all the Bayes factors. 

Example:
```
snakemake --cores 32 -s prep_chunks.smk
bash run_chunks.sh
```
 
