#!/bin/bash

#given a summary_output_betai_reg.out tsv from stdin, output as many tsvs as there are covariables
tsv=$1
out=$2
cat $tsv | awk '$1 == 1{print > ("betai_latitude.tsv")}'