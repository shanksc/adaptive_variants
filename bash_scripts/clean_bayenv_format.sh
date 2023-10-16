#!/bin/bash

#usage <rows to remove (1-based)> <snps file>

awk 'NR==FNR { nums[$0]; next} !(FNR in nums)' $1 $2



