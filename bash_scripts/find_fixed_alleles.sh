#!/bin/bash

#remove fixed alleles

#check for fixed alleles and print row number if fixed
#this solution isn't very good since we're reading lines into arr
awk '{
	for(i=1;i<=NF;i++) {sums[NR]+=$i}
} END {
	for (l=1;l<=NR;l+=2) {
		fst = sums[l]
		snd = sums[l+1]
		if (fst == 0 || snd = 0) {
			print l 
			print l+1
		}	       	
	}
}
'
