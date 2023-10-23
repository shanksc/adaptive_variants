#!/bin/bash

#remove fixed alleles

#check for fixed alleles and print row number if fixed
#this solution isn't very good since we're reading lines into arr
awk '{
	fst=0
	snd=0
	for(i=1;i<=NF;i+=2) {
		fst+=$i
		snd+=$(i+1)
	}
	if (fst == 0 || snd = 0) {
		print NR
	}
}
'
