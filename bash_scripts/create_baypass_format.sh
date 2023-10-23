#!/bin/bash

#usage: <dir of populations> <vcf.gz> <output>


#run in order of latittude for populations
for i in ${1}*
do
	echo $i >> $3.order
	#still need to remove fixed alleles	
	(join <(bcftools query -l $2 | sort) <(sort $i)) > ${i}.intersect
	#cut has tab as default delimiter
	#info should always be identical
	#Ou option speeds up process when we're piping to other bcftool command
	bcftools view -Ou -S ${i}.intersect $2 | bcftools query -f'%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | awk '{
	refs=0
	alts=0
	for (i=5; i<=NF; i++) {
		if (substr($i,1,1) ~ 0) refs+=1
		else alts += 1
		if (substr($i,3,3) ~ 0) refs+=1
		else alts += 1
	}
	printf "%d %d\n", refs, alts
}' > ${i}.temp &
done	
wait
#wait for everything to finish

#write info file
#this could change depending on what information we want down the line
bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\n' $2 > $3.INFO 

#combine temp files into tsv 
#paste $(ls $1 | grep temp | sed "s_.*_${1}/&_" | tr '\n' ' ') | column -t -s $'\t' > $3.SNPS
#paste $(ls $1 | grep temp | sed "s_.*_${1}/&_" | tr '\n' ' ') | column -t --output-separator $' ' > $3.SNPS
#BUG WITH SED CHAR 18, likely due to slashes?
#paste -d ' ' $(ls $1 | grep temp | sed "s_.*_${1}/&_") > $3.SNPS
#now just use find
#note that find isn't sorted like ls, so we have to pipe out to sort 
paste -d ' ' $(find $1 | sort | grep temp) > $3.SNPS

#delete temporary files
for i in ${1}*.temp
do
	rm $i
done
for i in ${1}*.intersect
do
	rm $i
done
