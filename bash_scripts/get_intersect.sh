#!/bin/bash

#usage: <dir of populations> <vcf.gz> <output>


#run in order of latittude for populations
for i in ${1}*.txt
do
	echo $i >> $3.order	
	(join <(bcftools query -l $2 | sort) <(sort $i)) > ${i}.intersect
	#cut has tab as default delimiter
	#info should always be identical
	
	#Ou option speeds up process when we're piping to other bcftool command
	#bcftools view -Ou -S ${i}.intersect $2 | bcftools query -f'%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | awk '{
	
	#echo ${i}
	#bcftools view -Ou --threads 32 -S ${i}.intersect $2 | bcftools +fill-tags -Ou -- -t AF | bcftools query -H -f'%INFO/AF\n' > ${i}.AF &
	#bcftools +fill-tags ${i}.temp.vcf.gz -t AF > ${i}.AF.vcf.gz 
	#bcftools query -H -f'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' ${i}.AF.vcf.gz > ${i}.AF
done	
wait
#wait for everything to finish

#write info file
#this could change depending on what information we want down the line
#added -H parameter to add header
#bcftools query -H -f'%CHROM\t%POS\t%REF\t%ALT\n' $2 > $3.INFO 

#combine temp files into tsv 
#paste $(ls $1 | grep temp | sed "s_.*_${1}/&_" | tr '\n' ' ') | column -t -s $'\t' > $3.SNPS
#paste $(ls $1 | grep temp | sed "s_.*_${1}/&_" | tr '\n' ' ') | column -t --output-separator $' ' > $3.SNPS
#BUG WITH SED CHAR 18, likely due to slashes?
#paste -d ' ' $(ls $1 | grep temp | sed "s_.*_${1}/&_") > $3.SNPS
#now just use find
#note that find isn't sorted like ls, so we have to pipe out to sort 
#make cut down versions to just the AF field
#paste -d '\t' $(ls ${1}/*txt.AF | sort | grep txt.AF) > $3

#delete temporary files
#for i in ${1}*.temp
#do
#	rm $i
#done
#for i in ${1}*.intersect
#do
#	rm $i
#done

