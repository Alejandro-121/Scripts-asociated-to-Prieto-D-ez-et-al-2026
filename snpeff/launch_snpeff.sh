#!/bin/bash

# load the conda env
# conda activate my_snpeff

mkdir -p out_vcf

parallel bgzip {} ::: *vcf
 
# exit if conda fails to load
for file in *.gz
do
	# get prefix
	prefix_0=$(basename "$file" .gz)
	prefix=$(basename "$prefix_0" .vcf)

	echo this is the prefix "$prefix"
	tabix -p vcf ${file}.vcf.gz

	# create a dir to save out
	mkdir ${prefix}_dir
 
	# launch snpEF           database     sample   
	java -jar snpEff.jar -v R64-1-1_sgd ${prefix}.vcf.gz > anotated_${prefix}.vcf -stats ${prefix}.txt -csvStats ${prefix}.csv
	

	# copy the file to out dir
	cp anotated_${prefix}.vcf out_vcf
	mv anotated_${prefix}.vcf ${prefix}_dir 
	mv ${prefix}.txt ${prefix}_dir
	mv ${prefix}.csv ${prefix}_dir
	mv ${prefix}.genes.txt ${prefix}_dir
	# mv the compresed file also
	mv ${prefix}.vcf.gz ${prefix}_dir
done
