#!/bin/bash
mkdir -p out_vcf

parallel bgzip {} ::: *.vcf
wait

for file in *.vcf.gz; do
    prefix=$(basename "$file" .vcf.gz)
    echo "Processing: $prefix"

    tabix -p vcf "$file"                          
    mkdir -p "${prefix}_dir"

    java -jar snpEff.jar -v R64-1-1_sgd \    
        -stats "${prefix}.txt" \
        -csvStats "${prefix}.csv" \
        "$file" > "anotated_${prefix}.vcf"

    cp "anotated_${prefix}.vcf" out_vcf/
    mv "anotated_${prefix}.vcf" "${prefix}_dir/"
    mv "${prefix}.txt"       "${prefix}_dir/"
    mv "${prefix}.csv"       "${prefix}_dir/"
    mv "${prefix}.genes.txt" "${prefix}_dir/"
    mv "$file"               "${prefix}_dir/"
done
