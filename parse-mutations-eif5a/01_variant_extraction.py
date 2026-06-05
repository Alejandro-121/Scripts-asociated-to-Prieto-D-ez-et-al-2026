#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 10:08:24 2024

@author: alejandro
"""

import re
import pandas as pd
import os

# os.chdir("/home/alejandro/Documentos/nextcloud/samoa/vcf/annotated/indels/")

# Define the sample names based on the VCF header
sample_names = ["2-1", "2-3", "sup.11", "sup.15", "sup.1", "sup.22", "sup.23", "sup.25", "sup.27", "sup.2", "WT"]

# Define the base columns for the table, including one for each sample
columns = ["CHROM", "POS", "REF", "ALT", "Effect", "Gene", "Transcript"] + [f"{sample}_mutation" for sample in sample_names] + [f"{sample}_allele" for sample in sample_names]

# Store rows as dictionaries for easy DataFrame creation later
rows = []

# Function to parse the ANN field for mutation effect annotations
def parse_ann_field(ann_field):
    effects = []
    annotations = ann_field.split(',')
    for ann in annotations:
        fields = ann.split('|')
        if len(fields) > 4:
            effects.append({
                "Effect": fields[1],
                "Gene": fields[3],
                "Transcript": fields[6],
            })
    return effects

# Process each line in the VCF
with open("merged.vcf", "r") as vcf:
    for line in vcf:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt_alleles = fields[4].split(",")  # Handle multiple ALT alleles
        info = fields[7]
        
        # Extract the ANN field for mutation effects
        ann_match = re.search(r'ANN=([^;]+)', info)
        if not ann_match:
            continue
        ann_field = ann_match.group(1)
        effects = parse_ann_field(ann_field)
        
        # Extract sample genotypes
        samples = fields[9:]
        
        # Process each ALT allele separately, creating a new row for each one
        for alt_index, alt_allele in enumerate(alt_alleles, start=1):
            # Create a separate row for each ALT allele
            for effect in effects:
                row = {
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt_allele,
                    "Effect": effect["Effect"],
                    "Gene": effect["Gene"],
                    "Transcript": effect["Transcript"],
                }
                
                # Initialize mutation and allele columns for each sample
                row.update({f"{sample}_mutation": 0 for sample in sample_names})
                row.update({f"{sample}_allele": ref for sample in sample_names})  # Default to REF allele

                # Update each sample's mutation presence and alleles
                for sample_name, sample_data in zip(sample_names, samples):
                    genotype = sample_data.split(":")[0]  # Extract GT field
                    alleles_present = []
                    if genotype not in ["0", "0/0", "."]:  # Check for non-reference alleles
                        alleles = genotype.replace("|", "/").split("/")  # Handle phased/unphased
                        for allele in alleles:
                            if allele.isdigit() and int(allele) == alt_index:  # Check if ALT matches this row's allele
                                alleles_present.append(alt_allele)
                            elif allele == "0":
                                alleles_present.append(ref)  # REF allele
                        # Mark mutation presence and allele(s)
                        if alt_allele in alleles_present:
                            row[f"{sample_name}_mutation"] = 1
                            row[f"{sample_name}_allele"] = "/".join(alleles_present)

                rows.append(row)

# Convert to DataFrame and save to CSV
df = pd.DataFrame(rows, columns=columns)
df.to_csv("mutation_wide_table_with_alleles.csv", index=False)
