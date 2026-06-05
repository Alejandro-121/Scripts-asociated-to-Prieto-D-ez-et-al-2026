#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 09:16:40 2024

@author: alejandro
"""

import subprocess
import os

######################
# documentation used #
######################

# Dabao's shared pipeline
# gatk videos and documentaion
# https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/

##########################
# dependcies call module #
##########################

# bwa
# gatk
# picard
# samtools

##########################
##########################
## fuctions of pipeline ##
##########################
##########################



#####################
# map to reference  #
#####################

# index
def index(refGen):
    """
    Calls bwa index to index ref genome creating all index file .amb, .bwt
    .sa .pac 

    Parameters
    ----------
    refGen : str
        Path to ref genome.

    Returns
    -------
    None. All index file with the same name as imput appendig the index
    extension. Like test.fa --> test.fa.amb

    """
    # -p prefix -a is the algorithim to use "is" is apropiate for short genomes
    # like s.cer
    
    subprocess.call(["gatk", "CreateSequenceDictionary", "-R", refGen])
    subprocess.call(["samtools", "faidx", refGen])
    subprocess.call(["bwa", "index", "-p", refGen, "-a", "is", refGen])


# mapping
def bwa_mem(pref, wd, numth, refGen, read1, read2 = None):
    """
    Funcion to call bwa-mem to create an alingment agains a ref genome
    creating a sam file
    
    Parameters
    ----------
    pref: str
        prefix of file in use will be used to create output file name.
        .sam extension will be added
    wd : str
        Path to working dir.
    numth : str (can also accept coercible to str like int)
        Number of threads to use by bwa
    refGen : str
        Path to reference genome
    read1 : str
        Path to read 1 or only read in single end.
    read2 : str
        Path to read 2 (paired end). Optional default None 

    Returns
    -------
    None. The output is writen to file sistem (alinged sam file)

    """
    # append extension
    bwaOutFile = pref + ".sam"
    # create stout object to write to
    bwaOutFile = open(wd + bwaOutFile, "w")
    # coerce to str
    numth = str(numth)
    
    # define readgroups needed for gatk, 
    # some reduntant info is placed, can be used for more complex experiments
    readGroup = "@RG" + "\\tID:" + pref + "\\tSM:" + pref + "\\tLB:" + pref + \
    "\\tPU:" + "ILLUMINA" + "\\tPL:" + "ILLUMINA" 


    if read2:
   # here stout is needed to redict de output bwa spits outs the alingment
   # to terminal usualy redirect to a file using ">"
        subprocess.call(["bwa", "mem", "-R", readGroup, "-t", numth, refGen,
                         read1, read2],
                        stdout=bwaOutFile, cwd=wd)
    else:
        subprocess.call(["bwa", "mem", "-R", readGroup,"-t", numth, refGen,
                         read1],
                        stdout=bwaOutFile, cwd=wd)
        
    
# filtering low quality????????????


    
def sam_pre_gatk(pref, comp = False):
    """
    Calls samtools, sorts, converts to bam and index.

    Parameters
    ----------
    pref : str
        prefix of file in use will be used to take .sam imput and create
        sorted.bam and index it.

    Returns
    -------
    None.
    Writes to file sistem a bamfile sorted and its index

    """
    samInFile = pref + ".sam"
    bamOutFile = pref + "sorted.bam"
    
    # sort
    subprocess.call(["samtools", "sort", samInFile, "-o", bamOutFile])
    # filter important for competitive mapping
    if comp == True:
        subprocess.call(["samtools", "view", "-h", "-q", "3", samInFile, samInFile])
    # index
    subprocess.call(["samtools", "index", bamOutFile])
    

###################
# mark duplicates #
###################
def mark_duplicates(pref):
    """
    Need as imput a .bam sorted (GATLK also takes .sam)
    Calls mark duplicates from GATK.
    Duplicates create false confidence in coverage of a region.
    Can be sequencing artifacts. This step marks them so next steps can take
    this information in consideration. It does not remove any read just in case
    those are needed at some point.

    Parameters
    ----------
    pref : str
        prefix of file in use will be used to take .sam imput and create
        depup.bam output.

    Returns
    -------
    The output is writen to file sistem .bam
    """
    bamInFile = pref + "sorted.bam"
    bamOutFile = pref + "sorted_dedup.bam"
    metrixFile = pref + "_Dupmetrix.txt"
    
    
    subprocess.call(["gatk", "MarkDuplicates", "-I", bamInFile, "-O",
                     bamOutFile, "-M", metrixFile])
    
    # index
    subprocess.call(["samtools", "index", bamOutFile])

####################################
# collect aligment summary metrics #
####################################

# not working on drago problem with picard path !!!!!!!!!!!!!!!!!!!!!!
#def collect_alingnmen_summary_metrics(refGen, pref,
#        #default path is the one used in HPC drago
#        path2picard = "$EBROOTPICARD/picard.jar"):
#    
#    I = pref + "sorted_dedup.bam" 
#    O1 = pref + "alignment_metrics.txt"
#    
#    subprocess.call(["java -jar " + path2picard,
#                     "CollectAlignmentSummaryMetrics", "R=", refGen,
#                     "I=" + I, "O=" + O1])
       
    
##############################
# CALL VARIANTS use pre bqsr #
##############################

def call_variants_pre_bqsr(refGen, pref):
    I = pref + "sorted_dedup.bam"
    O = pref + "pre_bqsr_variants.vcf"
    
    subprocess.call(["gatk", "HaplotypeCaller", "-R",
                     refGen, "-I", I, "-O", O])

##################################
# Extract SNPs & indels pre bqsr #
##################################

def snp_indel_pre_bqsr(refGen, pref):
    V = pref + "pre_bqsr_variants.vcf"
    O1 = pref + "pre_bqsr_snps.vcf"
    O2 = pref + "pre_bqsr_indels.vcf"
    
    # snps
    subprocess.call(["gatk", "SelectVariants", "-R", refGen, "-V", V,
                     "--select-type-to-include", "SNP", "-O", O1])
    # variants
    subprocess.call(["gatk", "SelectVariants", "-R", refGen, "-V", V,
                     "--select-type-to-include", "INDEL", "-O", O2])

def filter_pre_bqsr_snps(refGen, pref,
                    # filtering parameters
                    QD = "2.0", FS = "60.0", MQ = "40.0", SOR = "4.0", 
                    MQ_rank = "-12.5", Read_pos = "-8.0"):
    # snps
    V = pref + "pre_bqsr_snps.vcf"
    O = pref + "pre_bqsr_filtered_snps.vcf"
    
    # filtering parameters
    QD = "QD < " + str(QD)
    FS = "FS > " + str(FS)
    MQ = "MQ < " + str(MQ)
    SOR = "SOR > " + str(SOR)
    MQ_rank = "MQRankSum < " + str(MQ_rank)
    Read_pos = "ReadPosRankSum < " + str(Read_pos)
    
    subprocess.call(["gatk", "VariantFiltration",
                     "-R", refGen, "-V", V, "-O", O,
                     # filtering parameters
                     "-filter-name", "QD_filter", "-filter", QD,
                     "-filter-name", "FS_filter", "-filter", FS,
                     "-filter-name", "MQ_filter", "-filter", MQ,
                     "-filter-name", "SOR_filter", "-filter", SOR,
                     "-filter-name", "MQRankSum_filter", "-filter", MQ_rank,
                     "-filter-name", "ReadPosRankSum_filter", "-filter", Read_pos
                     ])
    # known error ' undefined variable ReadPosRankSum
    # https://gatk.broadinstitute.org/hc/en-us/community/posts/360057857352-VariantFiltration-ReadPosRankSum-MQRankSum-undefined
    
def filter_pre_bqsr_indels(refGen, pref,
                           QD = "2.0", FS = "200.0", SOR = "10.0"):
    
    # snps
    V = pref + "pre_bqsr_indels.vcf"
    O = pref + "pre_bqsr_filtered_indels.vcf"
    
    # filtering parameters
    QD = "QD < " + str(QD)
    FS = "FS > " + str(FS)
    SOR = "SOR > " + str(SOR)
                            
    subprocess.call(["gatk", "VariantFiltration",
                     "-R", refGen, "-V", V, "-O", O,
                     "-filter-name", "QD_filter", "-filter", QD,
                     "-filter-name", "FS_filter", "-filter", FS,
                     "-filter-name", "SOR_filter", "-filter", SOR
                     ])
    
def select_variants(pref):
    V1 = pref + "pre_bqsr_filtered_snps.vcf"
    V2 = pref + "pre_bqsr_filtered_indels.vcf"
    
    # output prepared for bqsr
    O1 = pref + "bqsr_in_snps.vcf"
    O2 = pref + "bqsr_in_indels.vcf"
    
    # snps
    subprocess.call(["gatk", "SelectVariants", "--exclude-filtered",
                     "-V", V1, "-O", O1 
                     ])
    # indels
    subprocess.call(["gatk", "SelectVariants", "--exclude-filtered",
                     "-V", V2, "-O", O2 
                     ])
    
def bqsr(refGen, pref):
    I = pref + "sorted_dedup.bam"
    bqsr_s = pref + "bqsr_in_snps.vcf"
    bqsr_i = pref + "bqsr_in_indels.vcf"
    report = pref + "bqsr_report.txt"
    post_report = pref + "post_bqsr_report.txt"
    plots = pref + "plots_post_bqsr_report.pdf"
    bqsred = pref + "bqsred.bam"
    
    # calculate
    subprocess.call(["gatk", "BaseRecalibrator",
                     "-R", refGen, "-I", I,
                     "--known-sites", bqsr_s,
                     "--known-sites", bqsr_i,
                     "-O", report
                     ])
    # apply to bam file
    subprocess.call(["gatk", "ApplyBQSR",
                     "-R", refGen, "-I", I, "-bqsr", report, "-O", bqsred
                     ])
    
    # rerun bqsr to be able to produce recalibration report
    subprocess.call(["gatk", "BaseRecalibrator",
                     "-R", refGen, "-I", bqsred,
                     "--known-sites", bqsr_s,
                     "--known-sites", bqsr_i,
                     "-O", post_report
        ])
    
    # recalibration report
    subprocess.call(["gatk", "AnalyzeCovariates",
                     "-before", report,
                     "-after", post_report,
                     "-plots", plots
                     ])
    
def call_variants(refGen, pref):
    
    I = pref + "bqsred.bam"
    O = pref + "bqsred_variants.vcf"
    
    subprocess.call(["gatk", "HaplotypeCaller",
                    "-R", refGen, "-I", I, "-O", O])
    

def snp_indel(refGen, pref):
    V = pref + "bqsred_variants.vcf"
    O1 = pref + "bqsred_snps.vcf"
    O2 = pref + "bqsred_indels.vcf"
    
    # snps
    subprocess.call(["gatk", "SelectVariants", "-R", refGen, "-V", V,
                     "--select-type-to-include", "SNP", "-O", O1])
    # variants
    subprocess.call(["gatk", "SelectVariants", "-R", refGen, "-V", V,
                     "--select-type-to-include", "INDEL", "-O", O2])


def filter_snps(refGen, pref,
                    # filtering parameters
                    QD = "2.0", FS = "60.0", MQ = "40.0", SOR = "4.0", 
                    MQ_rank = "-12.5", Read_pos = "-8.0"):
    # snps
    V = pref + "bqsred_snps.vcf"
    O = pref + "bqsred_filtered_snps.vcf"
    
    # filtering parameters
    QD = "QD < " + str(QD)
    FS = "FS > " + str(FS)
    MQ = "MQ < " + str(MQ)
    SOR = "SOR > " + str(SOR)
    MQ_rank = "MQRankSum < " + str(MQ_rank)
    Read_pos = "ReadPosRankSum < " + str(Read_pos)
    
    subprocess.call(["gatk", "VariantFiltration",
                     "-R", refGen, "-V", V, "-O", O,
                     # filtering parameters
                     "-filter-name", "QD_filter", "-filter", QD,
                     "-filter-name", "FS_filter", "-filter", FS,
                     "-filter-name", "MQ_filter", "-filter", MQ,
                     "-filter-name", "SOR_filter", "-filter", SOR,
                     "-filter-name", "MQRankSum_filter", "-filter", MQ_rank,
                     "-filter-name", "ReadPosRankSum_filter", "-filter", Read_pos])


def filter_indels(refGen, pref,
                           QD = "2.0", FS = "200.0", SOR = "10.0"):
    
    # snps
    V = pref + "bqsred_indels.vcf"
    O = pref + "bqsred_filtered_indels.vcf"
    
    # filtering parameters
    QD = "QD < " + str(QD)
    FS = "FS > " + str(FS)
    SOR = "SOR > " + str(SOR)
                            
    subprocess.call(["gatk", "VariantFiltration",
                     "-R", refGen, "-V", V, "-O", O,
                     "-filter-name", "QD_filter", "-filter", QD,
                     "-filter-name", "FS_filter", "-filter", FS,
                     "-filter-name", "SOR_filter", "-filter", SOR
                     ])



########
# main #
########

# imput

#pref = "test"
#wd = os.getcwd() + "/"
#numth = 4
#refGen = "EC1118.fasta"
#read1 = "244_R1.trimmed.fastq.gz"
#read2 = "244_R2.trimmed.fastq.gz"


# pipeline

# index(refGen) # tested
# bwa_mem(pref, wd, str(numth), refGen, read1, read2) # tested
# sam_pre_gatk(pref) # tested
# mark_duplicates(pref) # tested

# call_variants_pre_bqsr(refGen, pref) #tested
# snp_indel_pre_bqsr(refGen, pref) # tested


# collect_alingnmen_summary_metrics(refGen, pref) # not working ### subtitute by a bash script
############################################################################
# decide filtering parametes from collect alingment summary metrics output # <-- important checkpoint
############################################################################

# adjust parameters!! #
# filter_pre_bqsr_indels(refGen, pref) # tested
# filter_pre_bqsr_snps(refGen, pref) # tested

# select_variants(pref) # tested
# bqsr(refGen, pref) # tested

#####################################
# real variant calling starts now!! #
#####################################

# call_variants(refGen, pref) # tested
# snp_indel(refGen, pref) # tested
# filter_snps(refGen, pref) 
# filter_indels(refGen, pref)
