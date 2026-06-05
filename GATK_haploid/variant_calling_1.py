#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 12:56:27 2024

Will run GATK (part 1). Prepared to run in an HPC (tested on HPC drago)

This code depends on variant_calling_fun.py, here are containded the custom
functions that are required

@author: alejandro
"""
####################
# import functions #
####################
import os
import argparse

# import custom funcions
import variant_calling_fun as vc


###################
# parse arguments #
###################


# program info
name_logo = r"""
   _____       _______ _  __        _            _ _              _____        _____ _______   __ 
  / ____|   /\|__   __| |/ /       (_)          | (_)            |  __ \ /\   |  __ \__   __| /_ |
 | |  __   /  \  | |  | ' /   _ __  _ _ __   ___| |_ _ __   ___  | |__) /  \  | |__) | | |     | |
 | | |_ | / /\ \ | |  |  <   | '_ \| | '_ \ / _ \ | | '_ \ / _ \ |  ___/ /\ \ |  _  /  | |     | |
 | |__| |/ ____ \| |  | . \  | |_) | | |_) |  __/ | | | | |  __/ | |  / ____ \| | \ \  | |     | |
  \_____/_/    \_\_|  |_|\_\ | .__/|_| .__/ \___|_|_|_| |_|\___| |_| /_/    \_\_|  \_\ |_|     |_|
                             | |     | |                                                          
                             |_|     |_|                                                          
"""
program_description = "Runs the first part of my variant calling pipeline."
epilog = """It is run until picard CollectAlignmentSummaryMetrics.
    It's output needs to be reviewed before continuing to chose bqsr paramaters.
    Also on hpc drago picard can not be called from a python subprocess.
    It will be called in the bash script used a wraper to launch this
    program."""


# init parser
parser = argparse.ArgumentParser(prog=name_logo,
                                 description=program_description, 
                                 epilog=epilog
                                 )
# pref in
pref_help = """"Prefix to use in all output. All output will be genrated
 appending pref + _ + names + .extensions"""
parser.add_argument("-p", "--pref", help=pref_help, required=True)
# number threads
numth_help = "Select number threads to be used by bwa. Default is 4"
parser.add_argument("-numth", "--numberThreads", help=numth_help,
                    required=False, default=4, type=int)
# are you running a competitive aligment ?
parser.add_argument("-comp", "--isComp",
                    help = "Select true if you are doing a competive aligment will filter MQ strictly, 0 -> false, 1 -> ture, ",
                    required=False, default=0, type=int)
# ref gen
parser.add_argument("-r", "--refGen", help = "Set reference genome (fasta file)",
                    required=True) 
# read1
parser.add_argument("-r1", "--read1", help = "read1 in paired end, ",
                    required=True)
# read2
parser.add_argument("-r2", "--read2", help = "read2 in paired end, Illumina")

# load args
args = parser.parse_args()
pref = args.pref
numth = args.numberThreads
comp = args.isComp
refGen = args.refGen
read1 = args.read1
read2 = args.read2

# comp to logic
if comp == 0:
    comp = False
else:
    comp = True

# set wd
wd = os.getcwd() + "/"


############################
# launch the main pipeline #
############################

vc.index(refGen) # tested
vc.bwa_mem(pref, wd, str(numth), refGen, read1, read2) # tested
vc.sam_pre_gatk(pref, comp)
vc.mark_duplicates(pref) # tested

vc.call_variants_pre_bqsr(refGen, pref) #tested
vc.snp_indel_pre_bqsr(refGen, pref) # tested
