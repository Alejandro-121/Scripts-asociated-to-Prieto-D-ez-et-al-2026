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
   _____       _______ _  __        _            _ _              _____        _____ _______ ___  
  / ____|   /\|__   __| |/ /       (_)          | (_)            |  __ \ /\   |  __ \__   __|__ \ 
 | |  __   /  \  | |  | ' /   _ __  _ _ __   ___| |_ _ __   ___  | |__) /  \  | |__) | | |     ) |
 | | |_ | / /\ \ | |  |  <   | '_ \| | '_ \ / _ \ | | '_ \ / _ \ |  ___/ /\ \ |  _  /  | |    / / 
 | |__| |/ ____ \| |  | . \  | |_) | | |_) |  __/ | | | | |  __/ | |  / ____ \| | \ \  | |   / /_ 
  \_____/_/    \_\_|  |_|\_\ | .__/|_| .__/ \___|_|_|_| |_|\___| |_| /_/    \_\_|  \_\ |_|  |____|
                             | |     | |                                                          
                             |_|     |_|                                                          
                                                   
"""
program_description = "Runs the second part of my variant calling pipeline."
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
# ref gen
parser.add_argument("-r", "--refGen", help="Set reference genome (fasta file)",
                    required=True)

# load args
args = parser.parse_args()
pref = args.pref
refGen = args.refGen


# set wd
wd = os.getcwd() + "/"

# adjust parameters!! #
vc.filter_pre_bqsr_indels(refGen, pref) # tested
vc.filter_pre_bqsr_snps(refGen, pref) # tested

vc.select_variants(pref) # tested
vc.bqsr(refGen, pref) # tested

#####################################
# real variant calling starts now!! #
#####################################

vc.call_variants(refGen, pref) # tested
vc.snp_indel(refGen, pref) # tested
vc.filter_snps(refGen, pref) # tested
vc.filter_indels(refGen, pref) # tested 
