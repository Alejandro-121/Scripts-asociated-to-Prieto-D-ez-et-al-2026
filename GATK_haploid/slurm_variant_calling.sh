#!/bin/bash
#SBATCH -n 4                # Number of cores requested
#SBATCH -N 1                # Number of nodes requested
#SBATCH -t 40:00            # Runtime in minutes.
#SBATCH --mem=8G            # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o hostname_%j.out  # Standard output goes to this file
#SBATCH -e hostname_%j.err  # Standard error goes to this file

################
# load modules #
################
module load GCCcore/11.2.0 # bwa, GATK & R module dependency
module load BWA/0.7.17
module load GATK/4.2.5.0-Java-11
module load picard/2.25.1-Java-11
module load GCC/11.2.0 # samtools dependency
module load SAMtools/1.14
module load OpenMPI/4.1.1 # R module dependency
module load R/4.1.2


#########
# Imput #
#########

# this parameters need to be redifined to be taken from a list in a array job
pref="test"
numth="4"
refGen="EC1118.fasta"
read1="244_R1.trimmed.fastq.gz"
read2="244_R2.trimmed.fastq.gz"

########################
# launch main pipeline #
########################

# part 1
python variant_calling_1.py -p "$pref" -numth "$numth" -r "$refGen" -r1 "$read1" -r2 "$read2"
# this does not work from python subprocess need the enviroment path
java -jar "$EBROOTPICARD/picard.jar" CollectAlignmentSummaryMetrics R="$refGen" I="${pref}sorted_dedup.bam" O="${pref}alignment_metrics.txt"
# part 2 
python variant_calling_2.py -p "$pref" -r "$refGen"

