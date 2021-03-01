#!/bin/bash

###########################################################################
# Specify arguments for qsub command.
# Note that bash does not interpret qsub parameters, which are bash comments.
# Bash will not expand variables in qsub parameters.
# Shell.
#$ -S /bin/bash
# Name of job.
#$ -N waller_ldsc
# Contact.
#$ -M tcameronwaller@gmail.com
### #$ -m abes
#$ -m e
# Standard output and error.
# Specify as arguments when calling qsub.
### -o "./out"
### -e "./error"
# Queue.
# "1-hour", "1-day", "4-day", "7-day", "30-day", "lg-mem"
#$ -q 1-hour
# Priority 0-15.
### -p -10
# Memory per iteration.
# Segmentation errors commonly indicate a memory error.
#$ -l h_vmem=10G
# Concurrent threads; assigns value to variable NSLOTS.
#$ -pe threaded 16
# Range of indices.
# Specify as argument when calling qsub.
# Array batch indices cannot start at zero.
### -t 1-100:1
# Limit on concurrent processes.
#$ -tc 10

###########################################################################
###########################################################################
###########################################################################
# This script organizes an array batch job for the Sun Grid Engine.
# mForge
# 120 worker nodes, each with 32 CPUs and 512 Gigabytes memory
# 3 worker nodes, each with 64 CPUs and 1.5 Terabytes memory
###########################################################################
###########################################################################
###########################################################################


###########################################################################
# Organize variables.
index=$((SGE_TASK_ID-1))
path_instances=$1
path_dock=$2
count=$3

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_ldsc=$(<"./tools_ldsc.txt")
path_access="$path_dock/access"
path_disequilibrium="$path_dock/access/disequilibrium"
path_baseline="$path_dock/access/baseline"
path_weights="$path_dock/access/weights"
path_frequencies="$path_dock/access/frequencies"
path_alleles="$path_dock/access/alleles"
path_metabolites="$path_dock/access/metabolites"
path_metabolite_summaries="$path_dock/access/metabolites/metabolites_meta"
path_heritability_metabolites="$path_dock/heritability/metabolites"

###########################################################################
# Execute procedure.

# Read instance.
readarray -t files < $path_instances
file=${files[$index]}
path_file="$path_metabolite_summaries/$file"

#file="$(basename $path_file)"
identifier="$(cut -d'.' -f1 <<<$file)"
#echo "metabolite identifier: " $identifier

# Extract and organize information from summary.
# Write information to new, temporary file.
cd $path_heritability_metabolites
echo "SNP A1 A2 N BETA P" > ${identifier}_new.txt
zcat $path_file | awk 'NR > 1 {print $1, $2, $3, $16, $8, $10}' >> ${identifier}_new.txt
#head -10 ${identifier}_new.txt

$path_ldsc/munge_sumstats.py \
--sumstats ${identifier}_new.txt \
--out ${identifier}_munge \
--merge-alleles $path_alleles/w_hm3.snplist

# Partitioned heritability by stratified LD score regression.
#https://github.com/bulik/ldsc/wiki/Partitioned-Heritability-from-Continuous-Annotations

$path_ldsc/ldsc.py \
--h2 $path_heritability_metabolites/${identifier}_munge.sumstats.gz \
--ref-ld-chr $path_baseline/baselineLD. \
--w-ld-chr $path_weights/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile-chr $path_frequencies/1000G_Phase3_frq/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--print-delete-vals \
--out ${identifier}_heritability

# Report.
echo "----------"
echo "task index: " $index
echo "identifier: " $identifier
hostname
date
