#!/bin/bash

#chmod u+x script.sh

###########################################################################
###########################################################################
###########################################################################
# This script organizes directories and iteration instances then submits
# script "regress_metabolite_heritability.sh" to the Sun Grid Engine.

# version check: 2

###########################################################################
###########################################################################
###########################################################################

# Organize paths.
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-07-02/scramble"
path_dock="$path_process/dock"
path_cohorts_models="${path_dock}/organization/cohorts_models"
path_gwas="${path_dock}/gwas/cohorts_models_maf_0_pfilter" # run a few GWAS with minimal MAF 0

# Initialize directories.
#rm -r $path_gwas
mkdir -p $path_gwas

# Note:
# Each GWAS (30,000 - 200,000 persons; 22 chromosomes) requires about 5-7 hours to run on the grid.

# General parameters.
threads=32 # 32, needs to match the "pe threaded" argument to scheduler
maf=0.0 # run on all SNPs and filter in subsequent analyses
chromosomes=22 # 22 # Count of chromosomes on which to run GWAS

#cohorts_models+=("white_bipolar_disorder_control_unadjust;table_white_bipolar_disorder_control;")
#path_table_phenotypes_covariates="${path_cohorts_models}/table_white_bipolar_disorder_control_body_mass_index.tsv"

# first priority... started at 10:22 on 2 July 2021
#path_table_phenotypes_covariates="${path_cohorts_models}/table_white_bipolar_disorder_case_body_mass_index.tsv"
#path_report="${path_gwas}/white_bipolar_disorder_case_unadjust_body_mass_index"

# second priority... started at ___ on 2 July 2021
path_table_phenotypes_covariates="${path_cohorts_models}/table_white_bipolar_disorder_control_body_mass_index.tsv"
path_report="${path_gwas}/white_bipolar_disorder_control_unadjust_body_mass_index"

phenotypes="body_mass_index"
covariates="genotype_pc_1,genotype_pc_2,genotype_pc_3,genotype_pc_4,genotype_pc_5,genotype_pc_6,genotype_pc_7,genotype_pc_8,genotype_pc_9,genotype_pc_10"

# Initialize directories.
rm -r $path_report
mkdir -p $path_report

# Submit array batch to Sun Grid Engine.
# Array batch indices cannot start at zero.
echo "----------------------------------------------------------------------"
echo "Submit array batch to Sun Grid Engine."
echo "----------------------------------------------------------------------"
qsub -t 1-${chromosomes}:1 \
-o "$path_report/out.txt" -e "$path_report/error.txt" \
$path_scripts_record/2_run_batch_gwas.sh \
$path_table_phenotypes_covariates \
$path_report \
$phenotypes \
$covariates \
$threads \
$maf \
