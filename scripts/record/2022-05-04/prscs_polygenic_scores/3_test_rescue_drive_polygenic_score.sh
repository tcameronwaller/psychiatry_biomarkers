#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Rescue a few instances from batch that failed.
###########################################################################
###########################################################################
###########################################################################

# Parameters.
gwas_group="steroid_globulin_linear_1"
study="female_joint_1_steroid_globulin_imputation_log"
count_gwas_samples=174949 # count of samples in GWAS
file_name_prefix="test_shbg_female"

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_tools=$(<"./waller_tools.txt")
path_prs_cs="${path_tools}/prs_cs/"
path_snp_relevance_bim_prefix="${path_prs_cs}/PRScsx/test_data/test"

path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"
path_genetic_reference_prscs="${path_dock}/access/genetic_reference_prscs"
path_gwas_container="${path_dock}/gwas_format_prscs/${gwas_group}/${study}"
path_source_gwas_summary="${path_gwas_container}/gwas_format.txt.gz"
path_target_directory_score="${path_dock}/polygenic_score_test"

###########################################################################
# Execute procedure.

# Initialize target directory.
rm -r $path_target_score_directory
mkdir -p $path_target_score_directory

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_drive_ldsc_gwas_munge_heritability="${path_scripts_gwas_process}/polygenic_score/drive_test_prscs.sh"

# Parameters.
# chromosome=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,x,xy
population_ancestry="EUR"
chromosome=22
report="true" # "true" or "false"

/usr/bin/bash "${path_script_drive_test_prscs}" \
$path_genetic_reference_prscs \
$path_source_gwas_summary \
$count_gwas_samples \
$population_ancestry \
$path_snp_relevance_bim_prefix \
$path_target_directory_score \
$file_name_prefix \
$chromosome \
$report





#
