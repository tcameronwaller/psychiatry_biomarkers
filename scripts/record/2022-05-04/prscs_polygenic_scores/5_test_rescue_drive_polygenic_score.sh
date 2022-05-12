#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Rescue a few instances from batch that failed.
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize parameters.
gwas_group="steroid_globulin_linear_2"
study="female_male_priority_male_joint_1_steroid_globulin_imputation_log"
count_gwas_samples=329617 # count of samples in GWAS
file_name_prefix="test_shbg_female"
chromosome_x="true" # 'true' or 'false'
population_ancestry="EUR"
report="true" # "true" or "false"

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_tools=$(<"./waller_tools.txt")
path_prs_cs="${path_tools}/prs_cs/"
#path_snp_relevance_bim_container="${path_prs_cs}/PRScsx/test_data/test"

path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"
path_genetic_reference_prscs="${path_dock}/access/genetic_reference_prscs"
path_gwas_container="${path_dock}/gwas_format_prscs/${gwas_group}/${study}"
path_source_gwas_summary="${path_gwas_container}/gwas_format.txt"
path_snp_relevance_bim_container="${path_dock}/access/mayo_bipolar_genotype/snp_relevance_bim"
path_target_directory_score="${path_dock}/polygenic_score_test"

###########################################################################
# Execute procedure.

# Initialize target directory.
rm -r $path_target_directory_score
mkdir -p $path_target_directory_score

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_drive_test_prscs="${path_scripts_gwas_process}/prscs_polygenic_score/drive_estimate_prscs_allelic_effects.sh"

/usr/bin/bash "${path_script_drive_test_prscs}" \
$path_source_gwas_summary \
$count_gwas_samples \
$path_snp_relevance_bim_container \
$path_target_directory_score \
$file_name_prefix \
$path_genetic_reference_prscs \
$population_ancestry \
$chromosome_x \
$path_promiscuity_scripts \
$report





#
