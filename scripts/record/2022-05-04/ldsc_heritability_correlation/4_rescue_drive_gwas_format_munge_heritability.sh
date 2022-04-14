#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Rescue a few instances from batch that failed.
###########################################################################
###########################################################################
###########################################################################

# Parameters.
#cohorts_models="albumin_linear_1"
cohorts_models="testosterone_linear"

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"

path_gwas_concatenation_container="${path_dock}/gwas_concatenation/${cohorts_models}"
path_gwas_format_container="${path_dock}/gwas_ldsc_format/${cohorts_models}"
path_gwas_munge_container="${path_dock}/gwas_ldsc_munge/${cohorts_models}"
path_heritability_container="${path_dock}/heritability/${cohorts_models}"

#path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2022-05-04/ldsc_heritability_correlation"

################################################################################
# Organize argument variables.

# Parameters.

#study="male_age_high_joint_1_albumin_imputation" # successful rescue; TCW; 13 April 2022

#study="male_age_high_joint_1_testosterone_imputation" # unsuccessful rescue; LDSC Munge Failure; TCW; 13 April 2022
# "WARNING: median value of SIGNED_SUMSTATS is 0.24 (should be close to 0.0)."

study="male_age_high_unadjust_testosterone_imputation"

#study="male_age_middle_unadjust_testosterone_imputation"

name_gwas_concatenation_file="gwas_concatenation.txt.gz"
path_gwas_concatenation_compress="${path_gwas_concatenation_container}/${study}/${name_gwas_concatenation_file}"
regression_type="linear" # "linear" or "logistic"
response="coefficient" # "coefficient", "odds_ratio", or "z_score"; for linear GWAS, use "coefficient" unless "response_standard_scale" is "true", in which case "z_score"
response_standard_scale="false" # whether to convert reponse (effect, coefficient) to z-score standard scale ("true" or "false")
restore_target_study_directories="true" # whether to delete any previous directories for each study's format and munge GWAS ("true" or "false")

###########################################################################
# Execute procedure.

##############################################################################
# Format GWAS summary statistics for analysis in LDSC.
# Paths.
path_gwas_target_parent="${path_gwas_format_container}/${study}"
if [[ "$restore_target_study_directories" == "true" ]]; then
  rm -r $path_gwas_target_parent
fi
mkdir -p $path_gwas_target_parent
# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_drive_gwas_format="${path_promiscuity_scripts}/gwas_process/drive_gwas_format.sh"
path_script_gwas_format="${path_promiscuity_scripts}/gwas_process/format_gwas_ldsc/format_gwas_ldsc_plink_${regression_type}.sh"
##########
# Format adjustment.
# Parameters.
report="true" # "true" or "false"
/usr/bin/bash "${path_script_drive_gwas_format}" \
$path_gwas_concatenation_compress \
$path_gwas_target_parent \
$path_promiscuity_scripts \
$path_script_gwas_format \
$response_standard_scale \
$report

##############################################################################
# LDSC Munge and Heritability.
# Paths.
path_dock="$path_process/dock"
path_genetic_reference="${path_dock}/access/genetic_reference"
path_gwas_source_parent="${path_gwas_format_container}/${study}"
path_gwas_target_parent="${path_gwas_munge_container}/${study}"
path_heritability_parent="${path_heritability_container}/${study}"
if [[ "$restore_target_study_directories" == "true" ]]; then
  rm -r $path_gwas_target_parent # do not remove if same as format directory
  rm -r $path_heritability_parent
fi
mkdir -p $path_gwas_target_parent
mkdir -p $path_heritability_parent
# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_drive_ldsc_gwas_munge_heritability="${path_scripts_gwas_process}/drive_ldsc_gwas_munge_heritability.sh"
# Parameters.
report="true" # "true" or "false"
/usr/bin/bash "${path_script_drive_ldsc_gwas_munge_heritability}" \
$path_gwas_source_parent \
$path_gwas_target_parent \
$path_heritability_parent \
$path_genetic_reference \
$response \
$report





#
