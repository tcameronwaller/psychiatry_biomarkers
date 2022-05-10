#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Test procedure or rescue.
###########################################################################
###########################################################################
###########################################################################

# Parameters.
gwas_group="steroid_globulin_linear_1"

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"

path_gwas_concatenation_container="${path_dock}/gwas_concatenation/${gwas_group}"
path_gwas_format_container="${path_dock}/gwas_format_prscs/${gwas_group}"

#path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2022-05-04/ldsc_heritability_correlation"

################################################################################
# Organize argument variables.

# Parameters.

study="female_joint_1_steroid_globulin_imputation_log"

name_gwas_concatenation_file="gwas_concatenation.txt.gz"
path_gwas_concatenation_compress="${path_gwas_concatenation_container}/${study}/${name_gwas_concatenation_file}"
regression_type="linear" # "linear" or "logistic"
response="coefficient" # "coefficient", "odds_ratio", or "z_score"; for linear GWAS, use "coefficient" unless "response_standard_scale" is "true", in which case "z_score"
restore_target_study_directories="true" # whether to delete any previous directories for each study's format and munge GWAS ("true" or "false")

###########################################################################
# Execute procedure.

##############################################################################
# Format GWAS summary statistics for analysis in PRS-CS.
# Paths.
path_gwas_target_parent="${path_gwas_format_container}/${study}"
if [[ "$restore_target_study_directories" == "true" ]]; then
  rm -r $path_gwas_target_parent
fi
mkdir -p $path_gwas_target_parent
# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_drive_format_gwas="${path_promiscuity_scripts}/gwas_process/prscs_polygenic_score/drive_format_gwas_prscs.sh"
path_script_format_gwas="${path_promiscuity_scripts}/gwas_process/format_gwas_prscs/format_gwas_prscs_plink_${regression_type}.sh"
##########
# Format adjustment.
# Parameters.
report="true" # "true" or "false"
/usr/bin/bash "${path_script_drive_format_gwas}" \
$path_gwas_concatenation_compress \
$path_gwas_target_parent \
$path_promiscuity_scripts \
$path_script_format_gwas \
$report

#
