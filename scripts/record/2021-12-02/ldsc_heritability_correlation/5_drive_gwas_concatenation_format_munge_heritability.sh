#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize argument variables.

study=${1} # text name of GWAS study and also name of study's parent directory
path_gwas_concatenation_compress=${2} # full path to file for source GWAS summary statistics after concatenation and compression
response=${3} # whether GWAS response is beta coefficient ("coefficient"), odds ratio ("odds_ratio"), or z-scores ("z_score")
response_standard_scale=${4} # whether to convert response (coefficient) to z-score standard scale
path_gwas_format_container=${5} # full path to container directory of format GWAS summary statistics for each study
path_gwas_munge_container=${6} # full path to container directory of munge GWAS summary statistics for each study
path_heritability_container=${7} # full path to parent directory for heritability reports
path_process=${8} # full path to directory for all processes relevant to current project

###########################################################################
# Execute procedure.

##############################################################################
# Format GWAS summary statistics for analysis in LDSC.
# Paths.
path_gwas_target_parent="${path_gwas_format_container}/${study}"
rm -r $path_gwas_target_parent
mkdir -p $path_gwas_target_parent
# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_drive_gwas_format="${path_promiscuity_scripts}/gwas_process/drive_gwas_format.sh"
path_script_gwas_format="${path_promiscuity_scripts}/gwas_process/format_gwas_ldsc/format_gwas_ldsc_plink_linear.sh"
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
rm -r $path_gwas_target_parent # do not remove if same as format directory
rm -r $path_heritability_parent
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
