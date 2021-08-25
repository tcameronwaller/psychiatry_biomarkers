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
path_gwas_source_container=${2} # full path to parent directories of GWAS summary statistics for each study
path_gwas_target_container=${3} # full path to parent directories of GWAS summary statistics for each study

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths

path_process=$(<"./process_sexy_alcohol.txt")
path_dock="$path_process/dock"

# General paths.
#path_gwas_source_container="${path_dock}/gwas_complete/cohorts_models"
#path_gwas_target_container="${path_dock}/gwas_process/cohorts_models"

###########################################################################
# Execute procedure.

##############################################################################
# Concatenation and format.
# Paths.
path_gwas_source_parent="${path_gwas_source_container}/${study}"
path_gwas_target_parent="${path_gwas_target_container}/${study}"
rm -r $path_gwas_target_parent
mkdir -p $path_gwas_target_parent
# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_drive_gwas_concatenation_format="${path_scripts_gwas_process}/drive_gwas_concatenation_format.sh"
# Parameters.
report="true" # "true" or "false"
/usr/bin/bash "$path_script_drive_gwas_concatenation_format" \
$path_gwas_source_parent \
$path_gwas_target_parent \
$path_promiscuity_scripts \
$report

##############################################################################
# LDSC Munge and Heritability.
# Paths.
path_genetic_reference="${path_dock}/access/genetic_reference"
path_gwas_source_parent="${path_gwas_target_container}/${study}"
path_gwas_target_parent="${path_gwas_target_container}/${study}"
path_heritability_parent="${path_dock}/heritability/cohorts_models/${study}"
rm -r $path_heritability_parent
mkdir -p $path_gwas_target_parent
mkdir -p $path_heritability_parent
# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_drive_ldsc_gwas_munge_heritability="${path_scripts_gwas_process}/drive_ldsc_gwas_munge_heritability.sh"
# Parameters.
report="true" # "true" or "false"
/usr/bin/bash "$path_script_drive_ldsc_gwas_munge_heritability" \
$path_gwas_source_parent \
$path_gwas_target_parent \
$path_heritability_parent \
$path_genetic_reference \
$report
