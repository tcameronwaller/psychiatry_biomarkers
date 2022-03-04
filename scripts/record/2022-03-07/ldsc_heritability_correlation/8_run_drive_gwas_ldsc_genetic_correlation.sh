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

comparison_instance=${1} # text list of information for a single comparison instance
path_genetic_correlation_container=${2} # full path to parent directory for genetic correlation comparisons
path_genetic_reference=${3} # full path to directory for all processes relevant to current project
path_process=${4} # full path to directory for all processes relevant to current project
report=${5} # full path to directory for all processes relevant to current project

##############################################################################
# Extract details for comparison.
IFS=";" read -r -a array <<< "${comparison_instance}"
comparison_container="${array[0]}"
study_primary="${array[1]}"
path_gwas_primary_munge_suffix="${array[2]}"
study_secondary="${array[3]}"
path_gwas_secondary_munge_suffix="${array[4]}"
echo "----------"
echo "comparison container: ${comparison_container}"
echo "primary study: ${study_primary}"
echo "path: ${path_gwas_primary_munge_suffix}"
echo "secondary study: ${study_secondary}"
echo "path: ${path_gwas_secondary_munge_suffix}"

##############################################################################
# LDSC Genetic Correlation.
# Paths.
path_genetic_correlation_parent="${path_genetic_correlation_container}/${comparison_container}/${study_primary}/${study_secondary}"
rm -r $path_genetic_correlation_parent
mkdir -p $path_genetic_correlation_parent
# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_drive_ldsc_gwas_genetic_correlation="${path_scripts_gwas_process}/drive_ldsc_gwas_genetic_correlation.sh"
/usr/bin/bash "$path_script_drive_ldsc_gwas_genetic_correlation" \
$path_gwas_primary_munge_suffix \
$path_gwas_secondary_munge_suffix \
$path_genetic_correlation_parent \
$path_genetic_reference \
$report
