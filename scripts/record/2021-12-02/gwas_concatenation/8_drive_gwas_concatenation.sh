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
pattern_gwas_report_file=${2} # string glob pattern by which to recognize PLINK2 GWAS report files
path_gwas_source_container=${3} # full path to parent directories of GWAS summary statistics for each study
path_gwas_target_container=${4} # full path to parent directories of GWAS summary statistics for each study
path_scripts_record=${5} # full path to directory of scripts for a specific analysis report date
path_process=${6} # full path to directory for all processes relevant to current project
chromosome_x=${7} # whether to collect GWAS summary statistics report for Chromosome X

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
path_script_drive_gwas_concatenation="${path_scripts_gwas_process}/drive_gwas_concatenation_compression.sh"
# Parameters.
report="true" # "true" or "false"
/usr/bin/bash "$path_script_drive_gwas_concatenation" \
$pattern_gwas_report_file \
$path_gwas_source_parent \
$path_gwas_target_parent \
$path_promiscuity_scripts \
$chromosome_x \
$report



##############################################################################
