#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 17 May 2024
# Date, last execution: 17 May 2024
# Date, review: __ May 2024
################################################################################
# Note

# TODO: TCW; 4 April 2023
# TODO: write and call a python script to extract the LDSC estimates of heriability and genetic correlation
# TODO: follow the pattern of the QQ-Plot scripts.

################################################################################



################################################################################
# Organize paths.

# Identifiers or designators of parameter version, preparation batch, and
# analysis batch.
identifier_analysis="gwas_2023-12-30_ldsc_2024-01-08_extra_2024-05-15"
identifier_parameter="tcw_2023-12-30_dbsnp_rsid"

# Directories.
cd ~/paths
path_waller_tools=$(<"./waller_tools.txt")
path_directory_environment="${path_waller_tools}/python/environments/main"
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_temporary="${path_directory_process}/_temporary_blargish_tcw_2024-05-17_"
path_directory_dock="${path_directory_process}/dock"

path_directory_group_parent="${path_directory_dock}/${identifier_analysis}"

path_directory_source_h2="${path_directory_group_parent}/5_gwas_heritability_ldsc"
#path_directory_source_h2_no_liability="${path_directory_group_parent}/5_gwas_heritability_ldsc_no_liability"
path_directory_product_h2="${path_directory_group_parent}/5_gwas_heritability_ldsc"

#path_directory_source_rg="${path_directory_group_parent}/6_gwas_correlation_ldsc_all"
path_directory_source_rg="${path_directory_group_parent}/6_gwas_correlation_ldsc_primary"
#path_directory_source_rg="${path_directory_group_parent}/6_gwas_correlation_ldsc_secondary"
#path_directory_source_rg="${path_directory_group_parent}/6_gwas_correlation_ldsc_primary_secondary"
#path_directory_product_rg="${path_directory_group_parent}/6_gwas_correlation_ldsc_all_extraction"
path_directory_product_rg="${path_directory_group_parent}/6_gwas_correlation_ldsc_primary_extraction"
#path_directory_product_rg="${path_directory_group_parent}/6_gwas_correlation_ldsc_secondary_extraction"
#path_directory_product_rg="${path_directory_group_parent}/6_gwas_correlation_ldsc_primary_secondary_extraction"

# Scripts.
path_file_script_extract_h2_rg="${path_directory_process}/partner/scripts/ldsc/extract_ldsc_heritability_correlation.sh"

# Initialize directories.
rm -r $path_directory_temporary # caution
rm -r $path_directory_product_h2 # caution
rm -r $path_directory_product_rg # caution
mkdir -p $path_directory_temporary
mkdir -p $path_directory_product_h2
mkdir -p $path_directory_product_rg

################################################################################
# Organize parameters.

name_file_source_prefix="none" # must not be empty string
name_file_source_suffix=".log" # must not be empty string
name_file_source_not="...place_holder" # exclude any files that include this character string in file name
report="true"

################################################################################
# Execute procedure.

##########
# SNP heritability.
type_analysis="heritability"
traversal="false" # whether to extract from all files in child source directories, preserving names of child directories
name_file_product="table_heritability"
#name_file_product="table_heritability_no_liability"
/usr/bin/bash $path_file_script_extract_h2_rg \
$type_analysis \
$path_directory_source_h2 \
$traversal \
$name_file_source_prefix \
$name_file_source_suffix \
$name_file_source_not \
$name_file_product \
$path_directory_product_h2 \
$path_directory_process \
$path_directory_temporary \
$path_directory_environment \
$report

##########
# Genetic correlation.
type_analysis="correlation"
traversal="true" # whether to extract from all files in child source directories, preserving names of child directories
name_file_product="none"
/usr/bin/bash $path_file_script_extract_h2_rg \
$type_analysis \
$path_directory_source_h2 \
$traversal \
$name_file_source_prefix \
$name_file_source_suffix \
$name_file_source_not \
$name_file_product \
$path_directory_product_rg \
$path_directory_process \
$path_directory_temporary \
$path_directory_environment \
$report

##########
# Remove temporary, intermediate files.
rm -r $path_directory_temporary

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "7_extract_ldsc_heritability_correlation.sh"
  echo "----------"
fi



#
