#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_temporary=$(<"./processing_bipolar_metabolism.txt")

path_waller="$path_temporary/waller"
path_bipolar_metabolism="$path_waller/bipolar_metabolism"
path_scripts_organization="$path_waller/bipolar_metabolism/scripts/organization"
path_scripts_record="$path_waller/bipolar_metabolism/scripts/record/2021-03-23"
path_promiscuity_scripts="$path_waller/promiscuity/scripts"

path_dock="$path_waller/dock"
path_genetic_correlation="$path_dock/genetic_correlation"

###########################################################################
# Execute procedure.
# Organize information in format for LDSC.

# Yengo et al, Human Molecular Genetics, 2018 (PubMed:30124842)
# phenotype: body mass index
study="30124842_yengo_2018"
source_file="Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz"
path_source_directory="${path_gwas_summaries}/${study}"
path_source_file="${path_source_directory}/${source_file}"
path_product_directory="${path_genetic_correlation}/${study}"
path_script_gwas_organization="${path_scripts_organization}/organize_gwas_ldsc_${study}.sh"
path_gwas_format="${path_product_directory}/gwas_format_${study}.txt.gz"
report="true" # "true" or "false"
mkdir -p $path_product_directory
/usr/bin/bash "$path_script_gwas_organization" \
$study \
$source_file \
$path_source_file \
$path_gwas_format \
$path_product_directory \
$path_promiscuity_scripts \
$report

# Pulit et al, Human Molecular Genetics, 2018 (PubMed:30239722)
# phenotype: waist-to-hip ratio, adjusted for BMI
study="30239722_pulit_2018"
source_file="whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz"
path_source_directory="${path_gwas_summaries}/${study}"
path_source_file="${path_source_directory}/${source_file}"
path_product_directory="${path_genetic_correlation}/${study}"
path_script_gwas_organization="${path_scripts_organization}/organize_gwas_ldsc_${study}.sh"
path_gwas_format="${path_product_directory}/gwas_format_${study}.txt.gz"
report="true" # "true" or "false"
mkdir -p $path_product_directory
/usr/bin/bash "$path_script_gwas_organization" \
$study \
$source_file \
$path_source_file \
$path_gwas_format \
$path_product_directory \
$path_promiscuity_scripts \
$report
