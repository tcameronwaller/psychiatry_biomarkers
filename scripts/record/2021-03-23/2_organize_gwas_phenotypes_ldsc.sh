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
path_genetic_reference="$path_dock/access/genetic_reference"
path_gwas="$path_dock/gwas"
path_heritability="$path_dock/heritability"
path_genetic_correlation="$path_dock/genetic_correlation"

###########################################################################
# Execute procedure.

# Initialize directories.
rm -r $path_gwas
rm -r $path_heritability
rm -r $path_genetic_correlation
if [ ! -d $path_gwas ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_gwas
fi
if [ ! -d $path_heritability ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_heritability
fi
if [ ! -d $path_genetic_correlation ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_genetic_correlation
fi

# Organize information in format for LDSC.

# Yengo et al, Human Molecular Genetics, 2018 (PubMed:30124842)
# phenotype: body mass index
study="30124842_yengo_2018"
source_file="Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz"
path_source_directory="${path_gwas_summaries}/${study}"
path_source_file="${path_source_directory}/${source_file}"
path_script_gwas_organization="${path_scripts_organization}/organize_gwas_ldsc_${study}.sh"
report="true" # "true" or "false"
/usr/bin/bash "$path_script_gwas_organization" \
$study \
$source_file \
$path_source_file \
$path_genetic_reference \
$path_gwas \
$path_heritability \
$path_genetic_correlation \
$path_promiscuity_scripts \
$report

# Pulit et al, Human Molecular Genetics, 2018 (PubMed:30239722)
# phenotype: waist-to-hip ratio, adjusted for BMI
study="30239722_pulit_2018"
source_file="whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz"
path_source_directory="${path_gwas_summaries}/${study}"
path_source_file="${path_source_directory}/${source_file}"
path_script_gwas_organization="${path_scripts_organization}/organize_gwas_ldsc_${study}.sh"
report="true" # "true" or "false"
/usr/bin/bash "$path_script_gwas_organization" \
$study \
$source_file \
$path_source_file \
$path_genetic_reference \
$path_gwas \
$path_heritability \
$path_genetic_correlation \
$path_promiscuity_scripts \
$report
