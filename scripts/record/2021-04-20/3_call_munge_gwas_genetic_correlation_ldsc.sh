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
path_temporary=$(<"./processing_bipolar_metabolism.txt")
path_ldsc=$(<"./tools_ldsc.txt")

path_waller="$path_temporary/waller"
path_bipolar_metabolism="$path_waller/bipolar_metabolism"
path_scripts_record="$path_waller/bipolar_metabolism/scripts/record/2021-04-20"
path_promiscuity_scripts="$path_waller/promiscuity/scripts"
path_promiscuity_scripts_ldsc="$path_promiscuity_scripts/ldsc_genetic_heritability_correlation"
path_scripts_format="$path_waller/promiscuity/scripts/format_gwas_ldsc"

path_dock="$path_waller/dock"
path_genetic_reference="$path_dock/access/genetic_reference"
path_alleles="$path_genetic_reference/alleles"
path_disequilibrium="$path_genetic_reference/disequilibrium"
path_baseline="$path_genetic_reference/baseline"
path_weights="$path_genetic_reference/weights"
path_frequencies="$path_genetic_reference/frequencies"
path_gwas="$path_dock/gwas"
path_heritability="$path_dock/heritability"
path_genetic_correlation="$path_dock/genetic_correlation"


# Organize variables.
study_one="30124842_yengo_2018"
study_two="31043756_stahl_2019"
phenotype_one="body_mass_index"
phenotype_two="bipolar_disorder"
path_gwas_one="${path_gwas}/${study_one}/gwas_format.txt.gz"
path_gwas_two="${path_gwas}/${study_two}/gwas_format.txt.gz"
path_gwas_munge_one="${path_gwas}/${study_one}/gwas_munge"
path_gwas_munge_two="${path_gwas}/${study_two}/gwas_munge"
path_study_genetic_correlation="${path_genetic_correlation}/${study_one}"
mkdir -p $path_study_genetic_correlation
path_genetic_correlation_report="${path_study_genetic_correlation}/correlation_${study_two}"
path_genetic_correlation_report_suffix="${path_genetic_correlation_report}.log"
path_report=$path_genetic_correlation_report
/usr/bin/bash "$path_promiscuity_scripts_ldsc/munge_gwas_genetic_correlation_ldsc.sh" \
$phenotype_one \
$phenotype_two \
$path_gwas_one \
$path_gwas_two \
$path_gwas_munge_one \
$path_gwas_munge_two \
$path_report \
$path_alleles \
$path_disequilibrium \
$path_ldsc
