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
path_ldsc=$(<"./tools_ldsc.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")

path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-05-04"

path_dock="$path_process/dock"
path_genetic_reference="${path_dock}/access/genetic_reference"
path_gwas="${path_dock}/gwas"
path_genetic_correlation="${path_dock}/genetic_correlation"

###########################################################################
# Execute procedure.

# Define metabolite study.
metabolite_study="27005778_kettunen_2016"
path_gwas_metabolite="${path_gwas}/${metabolite_study}"
gwas_metabolite_munge_suffix="_gwas_munge.sumstats.gz"
pattern_gwas_metabolite_munge_suffix="*${gwas_metabolite_munge_suffix}"

# Define main phenotype studies.
phenotype_studies=()
phenotype_studies+=("30124842_yengo_2018")
phenotype_studies+=("30239722_pulit_2018")
phenotype_studies+=("30482948_walters_2018_all")
phenotype_studies+=("30482948_walters_2018_eur")
phenotype_studies+=("30482948_walters_2018_eur_unrel")
phenotype_studies+=("30718901_howard_2019")
phenotype_studies+=("29906448_ruderfer_2018_scz_vs_ctl")
phenotype_studies+=("29906448_ruderfer_2018_scz_bpd_vs_ctl")
phenotype_studies+=("29906448_ruderfer_2018_scz_vs_bpd")
phenotype_studies+=("29906448_ruderfer_2018_bpd_vs_ctl")
phenotype_studies+=("00000000_ripke_2021")
phenotype_studies+=("31043756_stahl_2019")
phenotype_studies+=("00000000_mullins_2021_all")
phenotype_studies+=("00000000_mullins_2021_bpd1")
phenotype_studies+=("00000000_mullins_2021_bpd2")

# Organize information in format for LDSC.
for phenotype_study in "${phenotype_studies[@]}"; do
  # Organize paths.
  path_gwas_phenotype="${path_gwas}/${phenotype_study}"
  path_gwas_phenotype_format_compress="${path_gwas_phenotype}/gwas_format.txt.gz"
  path_gwas_phenotype_munge_suffix="${path_gwas_phenotype}/gwas_munge.sumstats.gz"
  path_genetic_correlation_comparison="${path_genetic_correlation}/${phenotype_study}/${metabolite_study}"
  # Initialize directories.
  rm -r $path_genetic_correlation_comparison
  mkdir -p $path_genetic_correlation_comparison
  # Organize variables.
  report="true" # "true" or "false"
  /usr/bin/bash "${path_scripts_record}/6_organize_phenotype_metabolites_correlations.sh" \
  $metabolite_study \
  $phenotype_study \
  $path_gwas_metabolite \
  $gwas_metabolite_munge_suffix \
  $pattern_gwas_metabolite_munge_suffix \
  $path_gwas_phenotype \
  $path_gwas_phenotype_munge_suffix \
  $path_genetic_correlation_comparison \
  $path_genetic_reference \
  $path_promiscuity_scripts \
  $path_scripts_record \
  $path_ldsc \
  $report
done
