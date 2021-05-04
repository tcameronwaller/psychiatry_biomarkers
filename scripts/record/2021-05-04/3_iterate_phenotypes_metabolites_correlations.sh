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
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")

path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-05-04"

path_dock="$path_process/dock"

path_genetic_reference="${path_dock}/access/genetic_reference"
path_alleles="$path_genetic_reference/alleles"
path_disequilibrium="$path_genetic_reference/disequilibrium"
path_baseline="$path_genetic_reference/baseline"
path_weights="$path_genetic_reference/weights"
path_frequencies="$path_genetic_reference/frequencies"

path_gwas="${path_dock}/gwas"
path_heritability="${path_dock}/heritability"
path_genetic_correlation="${path_dock}/genetic_correlation"

###########################################################################
# Execute procedure.

# Define main phenotypes.

studies=()
studies+=("30124842_yengo_2018")
studies+=("30239722_pulit_2018")
studies+=("30482948_walters_2018_all")
studies+=("30482948_walters_2018_eur")
studies+=("30482948_walters_2018_eur_unrel")
studies+=("30718901_howard_2019")
studies+=("29906448_ruderfer_2018_scz_vs_ctl")
studies+=("29906448_ruderfer_2018_scz_bpd_vs_ctl")
studies+=("29906448_ruderfer_2018_scz_vs_bpd")
studies+=("29906448_ruderfer_2018_bpd_vs_ctl")
studies+=("00000000_ripke_2021")
studies+=("31043756_stahl_2019")
studies+=("00000000_mullins_2021_all")
studies+=("00000000_mullins_2021_bpd1")
studies+=("00000000_mullins_2021_bpd2")

# Organize information in format for LDSC.
for study in "${studies[@]}"; do

  # Organize paths.
  path_phenotype_gwas="${path_gwas}/${study}"
  path_gwas_format="${path_phenotype_gwas}/gwas_format.txt"
  path_gwas_format_compress="${path_gwas_format}.gz"
  path_gwas_munge="${path_phenotype_gwas}/gwas_munge"
  path_gwas_munge_suffix="${path_gwas_munge}.sumstats.gz"
  path_gwas_munge_log="${path_gwas_munge}.log"
  #path_phenotype_genetic_correlation="${path_genetic_correlation}/${study}"

  # Munge GWAS summary statistics for use in LDSC.
  $path_ldsc/munge_sumstats.py \
  --sumstats $path_gwas_format_compress \
  --out $path_gwas_munge \
  --merge-alleles $path_alleles/w_hm3.snplist \
  #--a1-inc

  phenotype_study=${study}
  path_phenotype_gwas_munge_suffix=${path_gwas_munge_suffix}
  report="true" # "true" or "false"
  /usr/bin/bash "${path_scripts_record}/4_organize_phenotype_metabolites_correlations.sh" \
  $phenotype_study \
  $path_phenotype_gwas \
  $path_phenotype_gwas_munge_suffix \
  $path_gwas \
  $path_heritability \
  $path_genetic_correlation \
  $path_genetic_reference \
  $path_promiscuity_scripts \
  $path_scripts_record \
  $path_ldsc \
  $report

  # Remove temporary files.
  rm $path_gwas_munge_suffix
  rm $path_gwas_munge_log

done
