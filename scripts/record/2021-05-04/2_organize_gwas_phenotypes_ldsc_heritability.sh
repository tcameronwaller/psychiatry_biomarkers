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
path_promiscuity_scripts_ldsc_heritability="${path_promiscuity_scripts}/ldsc_genetic_heritability_correlation"
path_scripts_format="${path_promiscuity_scripts}/format_gwas_ldsc"
#path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-05-04"

path_dock="$path_process/dock"
path_genetic_reference="${path_dock}/access/genetic_reference"
path_gwas="${path_dock}/gwas"
path_heritability="${path_dock}/heritability"
path_genetic_correlation="${path_dock}/genetic_correlation"

###########################################################################
# Execute procedure.

# Version check...

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

# Organize information about studies.

# Define multi-dimensional array of cohorts and covariates.
# Review:
# 3 May 2021: TCW reviewed study README's, confirmed source files, and confirmed format scripts

studies=()
studies+=("30124842_yengo_2018;${path_gwas_summaries}/30124842_yengo_2018/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz")
studies+=("30239722_pulit_2018;${path_gwas_summaries}/30239722_pulit_2018/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")
studies+=("30482948_walters_2018_all;${path_gwas_summaries}/30482948_walters_2018/pgc_alcdep.discovery.aug2018_release.txt.gz")
studies+=("30482948_walters_2018_eur;${path_gwas_summaries}/30482948_walters_2018/pgc_alcdep.eur_discovery.aug2018_release.txt.gz")
studies+=("30482948_walters_2018_eur_unrel;${path_gwas_summaries}/30482948_walters_2018/pgc_alcdep.eur_unrelated.aug2018_release.txt.gz")
studies+=("30718901_howard_2019;${path_gwas_summaries}/30718901_howard_2019/PGC_UKB_depression_genome-wide.txt.gz")
studies+=("29906448_ruderfer_2018_scz_vs_ctl;${path_gwas_summaries}/29906448_ruderfer_2018/sczvscont-sumstat.gz")
studies+=("29906448_ruderfer_2018_scz_bpd_vs_ctl;${path_gwas_summaries}/29906448_ruderfer_2018/BDSCZvsCONT.sumstats.gz")
studies+=("29906448_ruderfer_2018_scz_vs_bpd;${path_gwas_summaries}/29906448_ruderfer_2018/SCZvsBD.sumstats.gz")
studies+=("29906448_ruderfer_2018_bpd_vs_ctl;${path_gwas_summaries}/29906448_ruderfer_2018/BDvsCONT.sumstats.gz")
studies+=("00000000_ripke_2021;${path_gwas_summaries}/00000000_ripke_2021/PGC3_SCZ_wave3_public.v2.tsv.gz")
studies+=("31043756_stahl_2019;${path_gwas_summaries}/31043756_stahl_2019/daner_PGC_BIP32b_mds7a_0416a.gz")
studies+=("00000000_mullins_2021_all;${path_gwas_summaries}/00000000_mullins_2021/pgc-bip2021-all.vcf.tsv.gz")
studies+=("00000000_mullins_2021_bpd1;${path_gwas_summaries}/00000000_mullins_2021/pgc-bip2021-BDI.vcf.tsv.gz")
studies+=("00000000_mullins_2021_bpd2;${path_gwas_summaries}/00000000_mullins_2021/pgc-bip2021-BDII.vcf.tsv.gz")

# Organize information in format for LDSC.
for study_details in "${studies[@]}"; do
  IFS=";" read -r -a array <<< "${study_details}"
  study="${array[0]}"
  path_source_file="${array[1]}"
  path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${study}.sh"
  report="true" # "true" or "false"
  /usr/bin/bash "$path_promiscuity_scripts_ldsc_heritability/organize_gwas_phenotype_ldsc_heritability.sh" \
  $study \
  $path_source_file \
  $path_genetic_reference \
  $path_gwas \
  $path_heritability \
  $path_genetic_correlation \
  $path_script_gwas_format \
  $path_promiscuity_scripts \
  $path_ldsc \
  $report
done
