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
path_scripts_record="$path_waller/bipolar_metabolism/scripts/record/2021-04-06"
path_promiscuity_scripts="$path_waller/promiscuity/scripts"
path_promiscuity_scripts_ldsc="$path_promiscuity_scripts/ldsc_genetic_heritability_correlation"
path_scripts_format="$path_waller/promiscuity/scripts/format_gwas_ldsc"

path_dock="$path_waller/dock"
path_genetic_reference="$path_dock/access/genetic_reference"
path_gwas="$path_dock/gwas"
path_heritability="$path_dock/heritability"
path_genetic_correlation="$path_dock/genetic_correlation"

###########################################################################
# Execute procedure.

# Initialize directories.
#rm -r $path_gwas
#rm -r $path_heritability
#rm -r $path_genetic_correlation
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

if true; then
  # Ripke et al, pending, 2021 (PubMed:pending)
  # phenotype: schizophrenia
  study="00000000_ripke_2021"
  source_file="PGC3_SCZ_wave3_public.v2.tsv.gz"
  path_source_directory="${path_gwas_summaries}/${study}"
  path_source_file="${path_source_directory}/${source_file}"
  path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${study}.sh"
  report="true" # "true" or "false"
  /usr/bin/bash "$path_promiscuity_scripts_ldsc/organize_gwas_phenotype_ldsc_heritability.sh" \
  $study \
  $source_file \
  $path_source_file \
  $path_genetic_reference \
  $path_gwas \
  $path_heritability \
  $path_genetic_correlation \
  $path_script_gwas_format \
  $path_promiscuity_scripts \
  $report
fi

if false; then

  # Ruderfer et al, Cell, 2018 (PubMed:29906448)
  # phenotype: bipolar disorder and schizophrenia
  study_source="29906448_ruderfer_2018"
  study="29906448_ruderfer_2018_scz_vs_bpd"
  #source_file="BDSCZvsCONT.sumstats.gz" # "_scz_bpd_vs_ctl"
  source_file="SCZvsBD.sumstats.gz" # "_scz_vs_bpd"
  path_source_directory="${path_gwas_summaries}/${study_source}"
  path_source_file="${path_source_directory}/${source_file}"
  path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${study_source}.sh"
  report="true" # "true" or "false"
  /usr/bin/bash "$path_promiscuity_scripts_ldsc/organize_gwas_phenotype_ldsc_heritability.sh" \
  $study \
  $source_file \
  $path_source_file \
  $path_genetic_reference \
  $path_gwas \
  $path_heritability \
  $path_genetic_correlation \
  $path_script_gwas_format \
  $path_promiscuity_scripts \
  $report

  # PGC3 (pending) et al, (pending), 2021 (PubMed:pending)
  # phenotype: bipolar disorder types I and II
  study_source="00000000_pgc3_2021"
  study="00000000_pgc3_2021_bd2"
  #source_file="pgc-bip2021-all.vcf.tsv.gz" # "_all"
  #source_file="pgc-bip2021-BDI.vcf.tsv.gz" # "_bd1"
  source_file="pgc-bip2021-BDII.vcf.tsv.gz" # "_bd2"
  path_source_directory="${path_gwas_summaries}/${study_source}"
  path_source_file="${path_source_directory}/${source_file}"
  path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${study_source}.sh"
  report="true" # "true" or "false"
  /usr/bin/bash "$path_promiscuity_scripts_ldsc/organize_gwas_phenotype_ldsc_heritability.sh" \
  $study \
  $source_file \
  $path_source_file \
  $path_genetic_reference \
  $path_gwas \
  $path_heritability \
  $path_genetic_correlation \
  $path_script_gwas_format \
  $path_promiscuity_scripts \
  $report

  # Walters et al, Nature Neuroscience, 2018 (PubMed:30482948)
  # phenotype: alcohol dependence
  study="30482948_walters_2018"
  source_file="pgc_alcdep.discovery.aug2018_release.txt.gz"
  path_source_directory="${path_gwas_summaries}/${study}"
  path_source_file="${path_source_directory}/${source_file}"
  path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${study}.sh"
  report="true" # "true" or "false"
  /usr/bin/bash "$path_promiscuity_scripts_ldsc/organize_gwas_phenotype_ldsc_heritability.sh" \
  $study \
  $source_file \
  $path_source_file \
  $path_genetic_reference \
  $path_gwas \
  $path_heritability \
  $path_genetic_correlation \
  $path_script_gwas_format \
  $path_promiscuity_scripts \
  $report

  # Howard et al, Nature Neuroscience, 2019 (PubMed:30718901)
  # phenotype: major depression
  study="30718901_howard_2019"
  source_file="PGC_UKB_depression_genome-wide.txt.gz"
  path_source_directory="${path_gwas_summaries}/${study}"
  path_source_file="${path_source_directory}/${source_file}"
  path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${study}.sh"
  report="true" # "true" or "false"
  /usr/bin/bash "$path_promiscuity_scripts_ldsc/organize_gwas_phenotype_ldsc_heritability.sh" \
  $study \
  $source_file \
  $path_source_file \
  $path_genetic_reference \
  $path_gwas \
  $path_heritability \
  $path_genetic_correlation \
  $path_script_gwas_format \
  $path_promiscuity_scripts \
  $report

  # Stahl et al, Nature Genetics, 2019 (PubMed:31043756)
  # phenotype: bipolar disorder
  study="31043756_stahl_2019"
  source_file="daner_PGC_BIP32b_mds7a_0416a.gz"
  path_source_directory="${path_gwas_summaries}/${study}"
  path_source_file="${path_source_directory}/${source_file}"
  path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${study}.sh"
  report="true" # "true" or "false"
  /usr/bin/bash "$path_promiscuity_scripts_ldsc/organize_gwas_phenotype_ldsc_heritability.sh" \
  $study \
  $source_file \
  $path_source_file \
  $path_genetic_reference \
  $path_gwas \
  $path_heritability \
  $path_genetic_correlation \
  $path_script_gwas_format \
  $path_promiscuity_scripts \
  $report

  # Yengo et al, Human Molecular Genetics, 2018 (PubMed:30124842)
  # phenotype: body mass index
  study="30124842_yengo_2018"
  source_file="Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz"
  path_source_directory="${path_gwas_summaries}/${study}"
  path_source_file="${path_source_directory}/${source_file}"
  path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${study}.sh"
  report="true" # "true" or "false"
  /usr/bin/bash "$path_promiscuity_scripts_ldsc/organize_gwas_phenotype_ldsc_heritability.sh" \
  $study \
  $source_file \
  $path_source_file \
  $path_genetic_reference \
  $path_gwas \
  $path_heritability \
  $path_genetic_correlation \
  $path_script_gwas_format \
  $path_promiscuity_scripts \
  $report

  # Pulit et al, Human Molecular Genetics, 2018 (PubMed:30239722)
  # phenotype: waist-to-hip ratio, adjusted for BMI
  study="30239722_pulit_2018"
  source_file="whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz"
  path_source_directory="${path_gwas_summaries}/${study}"
  path_source_file="${path_source_directory}/${source_file}"
  path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${study}.sh"
  report="true" # "true" or "false"
  /usr/bin/bash "$path_promiscuity_scripts_ldsc/organize_gwas_phenotype_ldsc_heritability.sh" \
  $study \
  $source_file \
  $path_source_file \
  $path_genetic_reference \
  $path_gwas \
  $path_heritability \
  $path_genetic_correlation \
  $path_script_gwas_format \
  $path_promiscuity_scripts \
  $report
fi
