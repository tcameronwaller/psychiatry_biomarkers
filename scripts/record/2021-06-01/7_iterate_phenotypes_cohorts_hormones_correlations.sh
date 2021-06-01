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
path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-06-01"

path_dock="$path_process/dock"
path_genetic_reference="${path_dock}/access/genetic_reference"
path_gwas="${path_dock}/gwas"
path_gwas_cohorts_models="${path_gwas}/cohorts_models"
path_genetic_correlation="${path_dock}/genetic_correlation"

###########################################################################
# Execute procedure.

file_gwas_cohorts_models_munge_suffix="gwas_munge.sumstats.gz"

if true; then

  # Define main phenotype studies.
  phenotype_studies=()
  phenotype_studies+=("30124842_yengo_2018")
  phenotype_studies+=("30239722_pulit_2018")

  # Organize information in format for LDSC.
  for phenotype_study in "${phenotype_studies[@]}"; do
    # Organize paths.
    path_gwas_phenotype="${path_gwas}/${phenotype_study}"
    path_gwas_phenotype_format_compress="${path_gwas_phenotype}/gwas_format.txt.gz"
    path_gwas_phenotype_munge_suffix="${path_gwas_phenotype}/gwas_munge.sumstats.gz"

    # Organize variables.
    report="true" # "true" or "false"
    /usr/bin/bash "${path_scripts_record}/8_organize_phenotype_cohorts_hormones_correlations.sh" \
    $phenotype_study \
    $path_gwas_phenotype \
    $path_gwas_phenotype_munge_suffix \
    $path_gwas_cohorts_models \
    $file_gwas_cohorts_models_munge_suffix \
    $path_genetic_correlation \
    $path_genetic_reference \
    $path_promiscuity_scripts \
    $path_scripts_record \
    $path_ldsc \
    $report
  done
fi

# TODO: I still want to run these...
# TODO: need to make specific to the BMI BD analyses...


# Define specific pairs for genetic correlation.
if false; then

  # Organize paths.
  path_alleles="$path_genetic_reference/alleles"
  path_disequilibrium="$path_genetic_reference/disequilibrium"
  path_baseline="$path_genetic_reference/baseline"
  path_weights="$path_genetic_reference/weights"
  path_frequencies="$path_genetic_reference/frequencies"

  pairs=()
  pairs+=("female_premenopause_binary_testosterone_log;female_postmenopause_binary_testosterone_log")
  pairs+=("female_premenopause_ordinal_testosterone_log;female_postmenopause_ordinal_testosterone_log")
  pairs+=("female_premenopause_ordinal_testosterone_log;female_perimenopause_ordinal_testosterone_log")
  pairs+=("female_perimenopause_ordinal_testosterone_log;female_postmenopause_ordinal_testosterone_log")

  pairs+=("female_premenopause_ordinal_testosterone_log;male_testosterone_log")
  pairs+=("female_perimenopause_ordinal_testosterone_log;male_testosterone_log")
  pairs+=("female_postmenopause_ordinal_testosterone_log;male_testosterone_log")

  pairs+=("female_premenopause_ordinal_testosterone_free_log;female_postmenopause_ordinal_testosterone_free_log")
  pairs+=("female_premenopause_ordinal_testosterone_free_log;female_perimenopause_ordinal_testosterone_free_log")
  pairs+=("female_perimenopause_ordinal_testosterone_free_log;female_postmenopause_ordinal_testosterone_free_log")

  pairs+=("female_premenopause_ordinal_steroid_globulin_log;female_postmenopause_ordinal_steroid_globulin_log")
  pairs+=("female_premenopause_ordinal_steroid_globulin_log;female_perimenopause_ordinal_steroid_globulin_log")
  pairs+=("female_perimenopause_ordinal_steroid_globulin_log;female_postmenopause_ordinal_steroid_globulin_log")



  for pair in "${pairs[@]}"; do
    # Read information.
    IFS=";" read -r -a array <<< "${pair}"
    study_one="${array[0]}"
    study_two="${array[1]}"
    # Organize paths.
    path_gwas_one_munge_suffix="${path_gwas_cohorts_hormones}/${study_one}/${file_gwas_cohorts_models_munge_suffix}"
    path_gwas_two_munge_suffix="${path_gwas_cohorts_hormones}/${study_two}/${file_gwas_cohorts_models_munge_suffix}"

    # Organize paths.
    path_genetic_correlation_comparison="${path_genetic_correlation}/cohorts_hormones_pairs/${study_one}/${study_two}"
    path_genetic_correlation_report="${path_genetic_correlation_comparison}/correlation"
    path_genetic_correlation_report_suffix="${path_genetic_correlation_report}.log"
    # Initialize directories.
    rm -r $path_genetic_correlation_comparison
    mkdir -p $path_genetic_correlation_comparison
    # Estimate genetic correlation in LDSC.
    $path_ldsc/ldsc.py \
    --rg $path_gwas_one_munge_suffix,$path_gwas_two_munge_suffix \
    --ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
    --w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
    --out $path_genetic_correlation_report
  done
fi
