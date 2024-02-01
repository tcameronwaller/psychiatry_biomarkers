#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 6 August 2023
# Date, last execution: 31 January 2024
# Date, review: 31 January 2024
################################################################################
# Note

# For a reference on a previous script to set up a large count of comparisons
# between primary and secondary sets of GWAS summary statistics, look at the
# script "6_call_submit_gwas_ldsc_genetic_correlation.sh" in the directory
# "/.../sexy_alcohol/repository/scripts/record/2022-08-01/ldsc_heritability_correlation/".

# SLURM batch job: 5419161 (group: "secondaries-thyroid"; instances: 3,844; date: 31 January 2024)



################################################################################
# Organize paths.

# Identifiers or designators of parameter version, preparation batch, and
# analysis batch.
identifier_analysis="gwas_2023-12-30_ldsc_2024-01-08"
identifier_parameter="tcw_2023-12-30_dbsnp_rsid"

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"

path_directory_group_parent="${path_directory_dock}/${identifier_analysis}"
path_directory_reference="${path_directory_group_parent}/2_reference_ldsc"
path_directory_source_primary="${path_directory_group_parent}/4_gwas_munge_ldsc"
path_directory_source_secondary="${path_directory_group_parent}/4_gwas_munge_ldsc"
path_directory_product_parent="${path_directory_group_parent}/6_gwas_correlation_ldsc_secondary"
path_directory_product_child="${path_directory_product_parent}/thyroid_disorders_biomarkers"
path_directory_disequilibrium="${path_directory_reference}/disequilibrium/eur_w_ld_chr"
path_directory_batch="${path_directory_product_parent}/batch"

# Files.
path_file_batch_instances="${path_directory_batch}/batch_instances.txt"
#path_file_batch_out="${path_directory_batch}/batch_out.txt"
#path_file_batch_error="${path_directory_batch}/batch_error.txt"

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_directory_ldsc="${path_directory_partner_scripts}/ldsc"
path_file_script_ldsc_correlation="${path_directory_ldsc}/estimate_gwas_genetic_correlation_ldsc.sh"
path_file_script_ldsc_correlation_batch_1="${path_directory_ldsc}/ldsc_correlation_batch_1.sh"

# Initialize directories.
#rm -r $path_directory_product_parent # caution
rm -r $path_directory_batch # caution
mkdir -p $path_directory_product_parent
mkdir -p $path_directory_product_child
mkdir -p $path_directory_batch

# Initialize files.
rm $path_file_batch_instances

################################################################################
# Organize parameters.



##########
# Common parameters.
threads=2
report="true"



##########
# Secondary studies.

# Define array of secondary studies.
secondaries=()

# Thyroid disorders, thyroid biomarkers.
# Review: TCW; 31 January 2024
# Count of secondary studies and versions: 62

secondaries+=("37872160_williams_2023")
secondaries+=("37872160_williams_2023_dbsnp_rsid")
secondaries+=("36635386_chen_2023_thyroxine_total")
secondaries+=("36635386_chen_2023_thyroxine_total_dbsnp_rsid")
secondaries+=("36093044_mathieu_2022_hypothyroidism")
secondaries+=("34594039_sakaue_2021_multi_hypothyroidism")
secondaries+=("34594039_sakaue_2021_multi_hyperthyroidism")
secondaries+=("34594039_sakaue_2021_multi_hashimoto")
secondaries+=("34594039_sakaue_2021_multi_graves")
secondaries+=("34594039_sakaue_2021_eur_hypothyroidism")
secondaries+=("34594039_sakaue_2021_eur_hyperthyroidism")
secondaries+=("34594039_sakaue_2021_eur_hashimoto")
secondaries+=("34594039_sakaue_2021_eur_graves")
secondaries+=("34594039_sakaue_2021_gc_hypothyroidism")
secondaries+=("34594039_sakaue_2021_gc_hyperthyroidism")
secondaries+=("34594039_sakaue_2021_gc_hashimoto")
secondaries+=("34594039_sakaue_2021_gc_graves")
secondaries+=("34594039_sakaue_2021_multi_hypothyroidism_dbsnp_rsid")
secondaries+=("34594039_sakaue_2021_multi_hyperthyroidism_dbsnp_rsid")
secondaries+=("34594039_sakaue_2021_multi_hashimoto_dbsnp_rsid")
secondaries+=("34594039_sakaue_2021_multi_graves_dbsnp_rsid")
secondaries+=("34594039_sakaue_2021_eur_hypothyroidism_dbsnp_rsid")
secondaries+=("34594039_sakaue_2021_eur_hyperthyroidism_dbsnp_rsid")
secondaries+=("34594039_sakaue_2021_eur_hashimoto_dbsnp_rsid")
secondaries+=("34594039_sakaue_2021_eur_graves_dbsnp_rsid")
secondaries+=("34594039_sakaue_2021_gc_hypothyroidism_dbsnp_rsid")
secondaries+=("34594039_sakaue_2021_gc_hyperthyroidism_dbsnp_rsid")
secondaries+=("34594039_sakaue_2021_gc_hashimoto_dbsnp_rsid")
secondaries+=("34594039_sakaue_2021_gc_graves_dbsnp_rsid")
secondaries+=("33441150_dennis_2021_thyroxine_total")
secondaries+=("33441150_dennis_2021_thyroxine_free")
secondaries+=("33441150_dennis_2021_parathyrin")
secondaries+=("32769997_zhou_2020_thyroid_hormone")
secondaries+=("32769997_zhou_2020_thyroid_hormone_dbsnp_rsid")
secondaries+=("32581359_saevarsdottir_2020_thyroid_autoimmunity")
secondaries+=("32581359_saevarsdottir_2020_thyroid_autoimmunity_dbsnp_rsid")
secondaries+=("30367059_teumer_2018_thyroid_hormone_all")
secondaries+=("30367059_teumer_2018_thyroid_hormone_female")
secondaries+=("30367059_teumer_2018_thyroid_hormone_male")
secondaries+=("30367059_teumer_2018_thyroxine_free_all")
secondaries+=("30367059_teumer_2018_thyroxine_free_female")
secondaries+=("30367059_teumer_2018_thyroxine_free_male")
secondaries+=("30367059_teumer_2018_hypothyroidism")
secondaries+=("30367059_teumer_2018_hyperthyroidism")
secondaries+=("30367059_teumer_2018_thyroid_hormone_all_dbsnp_rsid")
secondaries+=("30367059_teumer_2018_thyroid_hormone_female_dbsnp_rsid")
secondaries+=("30367059_teumer_2018_thyroid_hormone_male_dbsnp_rsid")
secondaries+=("30367059_teumer_2018_thyroxine_free_all_dbsnp_rsid")
secondaries+=("30367059_teumer_2018_thyroxine_free_female_dbsnp_rsid")
secondaries+=("30367059_teumer_2018_thyroxine_free_male_dbsnp_rsid")
secondaries+=("30367059_teumer_2018_hypothyroidism_dbsnp_rsid")
secondaries+=("30367059_teumer_2018_hyperthyroidism_dbsnp_rsid")
secondaries+=("29875488_sun_2018_parathyrin")
secondaries+=("29875488_sun_2018_thyroid_peroxidase")
secondaries+=("29875488_sun_2018_parathyrin_dbsnp_rsid")
secondaries+=("29875488_sun_2018_thyroid_peroxidase_dbsnp_rsid")
secondaries+=("24586183_medici_2014_thyroid_peroxidase_antibody")
secondaries+=("24586183_medici_2014_thyroid_peroxidase_reactivity")
secondaries+=("00000000_neale_2020_hypothyroidism_self")
secondaries+=("00000000_neale_2020_hypothyroidism_icd")
secondaries+=("00000000_neale_2020_hyperthyroidism_self")
secondaries+=("00000000_neale_2020_hyperthyroidism_icd")



##########
# Organize multi-dimensional array of information about comparisons.
# [full path to base name of product file] ; \
# [full path to primary source file of LDSC munge GWAS summary statistics] ; \
# [full path to secondary source file of LDSC munge GWAS summary statistics]

comparisons=()
#comparisons+=(
#  "${path_directory_product}/hypothyroidism_against_hyperthyroidism;\
#  ${path_directory_source_primary}/30367059_teumer_2018_hypothyroidism.sumstats.gz;\
#  ${path_directory_source_secondary}/30367059_teumer_2018_hyperthyroidism.sumstats.gz"
#)

if true; then
  # Assemble array of batch instance details.
  for primary in "${secondaries[@]}"; do
    for secondary in "${secondaries[@]}"; do
      # Organize paths.
      name_comparison="${primary}_-_${secondary}"
      path_file_base_product="${path_directory_product_child}/${name_comparison}"
      path_file_source_primary="${path_directory_source_primary}/${primary}.sumstats.gz"
      path_file_source_secondary="${path_directory_source_secondary}/${secondary}.sumstats.gz"
      # Assemble parameters for comparison.
      comparisons+=("${path_file_base_product};${path_file_source_primary};${path_file_source_secondary}")
    done
  done
fi



################################################################################
# Report.

count_comparisons=${#comparisons[@]}

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Source directory:"
  echo $path_directory_source
  echo "count of comparisons: " $count_comparisons
  echo "first file: " ${comparisons[0]} # notice base-zero indexing
  echo "last file: " ${comparisons[$count_comparisons - 1]}
  echo "----------"
fi

sleep 5s

################################################################################
# Execute procedure.



##########
# Simple iteration.
if false; then
  for comparison in "${comparisons[@]}"; do
    # Separate fields from instance.
    # [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]
    IFS=";" read -r -a array <<< "${comparison}"
    path_file_base_product="${array[0]}"
    path_file_source_primary="${array[1]}"
    path_file_source_secondary="${array[2]}"
    # Estimate Genetic Correlation by LDSC.
    /usr/bin/bash $path_file_script_ldsc_correlation \
    $path_file_source_primary \
    $path_file_source_secondary \
    $path_file_base_product \
    $path_directory_disequilibrium \
    $threads \
    $report
    # Report.
    if [[ "$report" == "true" ]]; then
      echo "----------"
      echo "Script path:"
      echo $path_script
      echo "Product file path:"
      echo $path_file_base_product
      echo "Primary source file path:"
      echo $path_file_source_primary
      echo "Secondary source file path:"
      echo $path_file_source_secondary
      echo "----------"
    fi
  done
fi



##########
# Batch parallelization.
if true; then
  # Organize batch job instances.
  for comparison in "${comparisons[@]}"; do
    # Define parameters in array instance for batch job.
    echo $comparison >> $path_file_batch_instances
  done
  # Call first script in series for batch execution.
  /usr/bin/bash $path_file_script_ldsc_correlation_batch_1 \
  $path_file_batch_instances \
  $path_directory_batch \
  $path_directory_product_parent \
  $path_directory_disequilibrium \
  $path_directory_process \
  $threads \
  $report
fi



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "6_estimate_gwas_genetic_correlation_ldsc_secondary.sh"
  echo "----------"
fi



#
