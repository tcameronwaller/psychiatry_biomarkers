#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 6 August 2023
# Date, last execution: 9 January 2024
# Date, review: 9 January 2023
################################################################################
# Note

# For a reference on a previous script to set up a large count of comparisons
# between primary and secondary sets of GWAS summary statistics, look at the
# script "6_call_submit_gwas_ldsc_genetic_correlation.sh" in the directory
# "/.../sexy_alcohol/repository/scripts/record/2022-08-01/ldsc_heritability_correlation/".

# SLURM batch job: __, __ (group: "primaries"; instances: 6,400; date: 9 January 2024)



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
path_directory_product_parent="${path_directory_group_parent}/6_gwas_correlation_ldsc_primary"
path_directory_product_child="${path_directory_product_parent}/neuropsychiatry_substance_disorders"
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
rm -r $path_directory_product_parent # caution
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
# Primary studies.

# Define array of primary studies.
primaries=()

# Psychiatric and substance-use disorders
# Review: TCW; 9 January 2024
# Count of primary studies and versions: 80

primaries+=("36702997_demontis_2023_adhd")
primaries+=("36477530_saunders_2022_alcohol_all")
primaries+=("36477530_saunders_2022_alcohol_no_ukb")
primaries+=("36477530_saunders_2022_tobacco_all")
primaries+=("36477530_saunders_2022_tobacco_no_ukb")
primaries+=("36477530_saunders_2022_tobacco_ever_all")
primaries+=("36477530_saunders_2022_tobacco_ever_no_ukb")
primaries+=("36477530_saunders_2022_tobacco_age_all")
primaries+=("36477530_saunders_2022_tobacco_age_no_ukb")
primaries+=("36477530_saunders_2022_tobacco_cessation_all")
primaries+=("36477530_saunders_2022_tobacco_cessation_no_ukb")
primaries+=("35396580_trubetskoy_2022_all")
primaries+=("35396580_trubetskoy_2022_female")
primaries+=("35396580_trubetskoy_2022_male")
primaries+=("34099189_blokland_2022_mdd_female")
primaries+=("34099189_blokland_2022_mdd_male")
primaries+=("34099189_blokland_2022_rmdd_female")
primaries+=("34099189_blokland_2022_rmdd_male")
primaries+=("34099189_blokland_2022_scz_female")
primaries+=("34099189_blokland_2022_scz_male")
primaries+=("34099189_blokland_2022_bip_female")
primaries+=("34099189_blokland_2022_bip_male")
primaries+=("34099189_blokland_2022_mdd_sex")
primaries+=("34099189_blokland_2022_rmdd_sex")
primaries+=("34099189_blokland_2022_scz_sex")
primaries+=("34099189_blokland_2022_bip_sex")
primaries+=("34002096_mullins_2021_bd_all")
primaries+=("34002096_mullins_2021_bd_all_alt_1")
primaries+=("34002096_mullins_2021_bd_no_ukb")
primaries+=("34002096_mullins_2021_bd_no_ukb_alt_1")
primaries+=("34002096_mullins_2021_bd_1")
primaries+=("34002096_mullins_2021_bd_1_alt_1")
primaries+=("34002096_mullins_2021_bd_2")
primaries+=("34002096_mullins_2021_bd_2_alt_1")
primaries+=("33096046_johnson_2020_eur_all")
primaries+=("33096046_johnson_2020_eur_unrelated")
primaries+=("32747698_matoba_2020_europe")
primaries+=("32099098_polimanti_2020_eur_opioid_dep_exposed")
primaries+=("32099098_polimanti_2020_eur_opioid_dep_unexposed")
primaries+=("32099098_polimanti_2020_eur_opioid_exposure")
primaries+=("31748690_purves_2020_meta")
primaries+=("31748690_purves_2020_ukb")
primaries+=("31594949_nievergelt_2019_trans_all")
primaries+=("31594949_nievergelt_2019_europe_all")
primaries+=("31308545_watson_2019")
primaries+=("30818990_yu_2019")
primaries+=("30804558_grove_2019")
primaries+=("30718901_howard_2019_pgc_ukb")
primaries+=("30718901_howard_2019_pgc")
primaries+=("30643251_liu_2019_alcohol_all")
primaries+=("30643251_liu_2019_alcohol_no_ukb")
primaries+=("30643251_liu_2019_tobacco_all")
primaries+=("30643251_liu_2019_tobacco_no_ukb")
primaries+=("30643251_liu_2019_tobacco_ever_all")
primaries+=("30643251_liu_2019_tobacco_ever_no_ukb")
primaries+=("30643251_liu_2019_tobacco_age_all")
primaries+=("30643251_liu_2019_tobacco_age_no_ukb")
primaries+=("30643251_liu_2019_tobacco_cessation_all")
primaries+=("30643251_liu_2019_tobacco_cessation_no_ukb")
primaries+=("30482948_walters_2018_eur_all")
primaries+=("30482948_walters_2018_eur_all_alt_1")
primaries+=("30482948_walters_2018_eur_all_alt_2")
primaries+=("30482948_walters_2018_eur_unrel_meta")
primaries+=("30482948_walters_2018_eur_unrel_meta_alt_1")
primaries+=("30482948_walters_2018_eur_unrel_meta_alt_2")
primaries+=("30482948_walters_2018_eur_unrel_genotype")
primaries+=("30482948_walters_2018_eur_unrel_genotype_alt_1")
primaries+=("30482948_walters_2018_eur_unrel_genotype_alt_2")
primaries+=("30482948_walters_2018_female")
primaries+=("30482948_walters_2018_female_alt_1")
primaries+=("30482948_walters_2018_female_alt_2")
primaries+=("30482948_walters_2018_male")
primaries+=("30482948_walters_2018_male_alt_1")
primaries+=("30482948_walters_2018_male_alt_2")
primaries+=("30478444_demontis_2019_adhd")
primaries+=("29700475_wray_2018_pgc_ukb")
primaries+=("29700475_wray_2018_pgc")
primaries+=("29325848_martin_2018_adhd_female")
primaries+=("29325848_martin_2018_adhd_male")
primaries+=("28761083_arnold_2018")




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
  for primary in "${primaries[@]}"; do
    for secondary in "${primaries[@]}"; do
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
  echo "6_estimate_gwas_genetic_correlation_ldsc_primary.sh"
  echo "----------"
fi



#
