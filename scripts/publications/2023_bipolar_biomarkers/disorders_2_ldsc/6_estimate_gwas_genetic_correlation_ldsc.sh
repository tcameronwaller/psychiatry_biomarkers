#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 6 August 2023
# Date, last execution: 6 August 2023
# Date, review: 6 August 2023
################################################################################
# Note

# For a reference on a previous script to set up a large count of comparisons
# between primary and secondary sets of GWAS summary statistics, look at the
# script "6_call_submit_gwas_ldsc_genetic_correlation.sh" in the directory
# "/.../sexy_alcohol/repository/scripts/record/2022-08-01/ldsc_heritability_correlation/".


################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_reference="${path_directory_dock}/ldsc_gwas_disorders_tcw_2023-08-02/2_reference_ldsc"
path_directory_disequilibrium="${path_directory_reference}/disequilibrium/eur_w_ld_chr"

path_directory_source_primary="${path_directory_dock}/ldsc_gwas_disorders_tcw_2023-08-02/4_gwas_munge_ldsc"
path_directory_source_secondary="${path_directory_dock}/ldsc_gwas_biomarkers_tcw_2023-06-06/4_gwas_munge_ldsc"
path_directory_product_parent="${path_directory_dock}/ldsc_gwas_disorders_tcw_2023-08-02/6_gwas_correlation_ldsc"
path_directory_batch_logs="${path_directory_product_parent}/logs"

# Files.
path_file_batch_instances="${path_directory_product_parent}/batch_instances.txt"
path_file_batch_out="${path_directory_product_parent}/batch_out.txt"
path_file_batch_error="${path_directory_product_parent}/batch_error.txt"

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_directory_ldsc="${path_directory_partner_scripts}/ldsc"
path_file_script_ldsc_correlation="${path_directory_ldsc}/estimate_gwas_genetic_correlation_ldsc.sh"
path_file_script_ldsc_correlation_batch_1="${path_directory_ldsc}/ldsc_correlation_batch_1.sh"

# Initialize directories.
rm -r $path_directory_product_parent # caution
mkdir -p $path_directory_product_parent

# Initialize files.
rm $path_file_batch_instances

################################################################################
# Organize parameters.



##########
# Common parameters.
threads=4
report="true"



##########
# Primary studies.
# Define array of primary studies.
primaries=()
primaries+=("36477530_saunders_2022_alcohol")
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
primaries+=("34002096_mullins_2021_bd_no_ukb")
primaries+=("34002096_mullins_2021_bd_1")
primaries+=("34002096_mullins_2021_bd_2")
primaries+=("30643251_liu_2019_all")
primaries+=("30643251_liu_2019_no_ukb")
primaries+=("30482948_walters_2018_eur_all")
primaries+=("30482948_walters_2018_eur_unrel_meta")
primaries+=("30482948_walters_2018_eur_unrel_genotype")
primaries+=("30482948_walters_2018_female")
primaries+=("30482948_walters_2018_male")

##########
# Secondary studies.
# Define array of secondary studies.
secondaries=()
# Thyroid physiology.
secondaries+=("32769997_zhou_2020_thyroid_hormone")
secondaries+=("30367059_teumer_2018_thyroid_hormone_female")
secondaries+=("30367059_teumer_2018_thyroid_hormone_male")
secondaries+=("36635386_chen_2023_thyroxine_total")
secondaries+=("30367059_teumer_2018_thyroxine_free_all")
secondaries+=("30367059_teumer_2018_thyroxine_free_female")
secondaries+=("30367059_teumer_2018_thyroxine_free_male")
secondaries+=("29875488_sun_2018_thyroid_peroxidase")
secondaries+=("24586183_medici_2014_thyroid_peroxidase_antibody")
secondaries+=("24586183_medici_2014_thyroid_peroxidase_reactivity")
secondaries+=("33441150_dennis_2021_parathyrin")
secondaries+=("36093044_mathieu_2022_hypothyroidism")
secondaries+=("00000000_neale_2020_hyperthyroidism_self")
secondaries+=("32581359_saevarsdottir_2020_thyroid_autoimmunity")
secondaries+=("34594039_sakaue_2021_hashimoto")
secondaries+=("34594039_sakaue_2021_graves")

# Sex hormones.
secondaries+=("32042192_ruth_2020_testosterone_female")
secondaries+=("32042192_ruth_2020_testosterone_male")
secondaries+=("32042192_ruth_2020_testosterone_bioavailable_female")
secondaries+=("32042192_ruth_2020_testosterone_bioavailable_male")
secondaries+=("34255042_schmitz_2021_estradiol_female")
secondaries+=("34255042_schmitz_2021_estradiol_male")
secondaries+=("34822396_pott_2021_progesterone_female")
secondaries+=("34822396_pott_2021_progesterone_male")
secondaries+=("34822396_pott_2021_hydroxyprogesterone_female")
secondaries+=("34822396_pott_2021_hydroxyprogesterone_male")
secondaries+=("34822396_pott_2021_aldosterone_female")
secondaries+=("34822396_pott_2021_aldosterone_male")
secondaries+=("34822396_pott_2021_androstenedione_female")
secondaries+=("34822396_pott_2021_androstenedione_male")
secondaries+=("31169883_pott_2019_dheas_female")
secondaries+=("31169883_pott_2019_dheas_male")
secondaries+=("33441150_dennis_2021_follitropin")
secondaries+=("29875488_sun_2018_follistatin")
secondaries+=("33587031_sinnott-armstrong_2021_lutropin")
secondaries+=("32042192_ruth_2020_shbg_bmi_all")
secondaries+=("32042192_ruth_2020_shbg_bmi_female")
secondaries+=("32042192_ruth_2020_shbg_bmi_male")
secondaries+=("32042192_ruth_2020_shbg_all")
secondaries+=("32042192_ruth_2020_shbg_female")
secondaries+=("32042192_ruth_2020_shbg_male")



# Biomarkers.
secondaries+=("34226706_barton_2021_albumin")
secondaries+=("32059762_manousaki_2020_vitamin_d")
secondaries+=("36635386_chen_2023_cortisol")
secondaries+=("31169883_pott_2019_cortisol_female")
secondaries+=("31169883_pott_2019_cortisol_male")
secondaries+=("35459240_said_2022_c_reactive_protein")
secondaries+=("35078996_gudjonsson_2022_complement_c3")
secondaries+=("35078996_gudjonsson_2022_complement_c4")
secondaries+=("00000000_neale_2020_hemoglobin_glycation")



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
    for secondary in "${secondaries[@]}"; do
      # Organize paths.
      path_directory_product_child="${path_directory_product_parent}/${primary}"
      mkdir -p $path_directory_product_child
      #name_comparison="${primary}_-_${secondary}"
      name_comparison="${secondary}"
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
  echo "6_estimate_gwas_genetic_correlation_ldsc.sh"
  echo "----------"
fi



#
