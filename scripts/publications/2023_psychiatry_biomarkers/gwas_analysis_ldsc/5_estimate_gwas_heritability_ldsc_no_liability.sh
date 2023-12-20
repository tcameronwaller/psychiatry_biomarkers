#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 8 September 2023
# Date, last execution: 20 December 2023
# Date, review: 20 December 2023
################################################################################
# Note

# This script drives LDSC estimation of SNP heritability on the observed scale
# for all sets of GWAS summary statistics. There is no estimation of SNP
# heritability on the liability scale.

################################################################################
# Organize paths.

# Identifiers or designators of parameter version, preparation batch, and
# analysis batch.
identifier_analysis="gwas_2023-12-19_alcohol_sex_test_ldsc_2023-12-19"
identifier_parameter="tcw_2023-12-19_alcohol_sex_test"

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"

path_directory_group_parent="${path_directory_dock}/${identifier_analysis}"
path_directory_reference="${path_directory_group_parent}/2_reference_ldsc"
path_directory_source="${path_directory_group_parent}/4_gwas_munge_ldsc"
path_directory_product="${path_directory_group_parent}/5_gwas_heritability_ldsc_no_liability"
path_directory_disequilibrium="${path_directory_reference}/disequilibrium/eur_w_ld_chr"

# Files.
path_file_table_parameter="${path_directory_parameters}/table_gwas_translation_${identifier_parameter}.tsv"

# Files.

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_directory_ldsc="${path_directory_partner_scripts}/ldsc"
path_file_script="${path_directory_ldsc}/estimate_gwas_heritability_ldsc.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

################################################################################
# Organize parameters.

threads=4
report="true"

################################################################################
# Execute procedure.


##########
# Read lines from file and split fields within each line by space, tab, or new-line delimiters.
input=$path_file_table_parameter
while IFS=$' \t\n' read -r -a array
do

  # Extract values from individual columns within table's row.
  # These parameters match the format of the table as of 13 December 2023.
  raw_availability="${array[0]}"
  raw_inclusion="${array[1]}"
  raw_directory="${array[2]}"
  raw_name_study="${array[3]}"
  raw_phenotype="${array[4]}"
  raw_sex="${array[5]}"
  raw_name_file_source="${array[6]}"
  raw_suffix_file_source="${array[7]}"
  raw_bgzip="${array[8]}"
  raw_gzip="${array[9]}"
  raw_type="${array[10]}"
  raw_fill_observations="${array[11]}"
  raw_observations_total="${array[12]}"
  raw_fill_case_control="${array[13]}"
  raw_cases="${array[14]}"
  raw_controls="${array[15]}"
  raw_observations_effective="${array[16]}"
  raw_prevalence_sample="${array[17]}"
  raw_prevalence_population="${array[18]}"
  raw_script="${array[19]}"
  raw_note="${array[20]}"

  # Report.
  if [ $raw_inclusion == "1" ] && [ "$report" == "true" ]; then
    echo "----------"
    echo "field 0, availability: ${raw_availability}"
    echo "field 1, inclusion: ${raw_inclusion}"
    echo "field 2, directory: ${raw_directory}"
    echo "field 3, name_study: ${raw_name_study}"
    echo "field 4, phenotype: ${raw_phenotype}"
    echo "field 5, sex: ${raw_sex}"
    echo "field 6, file: ${raw_name_file_source}"
    echo "field 7, suffix: ${raw_suffix_file_source}"
    echo "field 8, bgzip: ${raw_bgzip}"
    echo "field 9, gzip: ${raw_gzip}"
    echo "field 10, type: ${raw_type}"
    echo "field 11, fill_observations: ${raw_fill_observations}"
    echo "field 12, observations_total: ${raw_observations_total}"
    echo "field 13, fill_case_control: ${raw_fill_case_control}"
    echo "field 14, cases: ${raw_cases}"
    echo "field 15, controls: ${raw_controls}"
    echo "field 16, observations_effective: ${raw_observations_effective}"
    echo "field 17, prevalence_sample: ${raw_prevalence_sample}"
    echo "field 18, prevalence_population: ${raw_prevalence_population}"
    echo "field 19, script: ${raw_script}"
    echo "field 20, note: ${raw_note}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  if [ $raw_inclusion == "1" ]; then
    # Organize paths.
    path_file_source="${path_directory_source}/${raw_name_study}.sumstats.gz"
    path_file_base_product="${path_directory_product}/${raw_name_study}"
    # Call LDSC.
    /usr/bin/bash "${path_file_script}" \
    $path_file_source \
    $path_file_base_product \
    $path_directory_disequilibrium \
    $threads \
    $report
  fi
done < "${input}"

if false; then
  # Extra instance.
  path_file_source_extra="${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.sumstats.gz"
  path_file_base_product_extra="${path_directory_product}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute"
  type_extra="logistic"
  # Call script.
  /usr/bin/bash "${path_file_script}" \
  $path_file_source_extra \
  $path_file_base_product_extra \
  $path_directory_disequilibrium \
  $threads \
  $report
fi


################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "5_estimate_gwas_heritability_ldsc_no_liability.sh"
  echo "----------"
fi



#
