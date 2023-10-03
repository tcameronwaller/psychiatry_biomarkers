#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 8 September 2023
# Date, last execution: 3 October 2023
# Date, review: 3 October 2023
################################################################################
# Note

# This script drives LDSC estimation of SNP heritability on the observed scale
# for all sets of GWAS summary statistics. There is no estimation of SNP
# heritability on the liability scale.

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_reference="${path_directory_dock}/ldsc_gwas_biomarkers_tcw_2023-09-29/2_reference_ldsc"
path_directory_source="${path_directory_dock}/ldsc_gwas_biomarkers_tcw_2023-09-29/4_gwas_munge_ldsc"
path_directory_product="${path_directory_dock}/ldsc_gwas_biomarkers_tcw_2023-09-29/5_gwas_heritability_ldsc_no_liability"
path_directory_disequilibrium="${path_directory_reference}/disequilibrium/eur_w_ld_chr"

# Files.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-09-29_biomarkers.tsv"

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


# Read lines from file and split fields within each line by space, tab, or new-line delimiters.
input=$path_file_translation
while IFS=$' \t\n' read -r -a array
do

  # Extract values from individual columns within table's row.
  raw_inclusion="${array[0]}"
  raw_directory="${array[1]}"
  raw_name_study="${array[2]}"
  raw_phenotype="${array[3]}"
  raw_sex="${array[4]}"
  raw_name_file_source="${array[5]}"
  raw_suffix_file_source="${array[6]}"
  raw_bgzip="${array[7]}"
  raw_gzip="${array[8]}"
  raw_type="${array[9]}"
  raw_fill_observations="${array[10]}"
  raw_observations="${array[11]}"
  raw_fill_case_control="${array[12]}"
  raw_cases="${array[13]}"
  raw_controls="${array[14]}"
  raw_prevalence_sample="${array[15]}"
  raw_prevalence_population="${array[16]}"
  raw_script="${array[17]}"
  raw_note="${array[18]}"

  # Report.
  if [ $raw_inclusion == "1" ] && [ "$report" == "true" ]; then
    echo "----------"
    echo "field 0, inclusion: ${array[0]}"
    echo "field 1, directory: ${array[1]}"
    echo "field 2, name: ${array[2]}"
    echo "field 3, phenotype: ${array[3]}"
    echo "field 4, sex: ${array[4]}"
    echo "field 5, file: ${array[5]}"
    echo "field 6, suffix: ${array[6]}"
    echo "field 7, bgzip: ${array[7]}"
    echo "field 8, gzip: ${array[8]}"
    echo "field 9, type: ${array[9]}"
    echo "field 10, fill_observations: ${array[10]}"
    echo "field 11, observations: ${array[11]}"
    echo "field 12, fill_case_control: ${array[12]}"
    echo "field 13, cases: ${array[13]}"
    echo "field 14, controls: ${array[14]}"
    echo "field 15, prevalence_sample: ${array[15]}"
    echo "field 16, prevalence_population: ${array[16]}"
    echo "field 17, script: ${array[17]}"
    echo "field 18, note: ${array[18]}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  if [[ $raw_inclusion == "1" ]]; then
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



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "5_estimate_gwas_heritability_ldsc_no_liability.sh"
  echo "----------"
fi



#
