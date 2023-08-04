#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 27 March 2023
# Date, last execution: 4 August 2023
# Date, review: 4 August 2023
################################################################################
# Note

# Perform this procedure judiciously.
# It is important not to calculate the effective observations repetitively.
# It is probably inadvisable to calculate effective observations in a collection
# of GWAS summary statistics before storage or before sharing with other users.
# Rather, calculate effective observations as appropriate immediately before
# analyses.

# Note: TCW; 3 August 2023

# 30482948_walters_2018
# All sets of GWAS summary statistics from study "30482948_walters_2018" already
# reported estimates of the effective sample size (effective observations).

# Other sets of GWAS summary statistics in the collection of 2 August 2023 did
# not already report estimates of the effective sample size.

# 1. Calculate effective observations for all logistic GWAS on dichotomous
# traits.
# 2. Replace the files from the "30482948_walters_2018" study.


################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"

path_directory_source="${path_directory_dock}/gwas_disorders_tcw_2023-08-02/4_gwas_clean_gwas2vcf"
path_directory_product="${path_directory_dock}/gwas_disorders_tcw_2023-08-02/5_gwas_effective_observations"

# Files.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-08-02.tsv"

# Files.

# Scripts.
path_file_script="${path_directory_process}/partner/scripts/gwas_clean/calculate_effective_observations_logistic_gwas.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

################################################################################
# Organize parameters.

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
    path_file_source="${path_directory_source}/${raw_name_study}.txt.gz"
    path_file_product="${path_directory_product}/${raw_name_study}.txt.gz"
    # Determine whether to report SNP heritability on observed or liability
    # scales.
    if [ $raw_type == "logistic" ]; then
      # Call script.
      /usr/bin/bash $path_file_script \
      $path_file_source \
      $path_file_product \
      $report
    elif [ $raw_type == "linear" ]; then
      # Copy source file to product file.
      cp $path_file_source $path_file_product
    fi
  fi
done < "${input}"



##########
# 30482948_walters_2018
# Replace GWAS summary statistics from study "30482948_walters_2018" so as to
# avoid redundant calculation of effective observations.

cp "${path_directory_source}/30482948_walters_2018_eur_all.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_all.txt.gz"

cp "${path_directory_source}/30482948_walters_2018_eur_unrel_meta.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_unrel_meta.txt.gz"

cp "${path_directory_source}/30482948_walters_2018_eur_unrel_genotype.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_unrel_genotype.txt.gz"

cp "${path_directory_source}/30482948_walters_2018_female.txt.gz" \
"${path_directory_product}/30482948_walters_2018_female.txt.gz"

cp "${path_directory_source}/30482948_walters_2018_male.txt.gz" \
"${path_directory_product}/30482948_walters_2018_male.txt.gz"



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "5_calculate_effective_observations_logistic_gwas.sh"
  echo "----------"
fi



#
