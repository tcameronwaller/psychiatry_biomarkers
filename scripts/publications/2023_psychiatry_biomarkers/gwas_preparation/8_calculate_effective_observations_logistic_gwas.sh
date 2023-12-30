#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 27 March 2023
# Date, last execution: __ December 2023
# Date, review: 30 December 2023
################################################################################
# Note

# Perform this procedure judiciously.
# It is important not to calculate the effective observations repetitively.
# It is probably inadvisable to calculate effective observations in a collection
# of GWAS summary statistics before storage or before sharing with other users.
# Rather, calculate effective observations as appropriate immediately before
# analyses.



################################################################################
# Organize paths.

# Identifiers or designators of parameter version and preparation batch.
identifier_preparation="gwas_2023-12-30"
identifier_parameter="tcw_2023-12-30_dbsnp_rsid"

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"

#path_directory_source="${path_directory_dock}/${identifier_preparation}/5_fill_dbsnp_rs_identifiers"
path_directory_source="${path_directory_dock}/${identifier_preparation}/7_filter_constrain_gwas_values"
path_directory_product="${path_directory_dock}/${identifier_preparation}/8_gwas_effective_observations"

# Files.
path_file_table_parameter="${path_directory_parameters}/table_gwas_translation_${identifier_parameter}.tsv"

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



##########
# Copy the GWAS summary statistics from the previous process.
# Most sets of GWAS summary statistics do not need extra processing.
# Subsequent processes on a few studies will replace the appropriate files.
cp $path_directory_source/*.txt.gz $path_directory_product

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
    path_file_source="${path_directory_source}/${raw_name_study}.txt.gz"
    path_file_product="${path_directory_product}/${raw_name_study}.txt.gz"
    # Determine whether to report SNP heritability on observed or liability
    # scales.
    if [ "$raw_type" == "logistic" ]; then
      # Call script.
      /usr/bin/bash $path_file_script \
      $path_file_source \
      $path_file_product \
      $report
    elif [ "$raw_type" == "linear" ]; then
      # Copy source file to product file.
      #cp $path_file_source $path_file_product
      echo "linear"
    fi
  fi
done < "${input}"

if false; then
  # Extra instance.
  path_file_source_extra="${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz"
  path_file_product_extra="${path_directory_product}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz"
  type_extra="logistic"
  # Call script.
  /usr/bin/bash $path_file_script \
  $path_file_source_extra \
  $path_file_product_extra \
  $report
fi

##########
# 30482948_walters_2018
# Replace GWAS summary statistics from study "30482948_walters_2018" so as to
# avoid redundant calculation of effective observations.

# Alternate 1 has filled total observations and no calculation of effective
# sample size.

# Alternate 2 preserves the original effective sample size.

rm "${path_directory_product}/30482948_walters_2018_eur_all_alt_1.txt.gz"
cp "${path_directory_source}/30482948_walters_2018_eur_all_alt_1.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_all_alt_1.txt.gz"

rm "${path_directory_product}/30482948_walters_2018_eur_all_alt_2.txt.gz"
cp "${path_directory_source}/30482948_walters_2018_eur_all_alt_2.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_all_alt_2.txt.gz"

rm "${path_directory_product}/30482948_walters_2018_eur_unrel_meta_alt_1.txt.gz"
cp "${path_directory_source}/30482948_walters_2018_eur_unrel_meta_alt_1.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_unrel_meta_alt_1.txt.gz"

rm "${path_directory_product}/30482948_walters_2018_eur_unrel_meta_alt_2.txt.gz"
cp "${path_directory_source}/30482948_walters_2018_eur_unrel_meta_alt_2.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_unrel_meta_alt_2.txt.gz"

rm "${path_directory_product}/30482948_walters_2018_eur_unrel_genotype_alt_1.txt.gz"
cp "${path_directory_source}/30482948_walters_2018_eur_unrel_genotype_alt_1.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_unrel_genotype_alt_1.txt.gz"

rm "${path_directory_product}/30482948_walters_2018_eur_unrel_genotype_alt_2.txt.gz"
cp "${path_directory_source}/30482948_walters_2018_eur_unrel_genotype_alt_2.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_unrel_genotype_alt_2.txt.gz"

rm "${path_directory_product}/30482948_walters_2018_female_alt_1.txt.gz"
cp "${path_directory_source}/30482948_walters_2018_female_alt_1.txt.gz" \
"${path_directory_product}/30482948_walters_2018_female_alt_1.txt.gz"

rm "${path_directory_product}/30482948_walters_2018_female_alt_2.txt.gz"
cp "${path_directory_source}/30482948_walters_2018_female_alt_2.txt.gz" \
"${path_directory_product}/30482948_walters_2018_female_alt_2.txt.gz"

rm "${path_directory_product}/30482948_walters_2018_male_alt_1.txt.gz"
cp "${path_directory_source}/30482948_walters_2018_male_alt_1.txt.gz" \
"${path_directory_product}/30482948_walters_2018_male_alt_1.txt.gz"

rm "${path_directory_product}/30482948_walters_2018_male_alt_2.txt.gz"
cp "${path_directory_source}/30482948_walters_2018_male_alt_2.txt.gz" \
"${path_directory_product}/30482948_walters_2018_male_alt_2.txt.gz"



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "8_calculate_effective_observations_logistic_gwas.sh"
  echo "----------"
fi



#
