#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 27 March 2023
# Date, last execution: 2 October 2023
# Date, review: 2 October 2023
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

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"

path_directory_source="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-29/6_filter_constrain_gwas_values"
path_directory_product="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-29/7_gwas_effective_observations"

# Files.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-09-29.tsv"

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



##########
# Copy the GWAS summary statistics from the previous process.
# Most sets of GWAS summary statistics do not need extra processing.
# Subsequent processes on a few studies will replace the appropriate files.
cp $path_directory_source/*.txt.gz $path_directory_product

# Read lines from file and split fields within each line by space, tab, or new-line delimiters.
input=$path_file_translation
while IFS=$' \t\n' read -r -a array
do
  # Report.
  if [[ "$report" == "true" ]]; then
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
  if [[ "${array[0]}" == "1" ]]; then
    # Define variables.
    name="${array[2]}"
    type="${array[9]}"
    # Organize paths.
    path_file_source="${path_directory_source}/${name}.txt.gz"
    path_file_product="${path_directory_product}/${name}.txt.gz"
    # Determine whether to report SNP heritability on observed or liability
    # scales.
    if [ "$type" == "logistic" ]; then
      # Call script.
      /usr/bin/bash $path_file_script \
      $path_file_source \
      $path_file_product \
      $report
    elif [ "$type" == "linear" ]; then
      # Copy source file to product file.
      #cp $path_file_source $path_file_product
      echo "linear"
    fi
  fi
done < "${input}"

# Extra instance.
path_file_source_extra="${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz"
path_file_product_extra="${path_directory_product}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz"
type_extra="logistic"
# Call script.
/usr/bin/bash $path_file_script \
$path_file_source_extra \
$path_file_product_extra \
$report



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "7_calculate_effective_observations_logistic_gwas.sh"
  echo "----------"
fi



#
