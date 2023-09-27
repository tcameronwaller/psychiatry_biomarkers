#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 23 December 2022
# Date, last execution: 27 September 2023
# Date, review: 27 September 2023
################################################################################
# Note

# The batch submitted at 11:37 Eastern Time on 26 September 2023 did not include
# the set of GWAS summary statistics below.
# 32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz
# Consider running the GWAS2VCF procedure separately for this set of GWAS
# summary statistics if the priority set for Saevarsdottir 2020 has problems.

################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"

path_directory_source="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-25/4_filter_constrain_gwas_values"
path_directory_product="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-25/5_gwas_clean_gwas2vcf"
#path_directory_product="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-25/5_gwas_clean_gwas2vcf_test"
path_directory_batch="${path_directory_product}/batch"

# Files.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-09-19_biomarkers.tsv"
#path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-09-19_biomarkers_test.tsv"
path_file_batch_instances="${path_directory_batch}/batch_instances.txt"
#path_file_batch_out="${path_directory_product_batch}/batch_out.txt"
#path_file_batch_error="${path_directory_product_batch}/batch_error.txt"

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_file_script_gwas2vcf_batch_1="${path_directory_partner_scripts}/gwas_clean/gwas2vcf_gwas_clean_batch_1.sh"

# Initialize directories.
rm -r $path_directory_product # caution
rm -r $path_directory_batch # caution
mkdir -p $path_directory_product
mkdir -p $path_directory_batch
cd $path_directory_product

# Initialize files.
rm $path_file_batch_instances

################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.

##########
# Organize batch instances.

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
  if [ $raw_inclusion == "1" ]; then
    # Organize paths.
    path_file_gwas_source="${path_directory_source}/${raw_name_study}.txt.gz"
    path_file_gwas_product="${path_directory_product}/${raw_name_study}.txt.gz"
    # Define and append a new batch instance.
    instance="${path_file_gwas_source};${path_file_gwas_product};${raw_type}"
    echo $instance >> $path_file_batch_instances
  fi
done < "${input}"

# Extra instance.
path_file_gwas_source_extra="${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz"
path_file_gwas_product_extra="${path_directory_product}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz"
type_extra="logistic"
# Define and append a new batch instance.
instance_extra="${path_file_gwas_source_extra};${path_file_gwas_product_extra};${type_extra}"
echo $instance_extra >> $path_file_batch_instances

##########
# Batch parallelization.
if true; then
  # Call first script in series for batch execution.
  /usr/bin/bash $path_file_script_gwas2vcf_batch_1 \
  $path_file_batch_instances \
  $path_directory_batch \
  $path_directory_process \
  $report
fi

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "5_clean_gwas_gwas2vcf.sh"
  echo "----------"
fi


#
