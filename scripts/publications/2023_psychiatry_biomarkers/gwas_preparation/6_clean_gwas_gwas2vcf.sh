#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 23 December 2022
# Date, last execution: 29 November 2023
# Date, review: 27 November 2023
################################################################################
# Note

# SLURM batch job: 2420655 (190 instances; 29 November 2023)

################################################################################



################################################################################
# Organize paths.

# Identifiers or designators of parameter version and preparation batch.
identifier_preparation="tcw_2023-12-14_test"
identifier_parameter="tcw_2023-12-14_test"

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"

path_directory_source="${path_directory_dock}/gwas_preparation_${identifier_preparation}/5_fill_dbsnp_rs_identifiers"
path_directory_product="${path_directory_dock}/gwas_preparation_${identifier_preparation}/6_gwas_clean_gwas2vcf"
path_directory_batch="${path_directory_product}/batch"

# Files.
path_file_table_parameter="${path_directory_parameters}/table_gwas_translation_${identifier_parameter}.tsv"
path_file_batch_instances="${path_directory_batch}/batch_instances.txt"
#path_file_batch_out="${path_directory_batch}/batch_out.txt"
#path_file_batch_error="${path_directory_batch}/batch_error.txt"

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
    path_file_gwas_source="${path_directory_source}/${raw_name_study}.txt.gz"
    path_file_gwas_product="${path_directory_product}/${raw_name_study}.txt.gz"
    # Define and append a new batch instance.
    instance="${path_file_gwas_source};${path_file_gwas_product};${raw_type}"
    echo $instance >> $path_file_batch_instances
  fi
done < "${input}"

# Extra instance.
if false; then
  path_file_gwas_source_extra="${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz"
  path_file_gwas_product_extra="${path_directory_product}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz"
  type_extra="logistic"
  # Define and append a new batch instance.
  instance_extra="${path_file_gwas_source_extra};${path_file_gwas_product_extra};${type_extra}"
  echo $instance_extra >> $path_file_batch_instances
fi

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
  echo "6_clean_gwas_gwas2vcf.sh"
  echo "----------"
fi


#
