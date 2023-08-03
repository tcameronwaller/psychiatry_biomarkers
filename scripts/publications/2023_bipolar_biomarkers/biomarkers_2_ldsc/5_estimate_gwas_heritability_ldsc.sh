#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date: 27 December 2022
################################################################################
################################################################################
################################################################################
# Note



################################################################################
################################################################################
################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_reference="${path_directory_dock}/ldsc_gwas_biomarkers_tcw_2023-05-25/2_reference_ldsc"
path_directory_source="${path_directory_dock}/ldsc_gwas_biomarkers_tcw_2023-05-25/4_gwas_munge_ldsc"
path_directory_product="${path_directory_dock}/ldsc_gwas_biomarkers_tcw_2023-05-25/5_gwas_heritability_ldsc"
path_directory_disequilibrium="${path_directory_reference}/disequilibrium/eur_w_ld_chr"

# Files.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-05-23.tsv"

# Files.

# Scripts.
path_directory_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
path_directory_ldsc="${path_directory_promiscuity_scripts}/ldsc"
path_file_script="${path_directory_ldsc}/estimate_gwas_heritability_observed_liability_scale_ldsc.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

################################################################################
# Organize parameters.

threads=8
report="true"

################################################################################
# Execute procedure.


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
    prevalence_sample="${array[15]}"
    prevalence_population="${array[16]}"
    # Organize paths.
    path_file_source="${path_directory_source}/${name}.sumstats.gz"
    path_file_base_product="${path_directory_product}/${name}"
    # Determine whether to report SNP heritability on observed or liability
    # scales.
    if [ "$type" == "logistic" ] && [ "$prevalence_sample" != "NA" ] && [ "$prevalence_population" != "NA" ]; then
      scale="liability"
    else
      scale="observed"
    fi
    # Call LDSC.
    /usr/bin/bash "${path_file_script}" \
    $path_file_source \
    $path_file_base_product \
    $path_directory_disequilibrium \
    $scale \
    $prevalence_sample \
    $prevalence_population
    $threads \
    $report
  fi
done < "${input}"

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "5_estimate_gwas_heritability_observed_liability_ldsc.sh"
  echo "----------"
fi



#
