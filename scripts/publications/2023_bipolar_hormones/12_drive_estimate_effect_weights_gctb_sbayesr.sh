#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 18 April 2022
# Date, last execution: 18 April 2023
# Review: TCW; 18 April 2023
################################################################################
# Note

# TODO: TCW; 18 April 2023
# TODO: It might or might not be helpful to organize all files for each set of GWAS summary statistics within
# TODO: a separate child directory...



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_tools=$(<"./waller_tools.txt")
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/gwas_extra_process_2023-04-18"
path_directory_product="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/effect_weights_sbayesr_2023-04-18"

# Files.
# Parameters for sets of GWAS summary statistics.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-02-24.tsv"
# Linkage Disequilibrium (LD) correlation matrices.
# sbayesr_1_ld_eur_ukb_50k_hm3_ss
path_file_ld_matrix_1_source="${path_directory_reference}/gctb/ukbEURu_hm3_sparse.zip"
path_file_ld_matrix_1_product="${path_directory_product}/ukbEURu_hm3_sparse.zip"
path_directory_ld_matrix_1="${path_directory_product}/ukbEURu_hm3_shrunk_sparse"
name_file_ld_matrix_1_prefix="ukbEURu_hm3_chr"
name_file_ld_matrix_1_suffix="_v3_50k.ldm.sparse"
# sbayesr_2_ld_eur_ukb_50k_hm3_ss_chi
# This LD matrix is a single file and is not divided by chromosome.
path_file_ld_matrix_2_source="${path_directory_reference}/gctb/ukbEURu_imp_v3_HM3_n50k.chisq10.zip"
path_file_ld_matrix_2_product="${path_directory_product}/ukbEURu_imp_v3_HM3_n50k.chisq10.zip"
path_directory_ld_matrix_2=$path_directory_product
path_file_base_ld_matrix_2="${path_directory_ld_matrix_2}/ukbEURu_imp_v3_HM3_n50k.chisq10.ldm.sparse"
# sbayesr_3_ld_eur_ukb_10k_hm3_band
path_file_ld_matrix_3_source="${path_directory_reference}/gctb/band_ukb_10k_hm3.zip"
path_file_ld_matrix_3_product="${path_directory_product}/band_ukb_10k_hm3.zip"
path_directory_ld_matrix_3="${path_directory_product}/band_ukb_10k_hm3"
name_file_ld_matrix_3_prefix="band_chr"
name_file_ld_matrix_3_suffix=".ldm.sparse"

# Scripts.
path_script_estimate_effects="${path_directory_process}/psychiatric_metabolism/scripts/publications/2023_bipolar_hormones/13_estimate_effect_weights_gctb_sbayesr.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.

suffix_file_gwas_source=".txt.gz"
suffix_file_gwas_product=".ma"
name_file_effect_suffix="_2023-04-18" # Suffix must not be an empty string.
report="true"

################################################################################
# Execute procedure.

##########
# Prepare Linkage Disequilibrium (LD) correlation matrices.
if true; then
  # 1.
  #cp $path_file_ld_matrix_1_source $path_file_ld_matrix_1_product
  #unzip $path_file_ld_matrix_1_product -d $path_directory_product
  # 2.
  #cp $path_file_ld_matrix_2_source $path_file_ld_matrix_2_product
  #unzip $path_file_ld_matrix_2_product -d $path_directory_product
  # 3.
  cp $path_file_ld_matrix_3_source $path_file_ld_matrix_3_product
  unzip $path_file_ld_matrix_3_product -d $path_directory_product
fi

##########
# Iterate on sets of GWAS summary statistics.
# Read lines from file and split fields within each line by space, tab, or new-line delimiters.
input=$path_file_translation
while IFS=$' \t\n' read -r -a array
do
  # Extract variables.
  inclusion="${array[0]}"
  name="${array[2]}"
  fill_observations="${array[10]}"
  # Execute procedure for current record's parameters.
  if [[ $inclusion == "1" ]]; then
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
    # Organize paths and parameters.
    path_file_gwas_source="${path_directory_source}/${name}${suffix_file_gwas_source}"
    path_file_gwas_product="${path_directory_product}/${name}${suffix_file_gwas_product}"
    name_file_effect_prefix="${name}_"
    # Determine whether counts of observations are specific to each variant.
    if [[ "$fill_observations" == "1" ]]; then
      observations_variant="0"
    else
      observations_variant="1"
    fi
    # Confirm whether the source file exists.
    if [ -f "$path_file_gwas_source" ]; then
      # Call script.
      /usr/bin/bash $path_script_estimate_effects \
      $path_file_gwas_source \
      $path_file_gwas_product \
      $path_directory_ld_matrix_3 \
      $name_file_ld_matrix_3_prefix \
      $name_file_ld_matrix_3_suffix \
      $path_directory_product \
      $name_file_effect_prefix \
      $name_file_effect_suffix \
      $observations_variant \
      $report
    fi
  fi
done < "${input}"



##########
# Remove temporary files.

# rm -r $path_directory_ld_matrix_1
# rm -r $path_directory_ld_matrix_2
# rm -r $path_directory_ld_matrix_3



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script completion:"
  echo $0 # Print full file path to script.
  echo "12_drive_estimate_effect_weights_gctb_sbayesr.sh"
  echo "----------"
fi




#
