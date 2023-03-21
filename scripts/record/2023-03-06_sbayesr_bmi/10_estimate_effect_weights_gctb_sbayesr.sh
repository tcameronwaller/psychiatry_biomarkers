#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date, first execution: __ February 2023
# Date, last execution: __ February 2023
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
path_directory_tools=$(<"./waller_tools.txt")
#path_plink2=$(<"./tools_plink2.txt")
path_gctb=$(<"./tools_gctb.txt")
path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_source="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/gwas_vcf_process_2"

path_directory_product_1="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_effects_1"
path_directory_product_2="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_effects_2"
path_directory_product_3="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_effects_3"

path_directory_ld_matrix_1="${path_directory_product_1}/ukbEURu_hm3_shrunk_sparse"
path_directory_ld_matrix_2=$path_directory_product_2
path_directory_ld_matrix_3="${path_directory_product_3}/band_ukb_10k_hm3"

# Files.
path_file_gwas_source="${path_directory_source}/BMI_GIANTUKB_EUR.txt.gz"
path_file_gwas_product="${path_directory_product_1}/BMI_GIANTUKB_EUR.ma"

# sbayesr_1_ld_eur_ukb_50k_hm3_ss
path_file_ld_matrix_1_source="${path_directory_reference}/gctb/ukbEURu_hm3_sparse.zip"
path_file_ld_matrix_1_product="${path_directory_product_1}/ukbEURu_hm3_sparse.zip"
name_file_ld_matrix_1_prefix="ukbEURu_hm3_chr"
name_file_ld_matrix_1_suffix="_v3_50k.ldm.sparse"

# sbayesr_2_ld_eur_ukb_50k_hm3_ss_chi
# This LD matrix is a single file and is not divided by chromosome.
path_file_ld_matrix_2_source="${path_directory_reference}/gctb/ukbEURu_imp_v3_HM3_n50k.chisq10.zip"
path_file_ld_matrix_2_product="${path_directory_product_2}/ukbEURu_imp_v3_HM3_n50k.chisq10.zip"
path_file_base_ld_matrix_2="${path_directory_ld_matrix_2}/ukbEURu_imp_v3_HM3_n50k.chisq10.ldm.sparse"

# sbayesr_3_ld_eur_ukb_10k_hm3_band
path_file_ld_matrix_3_source="${path_directory_reference}/gctb/band_ukb_10k_hm3.zip"
path_file_ld_matrix_3_product="${path_directory_product_3}/band_ukb_10k_hm3.zip"
name_file_ld_matrix_3_prefix="band_chr"
name_file_ld_matrix_3_suffix=".ldm.sparse"

# product
name_file_product_prefix="BMI_GIANTUKB_EUR_"
name_file_product_suffix="_tcw_2023-03-21" # Suffix must not be an empty string.
path_file_base_product_2="${path_directory_product_2}/BMI_GIANTUKB_EUR_tcw_2023-03-21"

# Scripts.
path_script_gwas_format="${path_directory_process}/promiscuity/scripts/gctb/constrain_translate_gwas_standard_to_gctb.sh"
path_script_submit_batch="${path_directory_process}/promiscuity/scripts/gctb/1_submit_batch_gctb_sbayesr_chromosomes.sh"
path_script_batch_run_sbayesr="${path_directory_process}/promiscuity/scripts/gctb/2_run_batch_gctb_sbayesr.sh"
path_script_run_sbayesr="${path_directory_process}/promiscuity/scripts/gctb/run_gctb_sbayesr.sh"

# Initialize directories.
#rm -r $path_directory_product_1
#rm -r $path_directory_product_2
#rm -r $path_directory_product_3
mkdir -p $path_directory_product_1
mkdir -p $path_directory_product_2
mkdir -p $path_directory_product_3
cd $path_directory_product_1



###########################################################################
# Organize parameters.

observations_variant="1"
chromosome_x="false"
threads=1
report="true"


# TODO: TCW; 27 February 2023
# TODO: need a new parameter for each GWAS sum stats <-- "observations_variant"
# TODO: derive from the "fill_observations" parameter in the GWAS format translation parameter table.



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "call_submit_batch_gctb_sbayesr.sh"
  echo "----------"
fi


###########################################################################
# Execute procedure.

##########
# 1. Prepare GWAS summary statistics.

if true; then
  /usr/bin/bash $path_script_gwas_format \
  $path_file_gwas_source \
  $path_file_gwas_product \
  $report
fi



##########
# 2. Prepare LD matrices.

if false; then
  # 1.
  cp $path_file_ld_matrix_1_source $path_file_ld_matrix_1_product
  unzip $path_file_ld_matrix_1_product -d $path_directory_product_1
  # 2.
  cp $path_file_ld_matrix_2_source $path_file_ld_matrix_2_product
  unzip $path_file_ld_matrix_2_product -d $path_directory_product_2
  # 3.
  cp $path_file_ld_matrix_3_source $path_file_ld_matrix_3_product
  unzip $path_file_ld_matrix_3_product -d $path_directory_product_3
fi



##########
# 3. Prepare and submit batch of jobs for processing on each chromosome.

if true; then
  # 1.
  cd $path_directory_product_1
  /usr/bin/bash $path_script_submit_batch \
  $path_file_gwas_product \
  $path_directory_ld_matrix_1 \
  $name_file_ld_matrix_1_prefix \
  $name_file_ld_matrix_1_suffix \
  $path_directory_product_1 \
  $name_file_product_prefix \
  $name_file_product_suffix \
  $observations_variant \
  $chromosome_x \
  $path_script_batch_run_sbayesr \
  $path_script_run_sbayesr \
  $path_gctb \
  $threads \
  $report
  # 3.
  cd $path_directory_product_3
  /usr/bin/bash $path_script_submit_batch \
  $path_file_gwas_product \
  $path_directory_ld_matrix_3 \
  $name_file_ld_matrix_3_prefix \
  $name_file_ld_matrix_3_suffix \
  $path_directory_product_3 \
  $name_file_product_prefix \
  $name_file_product_suffix \
  $observations_variant \
  $chromosome_x \
  $path_script_batch_run_sbayesr \
  $path_script_run_sbayesr \
  $path_gctb \
  $threads \
  $report
fi

if false; then
  # 2.
  cd $path_directory_product_2
  /usr/bin/bash $path_script_run_sbayesr \
  $path_file_gwas_product \
  $path_file_base_ld_matrix_2 \
  $path_file_base_product_2 \
  $observations_variant \
  $path_gctb \
  $threads \
  $report
fi



##########
# 4. Remove temporary files.

# rm -r $path_directory_ld_matrix



#
