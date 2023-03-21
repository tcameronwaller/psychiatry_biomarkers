#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 21 March 2023
# Date, last execution: 21 March 2023
# Review: TCW; __ March 2023
################################################################################
# Note



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_scores_ldpred2=$(<"./scores_mayo_bipolar_1_2_merge_ldpred2.txt")
path_directory_scores_ldpred2_2=$(<"./scores_mayo_bipolar_1_2_merge_ldpred2_2.txt")
path_directory_scores_prsice2=$(<"./scores_mayo_bipolar_1_2_merge_prsice2.txt")
path_directory_source_ldpred2="${path_directory_scores_ldpred2}/BMI_GIANTUKB_EUR"
path_directory_source_ldpred2_2="${path_directory_scores_ldpred2_2}/BMI_GIANTUKB_EUR"
path_directory_source_prsice2="${path_directory_scores_prsice2}/BMI_GIANTUKB_EUR"

path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_source_sbayesr_1="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_1_polygenic_scores"
path_directory_source_sbayesr_2="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_2_polygenic_scores"
path_directory_source_sbayesr_3="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_3_polygenic_scores"

path_directory_product="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/polygenic_scores_comparison/BMI_GIANTUKB_EUR_raw"

# Files.
path_file_source_ldpred2="${path_directory_source_ldpred2}/BMI_GIANTUKB_EUR.LDpred2.auto.gz"
path_file_source_ldpred2_2="${path_directory_source_ldpred2_2}/BMI_GIANTUKB_EUR.LDpred2.auto.gz"
path_file_source_prsice2="${path_directory_source_prsice2}/BMI_GIANTUKB_EUR.best.gz"
path_file_source_sbayesr_1="${path_directory_source_sbayesr_1}/BMI_GIANTUKB_EUR_combination.tsv"
path_file_source_sbayesr_2="${path_directory_source_sbayesr_2}/BMI_GIANTUKB_EUR_combination.tsv"
path_file_source_sbayesr_3="${path_directory_source_sbayesr_3}/BMI_GIANTUKB_EUR_combination.tsv"

path_file_product_compress_ldpred2="${path_directory_product}/scores_ldpred2.txt.gz"
path_file_product_compress_ldpred2_2="${path_directory_product}/scores_ldpred2_2.txt.gz"
path_file_product_compress_prsice2="${path_directory_product}/scores_prsice2.txt.gz"
path_file_product_ldpred2="${path_directory_product}/scores_ldpred2.txt"
path_file_product_ldpred2_2="${path_directory_product}/scores_ldpred2_2.txt"
path_file_product_prsice2="${path_directory_product}/scores_prsice2.txt"
path_file_product_sbayesr_1="${path_directory_product}/scores_sbayesr_1.txt"
path_file_product_sbayesr_2="${path_directory_product}/scores_sbayesr_2.txt"
path_file_product_sbayesr_3="${path_directory_product}/scores_sbayesr_3.txt"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

###########################################################################
# Organize parameters.



###########################################################################
# Execute procedure.

# Echo each command to console.
set -x

# Copy raw source files of polygenic scores.
cp $path_file_source_ldpred2 $path_file_product_compress_ldpred2
cp $path_file_source_ldpred2_2 $path_file_product_compress_ldpred2_2
cp $path_file_source_prsice2 $path_file_product_compress_prsice2
cp $path_file_source_sbayesr_1 $path_file_product_sbayesr_1
cp $path_file_source_sbayesr_2 $path_file_product_sbayesr_2
cp $path_file_source_sbayesr_3 $path_file_product_sbayesr_3

# Decompress files of polygenic scores.
gzip -dcvf $path_file_product_compress_ldpred2 > $path_file_product_ldpred2
gzip -dcvf $path_file_product_compress_ldpred2_2 > $path_file_product_ldpred2_2
gzip -dcvf $path_file_product_compress_prsice2 > $path_file_product_prsice2



#
