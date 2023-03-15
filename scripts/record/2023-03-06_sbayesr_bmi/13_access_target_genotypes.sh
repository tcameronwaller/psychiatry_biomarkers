#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 15 March 2023
# Date, last execution: 15 March 2023
# Review: TCW; ___
################################################################################
# Note



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_source_genotypes=$(<"./mayo_bipolar_disorder_genotypes_1_2_merge.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_product_genotypes="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/mayo_bipolar_disorder_genotypes_1_2_merge"

# Files.
name_file_genotypes_prefix="MERGED.maf0.dosR20.3.noDups.chr"
name_file_genotypes_suffix=".dose.vcf.gz"
chromosome="10"
path_file_genotypes_source="${path_directory_source_genotypes}/${name_file_genotypes_prefix}${chromosome}${name_file_genotypes_suffix}"
path_file_genotypes_product="${path_directory_product_genotypes}/${name_file_genotypes_prefix}${chromosome}${name_file_genotypes_suffix}"

# Initialize directories.
rm -r $path_directory_product_genotypes
#mkdir -p $path_directory_product_genotypes
cd $path_directory_product_genotypes


###########################################################################
# Execute procedure.

# Echo each command to console.
set -x

cp -r "$path_directory_source_genotypes" "$path_directory_product_genotypes"
#cp "$path_file_genotypes_source" "$path_file_genotypes_product"
#mv "$path_target_container/MERGED" "$path_target_mayo_bipolar_genotype"
