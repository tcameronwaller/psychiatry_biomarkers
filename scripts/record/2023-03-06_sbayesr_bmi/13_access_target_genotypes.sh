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
path_directory_source=$(<"./genotypes_mayo_bipolar_disorder_1_2_merge.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_product="${path_directory_dock}/genotypes_mayo_bipolar_disorder_1_2_merge"

# Files.
name_file_genotypes_prefix="MERGED.maf0.dosR20.3.noDups.chr"
name_file_genotypes_suffix=".dose.vcf.gz"
name_file_genotypes_suffix_index=".dose.vcf.gz.tbi"
#chromosome="10"
#path_file_genotypes_source="${path_directory_source}/${name_file_genotypes_prefix}${chromosome}${name_file_genotypes_suffix}"
#path_file_genotypes_product="${path_directory_product}/${name_file_genotypes_prefix}${chromosome}${name_file_genotypes_suffix}"

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

cp "${path_directory_source}/${name_file_genotypes_prefix}*${name_file_genotypes_suffix}*" $path_directory_product
#cp "$path_file_genotypes_source" "$path_file_genotypes_product"
#mv "$path_target_container/MERGED" "$path_target_mayo_bipolar_genotype"



#
