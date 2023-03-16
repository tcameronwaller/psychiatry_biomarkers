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
path_bcftools=$(<"./tools_bcftools.txt")
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_dbsnp="${path_directory_reference}/dbsnp/grch38_chromosome/"

path_directory_genotypes_source=$(<"./mayo_bipolar_disorder_genotypes_1_2_merge.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_genotypes_product="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/mayo_bipolar_disorder_genotypes_1_2_merge"
path_directory_genotypes_annotation="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/genotypes_mayo_bpd_1_2_merge_annotation"

# Files.
path_file_dbsnp="${path_directory_dbsnp}/GCF_000001405.40.gz"

name_file_genotypes_prefix="MERGED.maf0.dosR20.3.noDups.chr"
name_file_genotypes_suffix=".dose.vcf.gz"
name_file_genotypes_suffix_index=".dose.vcf.gz.tbi"
chromosome="10"

path_file_genotypes_source="${path_directory_genotypes_source}/${name_file_genotypes_prefix}${chromosome}${name_file_genotypes_suffix}"
path_file_genotypes_product="${path_directory_genotypes_product}/${name_file_genotypes_prefix}${chromosome}${name_file_genotypes_suffix}"
path_file_genotypes_annotation="${path_directory_genotypes_annotation}/${name_file_genotypes_prefix}${chromosome}${name_file_genotypes_suffix}"

# Scripts.
path_script_annotate="${path_directory_process}/promiscuity/scripts/bcftools/annotate_dbsnp_identifiers_vcf.sh"

# Initialize directories.
#rm -r $path_directory_genotypes_product
rm -r $path_directory_genotypes_annotation
#mkdir -p $path_directory_genotypes_product
mkdir -p $path_directory_genotypes_annotation
cd $path_directory_genotypes_product

###########################################################################
# Organize parameters.

threads=4
report="true"


###########################################################################
# Execute procedure.

# Echo each command to console.
set -x

#cp -r "$path_directory_genotypes_source" "$path_directory_genotypes_product" # <-- Use this for all chromosomes
#cp "$path_file_genotypes_source" "$path_file_genotypes_product"
#mv "$path_target_container/MERGED" "$path_target_mayo_bipolar_genotype"


##########
# Introduce dbSNP identifiers to genotype files.

if true; then
  /usr/bin/bash $path_script_annotate \
  $path_file_genotypes_product \
  $path_file_genotypes_annotation \
  $path_file_dbsnp \
  $threads \
  $path_bcftools \
  $report
fi



#
