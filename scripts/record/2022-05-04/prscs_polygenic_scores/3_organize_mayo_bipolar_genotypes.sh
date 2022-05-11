#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="${path_process}/dock"
path_mayo_bipolar_genotype="${path_dock}/access/mayo_bipolar_genotype"

path_genotype_product_bim_container="${path_dock}/access/mayo_bipolar_genotype/snp_relevance_bim"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_drive_convert_vcf_to_bim="${path_promiscuity_scripts}/utility/drive_convert_directory_all_vcf_to_plink_bim.sh"

###########################################################################
# Organize variables.

# Define patterns for file names.
pattern_genotype_source_vcf_file="MERGED.maf0.dosR20.3.noDups.chr*.dose.vcf.gz" # do not expand with full path yet

name_prefix_file_product_bim="test_"
report="true"

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
set -x

# Initialize directories and batch instances.
rm -r $path_genotype_product_bim_container
mkdir -p $path_genotype_product_bim_container

# Convert information from genotype files in VCF format to BIM format.
/usr/bin/bash "${path_script_drive_convert_vcf_to_bim}" \
$path_mayo_bipolar_genotype \
$pattern_genotype_source_vcf_file \
$path_genotype_product_bim_container \
$name_prefix_file_product_bim \
$path_promiscuity_scripts \
$report
