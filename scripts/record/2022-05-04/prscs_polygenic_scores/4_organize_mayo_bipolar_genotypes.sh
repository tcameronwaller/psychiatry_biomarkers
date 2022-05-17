#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_bcftools=$(<"./tools_bcftools.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="${path_process}/dock"
path_dbsnp_reference="${path_dock}/access/dbsnp_reference/GCF_000001405.39.gz"
path_mayo_bipolar_genotype="${path_dock}/access/mayo_bipolar_genotype"
path_genotype_dbsnp_annotation="${path_dock}/access/mayo_bipolar_genotype/dbsnp_annotation"
path_genotype_product_bim_container="${path_dock}/access/mayo_bipolar_genotype/snp_relevance_bim"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"

path_scripts_utility="${path_promiscuity_scripts}/utility"
path_script_dbsnp_rsid_to_vcf="${path_scripts_utility}/introduce_dbsnp_rsid_to_vcf.sh"
path_script_drive_convert_vcf_to_bim="${path_scripts_utility}/drive_convert_directory_all_vcf_to_plink_bim.sh"


###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
set -x

# TODO: TCW, 16 May 2022
# TODO: Before translating VCF to BIM...
# TODO: Introduce rsIDs to the VCF files

# TODO: temporarily run test on a single genotype file
# TODO: then write a "driver" script like I did for the VCF to BIM extraction


###########################################################################
# Annotate information in VCF.

# Initialize directory.
rm -r $path_genotype_dbsnp_annotation
mkdir -p $path_genotype_dbsnp_annotation

path_vcf_source="${path_mayo_bipolar_genotype}/MERGED.maf0.dosR20.3.noDups.chr7.dose.vcf.gz"
path_vcf_product="${path_genotype_dbsnp_annotation}/genotype_annotation_chromosome_7.vcf.gz"

report="true"

# Convert information from genotype files in VCF format to BIM format.
/usr/bin/bash "${path_script_dbsnp_rsid_to_vcf}" \
$path_dbsnp_reference \
$path_vcf_source \
$path_vcf_product \
$path_bcftools \
$report

###########################################################################
# Extract information from VCF to BIM.

# Define patterns for file names.
pattern_genotype_source_vcf_file="MERGED.maf0.dosR20.3.noDups.chr*.dose.vcf.gz" # do not expand with full path yet
report="true"


if false; then
  # Initialize directory.
  rm -r $path_genotype_product_bim_container
  mkdir -p $path_genotype_product_bim_container

  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_drive_convert_vcf_to_bim}" \
  $path_mayo_bipolar_genotype \
  $pattern_genotype_source_vcf_file \
  $path_genotype_product_bim_container \
  $path_promiscuity_scripts \
  $report
fi
