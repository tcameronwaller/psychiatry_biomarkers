#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_bcftools=$(<"./tools_bcftools.txt")
path_plink2=$(<"./tools_plink2.txt")

path_process=$(<"./process_psychiatric_metabolism.txt")
path_parameters="$path_process/dock/parameters"
path_translations_chromosome_prefix="$path_parameters/promiscuity/translations_chromosomes_chr_prefix.txt"
path_dock="${path_process}/dock"
path_dbsnp_reference="${path_dock}/access/dbsnp_reference/GCF_000001405.39.gz"
path_mayo_bipolar_genotype="${path_dock}/access/mayo_bipolar_genotype"

path_genotype_chr_format="${path_dock}/access/mayo_bipolar_genotype/chr_format"
path_genotype_dbsnp_annotation="${path_dock}/access/mayo_bipolar_genotype/dbsnp_annotation"
path_genotype_snp_bim="${path_dock}/access/mayo_bipolar_genotype/snp_relevance_bim"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"

path_scripts_utility="${path_promiscuity_scripts}/utility"
path_script_chr_prefix_in_vcf="${path_scripts_utility}/remove_chr_prefix_in_vcf.sh"
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
# Format genotype information in VCF.

# Initialize directory.
rm -r $path_genotype_chr_format
mkdir -p $path_genotype_chr_format

path_vcf_source="${path_mayo_bipolar_genotype}/MERGED.maf0.dosR20.3.noDups.chr21.dose.vcf.gz"
path_vcf_product="${path_genotype_chr_format}/genotype_chromosome_21.vcf.gz"
report="true"

# Convert information from genotype files in VCF format to BIM format.
/usr/bin/bash "${path_script_chr_prefix_in_vcf}" \
$path_translations_chromosome_prefix \
$path_vcf_source \
$path_vcf_product \
$path_bcftools \
$report \
$path_dbsnp_reference

###########################################################################
# Annotate genotype information in VCF.

if false; then
  # Initialize directory.
  rm -r $path_genotype_dbsnp_annotation
  mkdir -p $path_genotype_dbsnp_annotation

  path_vcf_source="${path_mayo_bipolar_genotype}/MERGED.maf0.dosR20.3.noDups.chr21.dose.vcf.gz"
  path_vcf_product="${path_genotype_dbsnp_annotation}/genotype_chromosome_21.vcf.gz"

  report="true"

  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_dbsnp_rsid_to_vcf}" \
  $path_dbsnp_reference \
  $path_vcf_source \
  $path_vcf_product \
  $path_bcftools \
  $report
fi

###########################################################################
# Extract information from VCF to BIM.

# Initialize directory.
rm -r $path_genotype_snp_bim
mkdir -p $path_genotype_snp_bim

# Temporary Testing...

# Scripts.
path_script_extract_vcf_to_bim="${path_scripts_utility}/extract_vcf_snps_to_plink_bim.sh"
path_genotype_source_vcf="${path_genotype_chr_format}/genotype_chromosome_21.vcf.gz"
#path_genotype_source_vcf="${path_genotype_dbsnp_annotation}/genotype_chromosome_21.vcf.gz"
name_file_product_bim="genotype_chromosome_21" # optional to add prefix or suffix here
report="true"

/usr/bin/bash "${path_script_extract_vcf_to_bim}" \
$path_genotype_source_vcf \
$path_genotype_snp_bim \
$name_file_product_bim \
$path_plink2 \
$report


if false; then

  # Define patterns for file names.
  pattern_genotype_source_vcf_file="MERGED.maf0.dosR20.3.noDups.chr*.dose.vcf.gz" # do not expand with full path yet
  report="true"

  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_drive_convert_vcf_to_bim}" \
  $path_mayo_bipolar_genotype \
  $pattern_genotype_source_vcf_file \
  $path_genotype_product_bim_container \
  $path_promiscuity_scripts \
  $report
fi
