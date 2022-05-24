#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_bcftools=$(<"./tools_bcftools.txt")
path_plink2=$(<"./tools_plink2.txt")

path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="${path_process}/dock"
path_parameters="${path_dock}/parameters"
path_translations_chromosomes_mayo="${path_parameters}/promiscuity/translations_chromosomes_mayo_bipolar.txt"
path_dbsnp_reference="${path_dock}/access/dbsnp_reference_format/GCF_000001405.39.gz"
path_mayo_bipolar_genotype="${path_dock}/access/mayo_bipolar_genotype"

path_genotype_chromosome="${path_dock}/access/mayo_bipolar_genotype/test_chromosome_identifier"
path_genotype_snp_rsid="${path_dock}/access/mayo_bipolar_genotype/test_snp_identifier"
path_genotype_snp_bim="${path_dock}/access/mayo_bipolar_genotype/test_snp_relevance_bim"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_chromosome_in_vcf="${path_promiscuity_scripts}/utility/bcftools/translate_chromosomes_in_vcf.sh"
path_script_dbsnp_rsid_to_vcf="${path_promiscuity_scripts}/utility/bcftools/introduce_dbsnp_rsid_to_vcf.sh"
path_script_extract_vcf_to_bim="${path_promiscuity_scripts}/utility/extract_vcf_snps_to_plink_bim.sh"
#path_script_drive_convert_vcf_to_bim="${path_promiscuity_scripts}/utility/drive_convert_directory_all_vcf_to_plink_bim.sh"


###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
set -x

################################################################################
# Translate chromosome and base-pair position coordinates from human genome
# assembly GRCh38 to GRCh37.





###########################################################################
# Format genotype information in VCF.

if false; then
  # Initialize directory.
  rm -r $path_genotype_chromosome
  mkdir -p $path_genotype_chromosome
  # Organize specific paths and parameters.
  path_vcf_source="${path_mayo_bipolar_genotype}/MERGED.maf0.dosR20.3.noDups.chr21.dose.vcf.gz"
  path_vcf_product="${path_genotype_chromosome}/genotype_chromosome_21.vcf.gz" # determine suffix from BCFTools argument
  threads=10
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_chromosome_in_vcf}" \
  $path_translations_chromosomes_mayo \
  $path_vcf_source \
  $path_vcf_product \
  $threads \
  $path_bcftools \
  $report
fi

###########################################################################
# Annotate genotype information in VCF.

if false; then
  # Initialize directory.
  rm -r $path_genotype_snp_rsid
  mkdir -p $path_genotype_snp_rsid
  # Organize specific paths and parameters.
  path_vcf_source="${path_genotype_chromosome}/genotype_chromosome_21.vcf.gz"
  path_vcf_product="${path_genotype_snp_rsid}/genotype_chromosome_21.vcf.gz"
  threads=10
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_dbsnp_rsid_to_vcf}" \
  $path_dbsnp_reference \
  $path_vcf_source \
  $path_vcf_product \
  $threads \
  $path_bcftools \
  $report
fi

###########################################################################
# Extract information from VCF to BIM.

if false; then
  # Initialize directory.
  rm -r $path_genotype_snp_bim
  mkdir -p $path_genotype_snp_bim
  # Organize specific paths and parameters.
  path_genotype_source_vcf="${path_genotype_snp_rsid}/genotype_chromosome_21.vcf.gz"
  name_file_product_bim="genotype_chromosome_21" # optional to add prefix or suffix here
  threads=10
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_extract_vcf_to_bim}" \
  $path_genotype_source_vcf \
  $path_genotype_snp_bim \
  $name_file_product_bim \
  $threads \
  $path_plink2 \
  $report
fi

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
