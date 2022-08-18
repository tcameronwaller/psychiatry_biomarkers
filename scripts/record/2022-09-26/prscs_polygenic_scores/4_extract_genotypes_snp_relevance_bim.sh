#!/bin/bash

###########################################################################
# Organize paths.
cd ~/paths
path_waller_tools=$(<"./waller_tools.txt")
path_plink2=$(<"./tools_plink2.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_process}/dock"

#host="ucsc"
host="ensembl"
path_directory_genotypes_vcf="${path_directory_dock}/access/mayo_bipolar_biobank_genotypes/grch37_${host}_format_annotation_vcf"
path_directory_genotypes_snp_bim="${path_directory_dock}/genotypes_snp_bim/mayo_bipolar_biobank_grch37_${host}_snp_bim"

# Scripts.
path_directory_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_drive_extract_vcf_to_bim="${path_directory_promiscuity_scripts}/utility/plink/drive_directory_all_extract_vcf_snps_to_plink_bim.sh"

# Initialize directory.
rm -r $path_directory_genotypes_snp_bim
mkdir -p $path_directory_genotypes_snp_bim

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
set -x

################################################################################
# Extract relevant genetic features (SNPs) to from VCF to BIM format.

# UCSC chain: TCW; 18 August 2022; complete
# Ensembl chain: TCW; ___ 2022; ___
if true; then
  # Organize specific paths and parameters.
  pattern_genotype_source_vcf_file="genotype_grch37_chromosome_*.vcf.gz" # do not expand with full path yet
  # Parameters.
  threads=16
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_drive_extract_vcf_to_bim}" \
  $path_directory_genotypes_vcf \
  $pattern_genotype_source_vcf_file \
  $path_directory_genotypes_snp_bim \
  $threads \
  $path_plink2 \
  $path_directory_promiscuity_scripts \
  $report
fi



#
