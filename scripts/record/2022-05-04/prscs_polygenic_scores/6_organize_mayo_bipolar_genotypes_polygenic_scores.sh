#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_waller_tools=$(<"./waller_tools.txt")
path_environment_prscs="${path_waller_tools}/python/environments/prs_cs"
path_prscsx="${path_waller_tools}/prs_cs/PRScsx/PRScsx.py"
path_plink2=$(<"./tools_plink2.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="${path_process}/dock"
host="ucsc"
#host="ensembl"
path_mayo_bipolar_genotype_format="${path_dock}/access/mayo_bipolar_genotype_grch37_${host}_format"
path_genotype_snp_bim_directory="${path_mayo_bipolar_genotype_format}/genotype_snp_relevance_bim"

path_source_gwas_summary_compression="${path_dock}/gwas_format_prscs/34002096_mullins_2021_all_no_mayo/gwas_format.txt.gz"
path_source_gwas_summary="${path_dock}/gwas_format_prscs/34002096_mullins_2021_all_no_mayo/gwas_format.txt"
# file pattern: "snp_MERGED.maf0.dosR20.3.noDups.chr*.dose.vcf.gz"
prefix_genotype_snp_bim_file="snp_MERGED.maf0.dosR20.3.noDups.chr" # do not expand with full path yet
suffix_genotype_snp_bim_file=".dose.vcf.gz" # omit the ".bim" suffix

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_submit_genotype_translate_assembly="${path_promiscuity_scripts}/utility/crossmap/1_submit_batch_directory_all_vcf_assembly_grch38_to_grch37.sh"
path_script_submit_genotype_format_annotation="${path_promiscuity_scripts}/utility/bcftools/1_submit_batch_directory_all_vcf_format_annotation.sh"
path_script_drive_extract_vcf_to_bim="${path_promiscuity_scripts}/utility/plink/drive_directory_all_extract_vcf_snps_to_plink_bim.sh"

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
set -x

################################################################################
# Translate chromosome and base-pair position coordinates from human genome
# assembly GRCh38 to GRCh37.

# UCSC chain:
# Ensembl chain:
if false; then
  # Initialize directory.
  rm -r $path_mayo_bipolar_genotype_assembly
  mkdir -p $path_mayo_bipolar_genotype_assembly
  # Organize specific paths and parameters.
  #gzip --decompress --stdout $path_human_genome_sequence_compress > $path_human_genome_sequence
  pattern_genotype_source_vcf_file="MERGED.maf0.dosR20.3.noDups.chr*.dose.vcf.gz" # do not expand with full path yet
  threads=16
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_submit_genotype_translate_assembly}" \
  $path_mayo_bipolar_genotype_raw \
  $pattern_genotype_source_vcf_file \
  $path_mayo_bipolar_genotype_assembly \
  $path_assembly_translation_chain \
  $path_human_genome_sequence \
  $threads \
  $path_promiscuity_scripts \
  $path_environment_crossmap \
  $path_bcftools \
  $report
fi
