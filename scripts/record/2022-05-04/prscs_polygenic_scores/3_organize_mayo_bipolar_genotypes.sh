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
path_human_genome_sequence="${path_dock}/access/human_genome_sequence/grch37/GRCh37.p13.genome.fa.gz"
#path_human_genome_assembly_chain="${path_dock}/access/human_genome_assembly_chain/ucsc/hg38ToHg19.over.chain.gz"
path_human_genome_assembly_chain="${path_dock}/access/human_genome_assembly_chain/ensembl/GRCh38_to_GRCh37.chain.gz"
path_dbsnp="${path_dock}/access/dbsnp/grch37_format/GCF_000001405.25.gz"
path_mayo_bipolar_genotype_raw="${path_dock}/access/mayo_bipolar_genotype_raw"
path_mayo_bipolar_genotype_assembly="${path_dock}/access/mayo_bipolar_genotype_grch37"
path_mayo_bipolar_genotype_format="${path_dock}/access/mayo_bipolar_genotype_format"
path_genotype_snp_relevance_bim="${path_dock}/access/mayo_bipolar_genotype_format/genotype_snp_relevance_bim"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_submit_genotype_translate_assembly="${path_promiscuity_scripts}/utility/crossmap/1_submit_batch_directory_all_vcf_assembly_grch38_to_grch37.sh"
path_script_submit_genotype_format_annotation="${path_promiscuity_scripts}/utility/bcftools/1_submit_batch_directory_all_vcf_format_annotation.sh"
path_script_drive_extract_vcf_to_bim="${path_promiscuity_scripts}/utility/drive_directory_all_extract_vcf_snps_to_plink_bim.sh"

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
set -x

################################################################################
# Translate chromosome and base-pair position coordinates from human genome
# assembly GRCh38 to GRCh37.

if true; then
  # Initialize directory.
  rm -r $path_mayo_bipolar_genotype_assembly
  mkdir -p $path_mayo_bipolar_genotype_assembly
  # Organize specific paths and parameters.
  pattern_genotype_source_vcf_file="MERGED.maf0.dosR20.3.noDups.chr*.dose.vcf.gz" # do not expand with full path yet
  threads=16
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_submit_genotype_translate_assembly}" \
  $path_mayo_bipolar_genotype_raw \
  $pattern_genotype_source_vcf_file \
  $path_mayo_bipolar_genotype_format \
  $path_translations_chromosomes_mayo \
  $path_dbsnp_reference \
  $path_promiscuity_scripts \
  $threads \
  $path_bcftools \
  $report
fi

###########################################################################
# Format and annotate genotype information in VCF.

if false; then
  # Initialize directory.
  rm -r $path_mayo_bipolar_genotype_format
  mkdir -p $path_mayo_bipolar_genotype_format
  # Organize specific paths and parameters.
  pattern_genotype_source_vcf_file="MERGED.maf0.dosR20.3.noDups.chr*.dose.vcf.gz" # do not expand with full path yet
  threads=16
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_submit_genotype_format_annotation}" \
  $path_mayo_bipolar_genotype_assembly \
  $pattern_genotype_source_vcf_file \
  $path_mayo_bipolar_genotype_format \
  $path_translations_chromosomes_mayo \
  $path_dbsnp \
  $path_promiscuity_scripts \
  $threads \
  $path_bcftools \
  $report
fi

###########################################################################
# Extract information from VCF to BIM.

if false; then
  # Initialize directory.
  rm -r $path_genotype_snp_relevance_bim
  mkdir -p $path_genotype_snp_relevance_bim
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
  pattern_genotype_source_vcf_file="snp_MERGED.maf0.dosR20.3.noDups.chr*.dose.vcf.gz" # do not expand with full path yet
  report="true"

  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_drive_convert_vcf_to_bim}" \
  $path_mayo_bipolar_genotype_format \
  $pattern_genotype_source_vcf_file \
  $path_genotype_snp_relevance_bim \
  $path_promiscuity_scripts \
  $report
fi
