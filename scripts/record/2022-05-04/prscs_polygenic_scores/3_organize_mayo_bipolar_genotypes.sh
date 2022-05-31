#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_waller_tools=$(<"./waller_tools.txt")
path_environment_crossmap="${path_waller_tools}/python/environments/crossmap"
path_bcftools=$(<"./tools_bcftools.txt")
path_plink2=$(<"./tools_plink2.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="${path_process}/dock"
path_parameters="${path_dock}/parameters"

path_directory_genotype_preparation_bcf="${path_dock}/genotype_mayo_bipolar/preparation_bcf"
path_file_list_vcf_bcf_files_combination="${path_directory_genotype_preparation_bcf}/list_files_chromosomes_combination.txt"
path_directory_genotype_combination_vcf="${path_dock}/genotype_mayo_bipolar/combination_vcf"
path_file_genotype_combination_vcf="${path_directory_genotype_combination_vcf}/genotype_combination.vcf.gz"

path_translations_chromosomes_mayo="${path_parameters}/promiscuity/translations_chromosomes_mayo_bipolar.txt"
path_human_genome_sequence_compress="${path_dock}/access/human_genome_sequence/grch37/GRCh37.p13.genome.fa.gz"
path_human_genome_sequence="${path_dock}/access/human_genome_sequence/grch37/GRCh37.p13.genome.fa"
#host="ucsc"
host="ensembl"
if [[ "$host" == "ucsc" ]]; then
  path_assembly_translation_chain="${path_dock}/access/human_genome_assembly_chain/ucsc/hg38ToHg19.over.chain.gz"
elif [[ "$host" == "ensembl" ]]; then
  path_assembly_translation_chain="${path_dock}/access/human_genome_assembly_chain/ensembl/GRCh38_to_GRCh37.chain.gz"
else
  echo "invalid specification of host for genetic reference information"
fi
path_dbsnp="${path_dock}/access/dbsnp/grch37_format/GCF_000001405.25.gz"
path_mayo_bipolar_genotype_raw="${path_dock}/access/mayo_bipolar_genotype_raw"
path_mayo_bipolar_genotype_assembly="${path_dock}/access/mayo_bipolar_genotype_grch37_${host}"
path_mayo_bipolar_genotype_format="${path_dock}/access/mayo_bipolar_genotype_grch37_${host}_format"
path_genotype_snp_relevance_bim="${path_mayo_bipolar_genotype_format}/genotype_snp_relevance_bim"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_prepare_combine_multiple_vcf="${path_promiscuity_scripts}/utility/bcftools/1_submit_batch_chromosomes_prepare_vcf_bcf_for_combination.sh"
path_script_combine_bcf_convert_vcf="${path_promiscuity_scripts}/utility/bcftools/combine_bcf_sort_convert_vcf.sh"

path_script_submit_genotype_translate_assembly="${path_promiscuity_scripts}/utility/crossmap/1_submit_batch_directory_all_vcf_assembly_grch38_to_grch37.sh"
path_script_submit_genotype_format_annotation="${path_promiscuity_scripts}/utility/bcftools/1_submit_batch_directory_all_vcf_format_annotation.sh"
path_script_drive_extract_vcf_to_bim="${path_promiscuity_scripts}/utility/plink/drive_directory_all_extract_vcf_snps_to_plink_bim.sh"

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
set -x

################################################################################
# Prepare individual genotype files, separate across chromosomes, for
# combination prior to mapping between genomic assemblies.

# batch submission: TCW; 27 May 2022
if false; then
  # Initialize directory.
  rm -r $path_directory_genotype_preparation_bcf
  mkdir -p $path_directory_genotype_preparation_bcf
  # Organize specific paths and parameters.
  prefix_file_source_genotype_vcf="MERGED.maf0.dosR20.3.noDups.chr" # do not expand with full path yet
  suffix_file_source_genotype_vcf=".dose.vcf.gz" # omit the ".bim" suffix
  chromosome_x="true"
  threads=16
  report="true"
  # Call script to test organization for combination of VCF files.
  /usr/bin/bash "${path_script_prepare_combine_multiple_vcf}" \
  $path_mayo_bipolar_genotype_raw \
  $prefix_file_source_genotype_vcf \
  $suffix_file_source_genotype_vcf \
  $chromosome_x \
  $path_directory_genotype_preparation_bcf \
  $path_file_list_vcf_bcf_files_combination \
  $threads \
  $path_promiscuity_scripts \
  $path_bcftools \
  $report
fi

################################################################################
# Combine individual genotype files across chromosomes.
# Write to genotype file in VCF format with BGZip compression.
# This combination genotype file is for mapping from GRCh38 to GRCh37.

if true; then
  # Initialize directory.
  rm -r $path_directory_genotype_combination_vcf
  mkdir -p $path_directory_genotype_combination_vcf
  # Organize specific paths and parameters.
  threads=16
  report="true"
  # Call script to test organization for combination of VCF files.
  /usr/bin/bash "${path_script_combine_bcf_convert_vcf}" \
  $path_file_list_vcf_bcf_files_combination \
  $path_file_genotype_combination_vcf \
  $threads \
  $path_bcftools \
  $report
fi

################################################################################
# Translate chromosome and base-pair position coordinates from human genome
# assembly GRCh38 to GRCh37.

# UCSC chain: TCW; 25 May 2022; running
# Ensembl chain: TCW; 25 May 2022; running
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

###########################################################################
# Format and annotate genotype information in VCF.

# UCSC chain: TCW; 25 May 2022; complete
# Ensembl chain: TCW; 25 May 2022; complete
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

# UCSC chain: TCW; 25 May 2022; complete
# Ensembl chain: TCW; 25 May 2022; complete
if false; then
  # Initialize directory.
  rm -r $path_genotype_snp_relevance_bim
  mkdir -p $path_genotype_snp_relevance_bim
  # Organize specific paths and parameters.
  pattern_genotype_source_vcf_file="snp_MERGED.maf0.dosR20.3.noDups.chr*.dose.vcf.gz" # do not expand with full path yet
  threads=16
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_drive_extract_vcf_to_bim}" \
  $path_mayo_bipolar_genotype_format \
  $pattern_genotype_source_vcf_file \
  $path_genotype_snp_relevance_bim \
  $threads \
  $path_plink2 \
  $path_promiscuity_scripts \
  $report
fi
