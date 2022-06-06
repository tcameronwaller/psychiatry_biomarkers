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
# Human genome reference sequences.
path_human_grch37_sequence_compression="${path_dock}/access/human_genome_sequence/grch37/GRCh37.p13.genome.fa.gz"
path_human_grch37_sequence="${path_dock}/access/human_genome_sequence/grch37/GRCh37.p13.genome.fa"
path_human_grch38_sequence_compression="${path_dock}/access/human_genome_sequence/grch38/GRCh38.p13.genome.fa.gz"
path_human_grch38_sequence="${path_dock}/access/human_genome_sequence/grch38/GRCh38.p13.genome.fa"
# Quality control preparation of genotypes before mapping from genome assembly GRCh38 to GRCh37.
path_mayo_bipolar_genotype_raw="${path_dock}/access/mayo_bipolar_genotype_raw"
path_directory_genotype_preparation_vcf="${path_dock}/genotype_mayo_bipolar/preparation_vcf"
# Mapping from genome assembly GRCh38 to GRCh37.
host="ucsc"
#host="ensembl"
if [[ "$host" == "ucsc" ]]; then
  path_assembly_translation_chain="${path_dock}/access/human_genome_assembly_chain/ucsc/hg38ToHg19.over.chain.gz"
elif [[ "$host" == "ensembl" ]]; then
  path_assembly_translation_chain="${path_dock}/access/human_genome_assembly_chain/ensembl/GRCh38_to_GRCh37.chain.gz"
else
  echo "invalid specification of host for genetic reference information"
fi
path_directory_genotype_assembly="${path_dock}/genotype_mayo_bipolar/assembly_vcf_grch37_${host}"
# Combine and sort records for genetic features across chromosomes.
path_directory_genotype_combination_vcf="${path_dock}/genotype_mayo_bipolar/combination_vcf"
path_file_list_vcf_files_combination="${path_directory_genotype_combination_vcf}/list_files_chromosomes_combination.txt"
path_file_genotype_combination_vcf="${path_directory_genotype_combination_vcf}/genotype_combination.vcf.gz"



#### work in progress... ###

# Split of genetic features by chromosomes.
path_directory_mayo_bipolar_genotype_split="${path_dock}/genotype_mayo_bipolar/assembly_grch37_${host}_split"
# Annotation of genetic features.
path_translations_chromosomes_mayo="${path_parameters}/promiscuity/translations_chromosomes_mayo_bipolar.txt"
path_dbsnp="${path_dock}/access/dbsnp/grch37_format/GCF_000001405.25.gz"
path_directory_mayo_bipolar_genotype_annotation="${path_dock}/genotype_mayo_bipolar/grch37_${host}_annotation"
# Extraction of relevant genetic features to file in BIM format.
path_genotype_snp_relevance_bim="${path_mayo_bipolar_genotype_format}/genotype_snp_relevance_bim"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_prepare_vcf="${path_promiscuity_scripts}/utility/bcftools/1_submit_batch_chromosomes_decompose_align_unique_sort.sh"
path_script_map_genome_assembly="${path_promiscuity_scripts}/utility/crossmap/1_submit_batch_chromosomes_map_genotype_genome_assembly.sh"
path_script_combine_vcf="${path_promiscuity_scripts}/utility/bcftools/1_submit_batch_single_combine_chromosomes_sort_vcf.sh"


#### work in progress... ###
path_script_split_vcf_chromosome="${path_promiscuity_scripts}/utility/bcftools/split_genotype_by_chromosome_vcf.sh"

path_script_submit_genotype_format_annotation="${path_promiscuity_scripts}/utility/bcftools/1_submit_batch_directory_all_vcf_format_annotation.sh"
path_script_drive_extract_vcf_to_bim="${path_promiscuity_scripts}/utility/plink/drive_directory_all_extract_vcf_snps_to_plink_bim.sh"

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
set -x

################################################################################
# Prepare individual genotype files, separate across chromosomes, for mapping
# between genomic assemblies (from GRCh38 to GRCh37).

# batch submission: TCW; at 21:11:31 on 05 June 2022; running
if false; then
  # Initialize directory.
  rm -r $path_directory_genotype_preparation_vcf
  mkdir -p $path_directory_genotype_preparation_vcf
  # Organize specific paths and parameters.
  gzip --decompress --stdout $path_human_grch38_sequence_compression > $path_human_grch38_sequence
  prefix_file_genotype_vcf_source="MERGED.maf0.dosR20.3.noDups.chr" # do not expand with full path yet
  suffix_file_genotype_vcf_source=".dose.vcf.gz" # omit the ".bim" suffix
  chromosome_x="true"
  prefix_file_genotype_vcf_product="genotype_grch38_chromosome_" # do not expand with full path yet
  suffix_file_genotype_vcf_product=".vcf.gz" # omit the ".bim" suffix
  threads=16
  report="true"
  # Call script to prepare genotype files in VCF format.
  /usr/bin/bash "${path_script_prepare_vcf}" \
  $path_mayo_bipolar_genotype_raw \
  $prefix_file_genotype_vcf_source \
  $suffix_file_genotype_vcf_source \
  $chromosome_x \
  $path_directory_genotype_preparation_vcf \
  $prefix_file_genotype_vcf_product \
  $suffix_file_genotype_vcf_product \
  $path_human_grch38_sequence \
  $threads \
  $path_promiscuity_scripts \
  $path_bcftools \
  $report
fi

################################################################################
# Map chromosome and base-pair position coordinates from human genome assembly
# GRCh38 to GRCh37.

# review: TCW; 5 June 2022; I think it's ready here... need to update driver scripts

# UCSC chain: TCW; at 09:01:31 on 06 June 2022; complete
# Ensembl chain: TCW; 09:34:52 on 06 June 2022; running
if false; then
  # Initialize directory.
  rm -r $path_directory_genotype_assembly
  mkdir -p $path_directory_genotype_assembly
  # Organize specific paths and parameters.
  #gzip --decompress --stdout $path_human_grch37_sequence_compression > $path_human_grch37_sequence
  prefix_file_genotype_vcf_source="genotype_grch38_chromosome_" # do not expand with full path yet
  suffix_file_genotype_vcf_source=".vcf.gz" # omit the ".bim" suffix
  chromosome_x="true"
  prefix_file_genotype_vcf_product="genotype_grch37_chromosome_" # do not expand with full path yet
  suffix_file_genotype_vcf_product=".vcf.gz" # omit the ".bim" suffix
  threads=16
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_map_genome_assembly}" \
  $path_directory_genotype_preparation_vcf \
  $prefix_file_genotype_vcf_source \
  $suffix_file_genotype_vcf_source \
  $chromosome_x \
  $path_directory_genotype_assembly \
  $prefix_file_genotype_vcf_product \
  $suffix_file_genotype_vcf_product \
  $path_assembly_translation_chain \
  $path_human_grch37_sequence \
  $threads \
  $path_promiscuity_scripts \
  $path_environment_crossmap \
  $path_bcftools \
  $report
fi

################################################################################
# Combine individual genotype files across chromosomes.
# Sort records for genetic features.
# Split into separate genotype files by chromosome.

# UCSC chain: TCW; at ___ on 06 June 2022; ___
# Ensembl chain: TCW; ___ on 06 June 2022;
if true; then
  # Initialize directory.
  rm -r $path_directory_genotype_combination_vcf
  mkdir -p $path_directory_genotype_combination_vcf
  # Organize specific paths and parameters.
  prefix_file_genotype_vcf_source="genotype_grch37_chromosome_" # do not expand with full path yet
  suffix_file_genotype_vcf_source=".vcf.gz" # omit the ".bim" suffix
  chromosome_x="true"
  threads=32
  report="true"
  # Call script to test organization for combination of VCF files.
  /usr/bin/bash "${path_script_combine_vcf}" \
  $path_directory_genotype_assembly \
  $prefix_file_genotype_vcf_source \
  $suffix_file_genotype_vcf_source \
  $chromosome_x \
  $path_file_list_vcf_files_combination \
  $path_file_genotype_combination_vcf \
  $threads \
  $path_promiscuity_scripts \
  $path_bcftools \
  $report
fi

################################################################################
# Within genotype information, split records for genetic features by chromosome.

# process start: TCW; ___ 2022;
if false; then
  # Initialize directory.
  rm -r $path_directory_mayo_bipolar_genotype_split
  mkdir -p $path_directory_mayo_bipolar_genotype_split
  # Organize specific paths and parameters.
  chromosome_x="true"
  prefix_file_product_vcf="genotype_grch37_chromosome_"
  threads=16
  report="true"
  # Call script to test organization for combination of VCF files.
  /usr/bin/bash "${path_script_split_vcf_chromosome}" \
  $path_file_mayo_genotype_assembly_vcf \
  $chromosome_x \
  $path_directory_mayo_bipolar_genotype_split \
  $prefix_file_product_vcf \
  $threads \
  $path_bcftools \
  $report
fi

###########################################################################
# Format and annotate genotype information in VCF.

# UCSC chain: TCW; ___ 2022;
# Ensembl chain: TCW; ___ 2022;
if false; then
  # Initialize directory.
  rm -r $path_directory_mayo_bipolar_genotype_annotation
  mkdir -p $path_directory_mayo_bipolar_genotype_annotation
  # Organize specific paths and parameters.
  pattern_genotype_source_vcf_file="genotype_grch37_chromosome_*.vcf.gz" # do not expand with full path yet
  threads=16
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_submit_genotype_format_annotation}" \
  $path_directory_mayo_bipolar_genotype_assembly \
  $pattern_genotype_source_vcf_file \
  $path_directory_mayo_bipolar_genotype_annotation \
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
