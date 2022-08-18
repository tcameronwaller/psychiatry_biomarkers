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

path_source_gwas_summary_compression="${path_dock}/gwas_format_prscs/34002096_mullins_2021_all_no_mayo/gwas_format.txt.gz"
path_source_gwas_summary="${path_dock}/gwas_format_prscs/34002096_mullins_2021_all_no_mayo/gwas_format.txt"
host="ucsc"
#host="ensembl"
path_mayo_bipolar_genotype_format="${path_dock}/access/mayo_bipolar_biobank_genotypes/grch37_${host}_format_annotation_vcf"
path_genotype_snp_bim_directory="${path_mayo_bipolar_genotype_format}/grch37_${host}_snp_relevance_bim"
path_genetic_reference_prscs="${path_dock}/access/genetic_reference_prscs"
path_product_allele_effect_directory="${path_dock}/prscs_allelic_effects_${host}"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_submit_prscs_effects="${path_promiscuity_scripts}/utility/prscs_polygenic_score/1_submit_batch_chromosomes_prscs_estimate_allelic_effects.sh"

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
set -x

################################################################################
# Estimate posterior allelic effects in PRS-CSX.

# UCSC chain:
# Ensembl chain:
if true; then
  # Initialize directory.
  rm -r $path_product_allele_effect_directory
  mkdir -p $path_product_allele_effect_directory
  # Organize specific paths and parameters.
  gzip --decompress --stdout $path_source_gwas_summary_compression > $path_source_gwas_summary
  count_gwas_samples=411778
  # file pattern: "genotype_grch37_chromosome_*.vcf.gz"
  prefix_genotype_snp_bim_file="genotype_grch37_chromosome_" # do not expand with full path yet
  suffix_genotype_snp_bim_file=".vcf.gz" # omit the ".bim" suffix
  population_ancestry="EUR"
  name_file_product_prefix="bpd"
  chromosome_x="false"
  threads=16
  report="true"
  # Estimate posterior allelic effects in PRS-CSX.
  /usr/bin/bash "${path_script_submit_prscs_effects}" \
  $path_source_gwas_summary \
  $count_gwas_samples \
  $path_genotype_snp_bim_directory \
  $prefix_genotype_snp_bim_file \
  $suffix_genotype_snp_bim_file \
  $path_genetic_reference_prscs \
  $population_ancestry \
  $path_product_allele_effect_directory \
  $name_file_product_prefix \
  $chromosome_x \
  $threads \
  $path_promiscuity_scripts \
  $path_environment_prscs \
  $path_prscsx \
  $report
fi



#
