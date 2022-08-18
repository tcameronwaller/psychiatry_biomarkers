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
path_directory_dock="${path_process}/dock"

path_directory_genetic_reference_prscs="${path_directory_dock}/access/genetic_reference_prscs"
host="ucsc"
#host="ensembl"
path_directory_genotypes_snp_bim="${path_directory_dock}/genotypes_snp_bim/mayo_bipolar_biobank_grch37_${host}_snp_bim"
path_file_gwas_summary_compression="${path_directory_dock}/gwas_format_prscs/34002096_mullins_2021_all_no_mayo/gwas_format.txt.gz"
path_file_gwas_summary="${path_directory_dock}/gwas_format_prscs/34002096_mullins_2021_all_no_mayo/gwas_format.txt"
path_directory_allele_effect="${path_directory_dock}/prscs_allelic_effects_bipolar_disorder"

# Scripts.
path_directory_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_submit_prscs_effects="${path_directory_promiscuity_scripts}/utility/prscs_polygenic_score/1_submit_batch_chromosomes_prscs_estimate_allelic_effects.sh"

# Initialize directory.
rm -r $path_directory_allele_effect
mkdir -p $path_directory_allele_effect

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
  # Organize specific paths and parameters.
  gzip --decompress --stdout $path_file_gwas_summary_compression > $path_file_gwas_summary
  count_gwas_samples=411778
  # file pattern: "genotype_grch37_chromosome_*.vcf.gz"
  prefix_genotype_snp_bim_file="genotype_grch37_chromosome_" # do not expand with full path yet
  suffix_genotype_snp_bim_file=".vcf.gz" # omit the ".bim" suffix
  population_ancestry="EUR"
  name_file_product_prefix="bpd"
  chromosome_x="false"
  threads=1
  report="true"
  # Estimate posterior allelic effects in PRS-CSX.
  /usr/bin/bash "${path_script_submit_prscs_effects}" \
  $path_file_gwas_summary \
  $count_gwas_samples \
  $path_directory_genotypes_snp_bim \
  $prefix_genotype_snp_bim_file \
  $suffix_genotype_snp_bim_file \
  $path_directory_genetic_reference_prscs \
  $population_ancestry \
  $path_directory_allele_effect \
  $name_file_product_prefix \
  $chromosome_x \
  $threads \
  $path_directory_promiscuity_scripts \
  $path_environment_prscs \
  $path_prscsx \
  $report
fi



#
