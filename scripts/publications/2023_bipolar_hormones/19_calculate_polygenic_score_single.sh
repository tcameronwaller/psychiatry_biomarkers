#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 15 March 2023
# Date, last execution: 5 April 2023
# Review: TCW; 5 April 2023
################################################################################
# Note

# When calculating polygenic scores on large cohorts of genotypes, this process
# requires considerable computational resources and multiple hours of time.

# Write a batch submission script to handle the parallelization across
# chromosomes.

# TODO: TCW; 4 April 2023
# PIPE: submit one SLURM batch for each set of SBayesR effects
# PIPE: each SLURM batch should include autosomes 1-22


################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
#path_directory_source_genotypes="${path_directory_dock}/genotypes_mayo_bipolar_disorder_1_2_merge"
path_directory_source_genotypes="${path_directory_dock}/genotypes_regeneron_mayo_bipolar_disorder_1_2_merge_mayo_control"
path_directory_source_effects="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/sbayesr_effects_grch38"
path_directory_product="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/sbayesr_plink_polygenic_scores_short_queue_16g"

# Files.
path_file_source_effects="${path_directory_source_effects}/BMI_GIANTUKB_EUR_1_grch38_standard.txt.gz"
name_file_product_prefix="32769997_zhou_2020_tsh_chromosome_"
name_file_product_suffix="_tcw" # only include suffix for base file name without the final suffix of format or type
#path_file_product="${path_directory_product}/${name_file_product_prefix}${chromosome}${name_file_product_suffix}"
#name_file_genotypes_prefix="MERGED.maf0.dosR20.3.noDups.chr"
#name_file_genotypes_suffix=".dose.vcf.gz"
name_file_genotypes_prefix="GWAS_MERGED_BPphe_wMCBBctrl.maf0.dosR20.3.noDups.noSM.chr"
name_file_genotypes_suffix=".dose.vcf.gz"

# Scripts.
path_script_submit_batch="${path_directory_process}/promiscuity/scripts/plink/1_submit_batch_polygenic_score_chromosomes.sh"

# Initialize directories.
rm -r $path_directory_product # Only remove previous if starting new.
mkdir -p $path_directory_product
cd $path_directory_product

###########################################################################
# Organize parameters.

report="true"

###########################################################################
# Execute procedure.

/usr/bin/bash $path_script_submit_batch \
$path_file_source_effects \
$path_directory_source_genotypes \
$name_file_genotypes_prefix \
$name_file_genotypes_suffix \
$path_directory_product \
$name_file_product_prefix \
$name_file_product_suffix \
$report



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script completion:"
  echo $0 # Print full file path to script.
  echo "19_calculate_polygenic_score_single.sh"
  echo "----------"
fi



#
