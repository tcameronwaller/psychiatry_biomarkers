#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 15 March 2023
# Date, last execution: 4 April 2023
# Review: TCW; 21 March 2023
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
path_directory_source_effects="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_effects_grch38"
path_directory_product_1="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_1_polygenic_scores"
path_directory_product_2="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_2_polygenic_scores"
path_directory_product_3="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_3_polygenic_scores"

# Files.
path_file_source_effects_1="${path_directory_source_effects}/BMI_GIANTUKB_EUR_1_grch38_standard.txt.gz"
path_file_source_effects_2="${path_directory_source_effects}/BMI_GIANTUKB_EUR_2_grch38_standard.txt.gz"
path_file_source_effects_3="${path_directory_source_effects}/BMI_GIANTUKB_EUR_3_grch38_standard.txt.gz"

#name_file_genotypes_prefix="MERGED.maf0.dosR20.3.noDups.chr"
#name_file_genotypes_suffix=".dose.vcf.gz"
name_file_genotypes_prefix="GWAS_MERGED_BPphe_wMCBBctrl.maf0.dosR20.3.noDups.noSM.chr"
name_file_genotypes_suffix=".dose.vcf.gz"
name_file_product_prefix="BMI_GIANTUKB_EUR_chromosome_"
name_file_product_suffix=""
#path_file_product="${path_directory_product}/${name_file_product_prefix}${chromosome}${name_file_product_suffix}"

# Scripts.
path_script_drive_calculate="${path_directory_process}/promiscuity/scripts/plink/calculate_polygenic_score.sh"

# Initialize directories.
#rm -r $path_directory_product_1 # Only remove previous if starting new.
#rm -r $path_directory_product_2 # Only remove previous if starting new.
#rm -r $path_directory_product_3 # Only remove previous if starting new.
mkdir -p $path_directory_product_1
mkdir -p $path_directory_product_2
mkdir -p $path_directory_product_3
cd $path_directory_product_1

###########################################################################
# Organize parameters.

chromosome_x="false"
threads=8
report="true"

###########################################################################
# Execute procedure.

# Iterate on relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
else
  #chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22") # temporarily disable
  chromosomes=("11") # temporarily finish chromosomes
fi
for chromosome in "${chromosomes[@]}"; do
  # Define paths and names to files for current chromosome.
  # Files.
  path_file_source_genotypes="${path_directory_source_genotypes}/${name_file_genotypes_prefix}${chromosome}${name_file_genotypes_suffix}"
  name_base_file_product="${name_file_product_prefix}${chromosome}${name_file_product_suffix}"
  # Calculate polygenic scores for a single chromosome.
  # 1.
  #cd $path_directory_product_1
  #/usr/bin/bash $path_script_drive_calculate \
  #$path_file_source_effects_1 \
  #$path_file_source_genotypes \
  #$path_directory_product_1 \
  #$name_base_file_product \
  #$threads \
  #$report
  # 2.
  #cd $path_directory_product_2
  #/usr/bin/bash $path_script_drive_calculate \
  #$path_file_source_effects_2 \
  #$path_file_source_genotypes \
  #$path_directory_product_2 \
  #$name_base_file_product \
  #$threads \
  #$report
  # 3.
  cd $path_directory_product_3
  /usr/bin/bash $path_script_drive_calculate \
  $path_file_source_effects_3 \
  $path_file_source_genotypes \
  $path_directory_product_3 \
  $name_base_file_product \
  $threads \
  $report
done



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "15_calculate_polygenic_scores.sh"
  echo "----------"
fi



#
