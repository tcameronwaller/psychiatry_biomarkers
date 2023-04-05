#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date, first execution: 10 March 2023
# Date, last execution: 21 March 2023
################################################################################
################################################################################
################################################################################
# Note

# TODO: TCW; 4 April 2023
# PIPE: iterate on a list of directories all within a main parent directory
# PIPE: for each child directory, combine the SBayesR effects from each chromosome


################################################################################
################################################################################
################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes

path_directory_source_1="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_effects_1"
path_directory_source_2="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_effects_2"
path_directory_source_3="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_effects_3"

path_directory_product_1="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_effects_1_combination"
path_directory_product_2="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_effects_2_combination"
path_directory_product_3="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-21/sbayesr_effects_3_combination"

# Files.
name_file_source_prefix="BMI_GIANTUKB_EUR_"
name_file_source_suffix="_tcw_2023-03-21.snpRes"
path_file_source_2="${path_directory_source_2}/BMI_GIANTUKB_EUR_tcw_2023-03-21.snpRes"
path_file_product_1="${path_directory_product_1}/BMI_GIANTUKB_EUR_tcw_2023-03-21.snpRes"
path_file_product_2="${path_directory_product_2}/BMI_GIANTUKB_EUR_tcw_2023-03-21.snpRes"
path_file_product_3="${path_directory_product_3}/BMI_GIANTUKB_EUR_tcw_2023-03-21.snpRes"

# Scripts.
path_script_combine="${path_directory_process}/promiscuity/scripts/gctb/combine_sbayesr_snp_effect_weights_chromosomes.sh"

# Initialize directories.
rm -r $path_directory_product_1
rm -r $path_directory_product_2
rm -r $path_directory_product_3
mkdir -p $path_directory_product_1
#mkdir -p $path_directory_product_2 # Do not create this directory so that the copy works properly.
mkdir -p $path_directory_product_3
cd $path_directory_product_1



###########################################################################
# Organize parameters.

chromosome_x="false"
report="true"

###########################################################################
# Execute procedure.

##########
# Combine SBayesR SNP effect weights across chromosomes (autosomes).

if true; then
  # 1.
  cd $path_directory_product_1
  /usr/bin/bash $path_script_combine \
  $path_directory_source_1 \
  $name_file_source_prefix \
  $name_file_source_suffix \
  $path_file_product_1 \
  $chromosome_x \
  $report
  # 2.
  cd $path_directory_product_2
  cp -r $path_directory_source_2 $path_directory_product_2
  # 3.
  cd $path_directory_product_3
  /usr/bin/bash $path_script_combine \
  $path_directory_source_3 \
  $name_file_source_prefix \
  $name_file_source_suffix \
  $path_file_product_3 \
  $chromosome_x \
  $report
fi



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "11_combine_effect_weights_gctb_sbayesr_chromosomes.sh"
  echo "----------"
fi



#
