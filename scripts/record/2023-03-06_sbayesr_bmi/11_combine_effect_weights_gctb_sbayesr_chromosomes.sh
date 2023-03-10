#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date, first execution: 10 March 2023
# Date, last execution: 10 March 2023
################################################################################
################################################################################
################################################################################
# Note



################################################################################
################################################################################
################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes

path_directory_source="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/effect_weights_sbayesr"
path_directory_product="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/combination_effect_weights_sbayesr"

# Files.
name_file_source_prefix="BMI_GIANTUKB_EUR_"
name_file_source_suffix="_tcw_2023-03-01.snpRes"
path_file_product="${path_directory_product}/BMI_GIANTUKB_EUR_tcw_2023-03-01_chromosomes.snpRes"

# Scripts.
path_script_combine="${path_directory_process}/promiscuity/scripts/gctb/combine_sbayesr_snp_effect_weights_chromosomes.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product



###########################################################################
# Organize parameters.

chromosome_x="false"
report="true"

###########################################################################
# Execute procedure.

##########
# Combine SBayesR SNP effect weights across chromosomes (autosomes).

if true; then
  /usr/bin/bash $path_script_combine \
  $path_directory_source \
  $name_file_source_prefix \
  $name_file_source_suffix \
  $path_file_product \
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
