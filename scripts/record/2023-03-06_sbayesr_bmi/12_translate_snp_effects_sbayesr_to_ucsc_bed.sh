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


# TODO: TCW; 10 March 2023
# 1. translate to BED format for CrossMap
# 2. map in CrossMap GRCh37 to GRCh38
# 3. translate to format for PLINK2 --score
# 4. apply PLINK2 --score



################################################################################
################################################################################
################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes

path_directory_source="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/combination_effect_weights_sbayesr"
path_directory_product="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/snp_effects_ucsc_bed"

# Files.
path_file_source="${path_directory_source}/BMI_GIANTUKB_EUR_tcw_2023-03-01_chromosomes.snpRes"
path_file_product="${path_directory_product}/BMI_GIANTUKB_EUR_tcw_2023-03-01_chromosomes.bed.gz"

# Scripts.
path_script_translate="${path_directory_process}/promiscuity/scripts/gctb/translate_snp_effects_sbayesr_to_ucsc_bed.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product



###########################################################################
# Organize parameters.

report="true"

###########################################################################
# Execute procedure.

##########
# Combine SBayesR SNP effect weights across chromosomes (autosomes).

if true; then
  /usr/bin/bash $path_script_translate \
  $path_file_source \
  $path_file_product \
  $report
fi



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "12_translate_snp_effects_sbayesr_to_ucsc_bed.sh"
  echo "----------"
fi



#
