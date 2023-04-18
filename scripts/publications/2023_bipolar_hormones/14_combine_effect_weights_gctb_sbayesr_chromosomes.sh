#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 18 April 2022
# Date, last execution: 18 April 2023
# Review: TCW; 18 April 2023
################################################################################
# Note

# TODO: TCW; 4 April 2023 <-- This is just 1 option...
# PIPE: iterate on a list of directories all within a main parent directory
# PIPE: for each child directory, combine the SBayesR effects from each chromosome


################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_source="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/effect_weights_sbayesr_2023-04-18"
path_directory_product="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/effect_weights_sbayesr_2023-04-18_combination"

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
# Organize multi-dimensional array of information about studies.
studies=()
# [name of study for file prefix without delimiter] ; [name file suffix with delimiter]

studies+=(
  "32769997_zhou_2020;\
  _2023-04-18.snpRes"
)
studies+=(
  "30367059_teumer_2018_tsh_female;\
  _2023-04-18.snpRes"
)
#studies+=(
#  "30367059_teumer_2018_tsh_male;\
#  _2023-04-18.snpRes"
#)
studies+=(
  "30367059_teumer_2018_ft4_all;\
  _2023-04-18.snpRes"
)
studies+=(
  "30367059_teumer_2018_ft4_female;\
  _2023-04-18.snpRes"
)
studies+=(
  "30367059_teumer_2018_ft4_male;\
  _2023-04-18.snpRes"
)


# Organize information for each study.
for study in "${studies[@]}"; do
  # Separate fields from instance.
  # [name of study for file prefix without delimiter] ; [name file suffix with delimiter]
  IFS=";" read -r -a array <<< "${study}"
  prefix="${array[0]}"
  suffix="${array[1]}"
  # Organize paths and parameters.
  name_file_source_prefix="${prefix}_"
  name_file_source_suffix=$suffix
  path_file_product="${path_directory_product}/${prefix}.snpRes"
  # Combine SBayesR SNP effect weights across chromosomes (autosomes).
  /usr/bin/bash $path_script_combine \
  $path_directory_source \
  $name_file_source_prefix \
  $name_file_source_suffix \
  $path_file_product \
  $chromosome_x \
  $report
done



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script completion:"
  echo $0 # Print full file path to script.
  echo "14_combine_effect_weights_gctb_sbayesr_chromosomes.sh"
  echo "----------"
fi



#
