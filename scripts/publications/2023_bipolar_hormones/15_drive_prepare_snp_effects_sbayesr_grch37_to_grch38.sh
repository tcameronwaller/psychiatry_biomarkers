#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 18 April 2022
# Date, last execution: 18 April 2023
# Review: TCW; 18 April 2023
################################################################################
# Note



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_source="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/effect_weights_sbayesr_2023-04-18_combination"
path_directory_product="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/sbayesr_effects_grch38"

# Scripts.
path_script_prepare="${path_directory_process}/psychiatric_metabolism/scripts/publications/2023_bipolar_hormones/16_prepare_snp_effects_sbayesr_grch37_to_grch38.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product



###########################################################################
# Organize parameters.

suffix_grch37_sbayesr=".snpRes"
suffix_grch37_bed=".bed.gz"
suffix_grch38_bed=".bed.gz"
suffix_grch38_standard=".txt.gz"
threads=1
report="true"

###########################################################################
# Execute procedure.

# Collect files.
#cd $path_directory_source
# Bash version 4.4 introduced the "-d" option for "readarray".
#readarray -d "" -t paths_files_source < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
paths_files_source=()
while IFS= read -r -d $'\0'; do
  paths_files_source+=("$REPLY")
done < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*${suffix_grch37_sbayesr}" -print0)
count_paths_files_source=${#paths_files_source[@]}

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Source directory:"
  echo $path_directory_source
  echo "count of files: " $count_paths_files_source
  echo "first file: " ${paths_files_source[0]} # notice base-zero indexing
  echo "last file: " ${paths_files_source[$count_paths_files_source - 1]}
  echo "Product directory:"
  echo $path_directory_product
  echo "----------"
fi

for path_file_source in "${paths_files_source[@]}"; do
  # Extract base name of file.
  name_base_file_source="$(basename $path_file_source $suffix_grch37_sbayesr)"
  # Translate SBayesR SNP effect weights format from SBayesR to CrossMap UCSC BED.
  # Translate SBayesR SNP effect weights in CrossMap from GRCh37 to GRCh38.
  # Translate SBayesR SNP effect weights format from CrossMap UCSC BED to team standard with special identifiers.
  # The format of variant (SNP) identifiers must match the target genotypes.
  /usr/bin/bash $path_script_prepare \
  $path_directory_source \
  $name_base_file_source \
  $suffix_grch37_sbayesr \
  $suffix_grch37_bed \
  $suffix_grch38_bed \
  $suffix_grch38_standard \
  $path_directory_product \
  $report
done



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script completion:"
  echo $0 # Print full file path to script.
  echo "15_drive_prepare_snp_effects_sbayesr_grch37_to_grch38.sh"
  echo "----------"
fi



#
