#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 23 Decemboer 2022
# Date, last execution: 8 January 2024
# Date, review: 8 January 2023
################################################################################
# Note



################################################################################
# Organize paths.

# Identifiers or designators of parameter version, preparation batch, and
# analysis batch.
identifier_analysis="gwas_2023-12-30_ldsc_2024-01-08"
identifier_parameter="tcw_2023-12-30_dbsnp_rsid"

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"

path_directory_group_parent="${path_directory_dock}/${identifier_analysis}"
path_directory_source="${path_directory_group_parent}/1_gwas_summaries_source"
path_directory_product="${path_directory_group_parent}/3_gwas_format_ldsc"

# Files.

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_file_script="${path_directory_partner_scripts}/ldsc/constrain_translate_gwas_standard_to_ldsc.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

################################################################################
# Organize parameters.

# Report.
report="true"

################################################################################
# Execute procedure.

# Collect files.
#cd $path_directory_source
# Bash version 4.4 introduced the "-d" option for "readarray".
#readarray -d "" -t paths_files_source < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
paths_files_source=()
while IFS= read -r -d $'\0'; do
  paths_files_source+=("$REPLY")
done < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
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
  name_file_source="$(basename $path_file_source)"
  path_file_product="${path_directory_product}/${name_file_source}" # hopefully unique
  # Translate GWAS summary statistics to format for LDSC.
  /usr/bin/bash "${path_file_script}" \
  $path_file_source \
  $path_file_product \
  $report
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "Script path:"
    echo $path_file_script
    echo "Source file path:"
    echo $path_file_source
    echo "Product file path:"
    echo $path_file_product
    echo "----------"
  fi
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "3_translate_gwas_to_ldsc_format.sh"
  echo "----------"
fi



#
