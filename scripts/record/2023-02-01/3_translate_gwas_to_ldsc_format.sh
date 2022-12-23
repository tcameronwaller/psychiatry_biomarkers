#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date: 23 December 2022
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
path_dock="${path_directory_process}/dock"
path_directory_source="${path_dock}/hormone_genetics/gwas_format_standard"
path_directory_product="${path_dock}/hormone_genetics/gwas_format_ldsc"
# Files.

# Scripts.
path_directory_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
path_file_script="${path_directory_promiscuity_scripts}/gwas_process/format_gwas_ldsc/constrain_translate_gwas_team.sh"

# Initialize directories.
#rm -r $path_directory_product
#mkdir -p $path_directory_product

################################################################################
# Organize parameters.

# Report.
report="true"

################################################################################
# Execute procedure.

#cd $path_directory_source
readarray -d "" -t paths_files_source < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.txt.gz" -print0)
for path_file_source in "${paths_files_source[@]}"; do
  echo $path_file_source
done

################################################################################

if false; then
  # Iterate on files.
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
fi

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "3_translate_gwas_to_ldsc_format.sh"
  echo "----------"
fi



#
