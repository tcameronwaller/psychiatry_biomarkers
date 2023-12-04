#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 27 Decemboer 2022
# Date, last execution: 4 December 2023
# Date, review: 4 December 2023
################################################################################
# Note

# The "prior_1" version of GWAS summary statistics from study
# "37872160_williams_2023" failed the LDSC munge procedure on 30 November 2023
# because these summary statistics originally lacked rsIDs for SNPs. Hopefully
# it will be possible to run a successful procedure in GWAS2VCF.

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"

path_directory_group_parent="${path_directory_dock}/ldsc_gwas_tcw_2023-11-26_2023-12-04"
path_directory_reference="${path_directory_group_parent}/2_reference_ldsc"
path_directory_source="${path_directory_group_parent}/3_gwas_format_ldsc"
path_directory_product="${path_directory_group_parent}/4_gwas_munge_ldsc"
path_file_alleles="${path_directory_reference}/alleles/w_hm3.snplist"

# Files.

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_directory_ldsc="${path_directory_partner_scripts}/ldsc"
path_file_script="${path_directory_ldsc}/munge_gwas_ldsc.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.

response="coefficient"
threads=8
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
  name_base_file_product="$(basename $path_file_source .txt.gz)"
  path_file_base_product="${path_directory_product}/${name_base_file_product}" # hopefully unique
  # Translate GWAS summary statistics to format for LDSC.
  /usr/bin/bash $path_file_script \
  $path_file_source \
  $path_file_base_product \
  $path_file_alleles \
  $response \
  $threads \
  $report
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "Script path:"
    echo $path_file_script
    echo "Source file path:"
    echo $path_file_source
    echo "Product file path:"
    echo $path_file_base_product
    echo "----------"
  fi
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "4_munge_gwas_ldsc.sh"
  echo "----------"
fi



#
