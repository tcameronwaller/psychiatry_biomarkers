#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 25 May 2023
# Date, last execution: 8 January 2024
# Date, review: 8 January 2023
################################################################################
# Note

# Most recent accession of LDSC reference data was on 2 March 2023.

################################################################################
# Organize paths.

# Identifiers or designators of parameter version, preparation batch, and
# analysis batch.
identifier_analysis="gwas_2023-12-30_ldsc_2024-01-08"
identifier_parameter="tcw_2023-12-30_dbsnp_rsid"

# Read private, local file paths.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_reference_ldsc="${path_directory_reference}/ldsc" # accession: TCW; 2 March 2023
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"

path_directory_group_parent="${path_directory_dock}/${identifier_analysis}"
path_directory_product_child="${path_directory_group_parent}/ldsc"
path_directory_product="${path_directory_group_parent}/2_reference_ldsc"

#path_directory_partner_scripts="${path_directory_process}/partner/scripts"
#path_script_access="${path_directory_partner_scripts}/ldsc/access_ldsc_genetic_references.sh"

# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_group_parent
cd $path_directory_group_parent

################################################################################
# Execute procedure.

#/usr/bin/bash "${path_script_access}" \
#$path_directory_product

cp -r $path_directory_reference_ldsc $path_directory_group_parent
mv $path_directory_product_child $path_directory_product

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "2_access_reference_ldsc.sh"
  echo "----------"
fi
