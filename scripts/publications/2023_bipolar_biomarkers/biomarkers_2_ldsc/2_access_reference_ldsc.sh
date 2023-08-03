#!/bin/bash

################################################################################
# Note.

# Most recent accession of LDSC reference data was on 2 March 2023.

################################################################################
# Organize paths.

# Read private, local file paths.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_reference_ldsc="${path_directory_reference}/ldsc" # accession: TCW; 2 March 2023
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_product_parent="${path_directory_dock}/ldsc_gwas_biomarkers_tcw_2023-05-25"
path_directory_product_child="${path_directory_product_parent}/ldsc"
path_directory_product="${path_directory_product_parent}/2_reference_ldsc"

path_directory_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
path_script_access="${path_directory_promiscuity_scripts}/ldsc/access_ldsc_genetic_references.sh"

# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product_parent
cd $path_directory_product_parent

################################################################################
# Execute procedure.

#/usr/bin/bash "${path_script_access}" \
#$path_directory_product

cp -r $path_directory_reference_ldsc $path_directory_product_parent
mv $path_directory_product_child $path_directory_product

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "2_access_reference_ldsc.sh"
  echo "----------"
fi
