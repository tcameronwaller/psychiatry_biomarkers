#!/bin/bash


# TODO: TCW; 2 March 2023
# TODO: Move the LDSC reference information to the "storage/reference" directory.

###########################################################################
# Organize paths.

# Read private, local file paths.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")

path_directory_dock="${path_directory_process}/dock"
path_directory_reference="${path_directory_dock}/ldsc_gwas_biomarkers_tcw_2023-05-25/2_reference_ldsc"

path_directory_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
path_script_access="${path_directory_promiscuity_scripts}/ldsc/access_ldsc_genetic_references.sh"

###########################################################################
# Execute procedure.
###########################################################################

/usr/bin/bash "${path_script_access}" \
$path_directory_reference

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "2_access_reference_ldsc.sh"
  echo "----------"
fi
