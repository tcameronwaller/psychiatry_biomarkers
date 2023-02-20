#!/bin/bash

###########################################################################
# Organize paths.

# Read private, local file paths.
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")

path_dock="${path_process}/dock"
path_directory_reference="${path_dock}/hormone_genetics_tcw_2023-02-17/reference_ldsc"

path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_access="${path_promiscuity_scripts}/ldsc/access_ldsc_genetic_references.sh"

###########################################################################
# Execute procedure.
###########################################################################

/usr/bin/bash "${path_script_access}" \
$path_directory_reference
