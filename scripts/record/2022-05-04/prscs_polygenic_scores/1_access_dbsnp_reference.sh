#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")

path_dock="${path_process}/dock"
path_dbsnp_reference="${path_dock}/access/dbsnp_reference"

path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_scripts_utility="${path_promiscuity_scripts}/utility"
path_script_access_dbsnp="${path_scripts_utility}/access_dbsnp_reference.sh"

###########################################################################
# Execute procedure.
###########################################################################

/usr/bin/bash "${path_script_access_dbsnp}" \
$path_dbsnp_reference
