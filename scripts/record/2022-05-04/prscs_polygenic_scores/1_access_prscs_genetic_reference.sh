#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")

path_dock="${path_process}/dock"
path_genetic_reference="${path_dock}/access/genetic_reference_prscs"

path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_promiscuity_scripts_polygenic_score="${path_promiscuity_scripts}/polygenic_score"

###########################################################################
# Execute procedure.
###########################################################################

/usr/bin/bash "${path_promiscuity_scripts_polygenic_score}/access_prscs_genetic_references.sh" \
$path_genetic_reference
