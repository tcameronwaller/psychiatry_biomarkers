#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")

path_dock="${path_process}/dock"
path_genetic_reference="${path_dock}/access/linkage_disequilibrium_ldsc"

path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_promiscuity_scripts_ldsc="${path_promiscuity_scripts}/gwas_process/ldsc_genetic_heritability_correlation"

###########################################################################
# Execute procedure.
###########################################################################

/usr/bin/bash "${path_promiscuity_scripts_ldsc}/access_ldsc_genetic_references.sh" \
$path_genetic_reference
