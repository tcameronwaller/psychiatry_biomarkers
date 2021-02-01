#!/bin/bash

#chmod u+x script.sh

###########################################################################
###########################################################################
###########################################################################
# This script organizes paths and parameters to access local information
# about persons and phenotypes in the UKBiobank.
###########################################################################
###########################################################################
###########################################################################


###########################################################################
# Organize script parameters.
#project="sexy_alcohol"
project="bipolar_metabolism"

###########################################################################
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_temporary=$(<"./processing_bipolar_metabolism.txt")
path_dock="$path_temporary/waller/dock"
path_variables="$path_dock/parameters/${project}/uk_biobank_access_variables.txt"
path_access="$path_dock/access/ukbiobank_phenotypes"
path_scripts="$path_temporary/waller/uk_biobank/scripts"

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

# Determine whether the temporary directory structure already exists.
if [ ! -d $path_access ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_dock
    mkdir -p $path_access
fi

# Specify UKBiobank phenotype variables of interest.
# Only use this method if not using the file in the parameter directory.
#echo "specify phenotype variables to access from UKBiobank"
#echo "31 22001 21022 21002 30600" | tr -s " " "\n" > $path_variables

# Access phenotype variables and auxiliary information from UKBiobank.
/usr/bin/bash "$path_scripts/access_ukbiobank_phenotypes.sh" \
$path_variables \
$path_access
