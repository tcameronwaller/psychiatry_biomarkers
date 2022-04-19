#!/bin/bash

################################################################################
# Activate Virtual Environment.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_tools=$(<"./waller_tools.txt")
path_ldsc=$(<"./tools_ldsc.txt")
path_environment_ldsc="${path_tools}/python/environments/ldsc"
source "${path_environment_ldsc}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python2
sleep 5s

path_process=$(<"./process_psychiatric_metabolism.txt")
path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2022-05-04/ldsc_heritability_correlation"

python2 $path_scripts_record/4_test_environment_memory.py

################################################################################
# Deactivate Virtual Environment.
deactivate
which python2


#
