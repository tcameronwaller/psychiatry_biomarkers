#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

# TODO: TCW; 25 July 2022
# TODO: "psychiatric_metabolism" and "sexy_alcohol" will become the main driver packages with the "interface.py".
# TODO: These main driver packages will have subpackages: "uk_biobank", "local_studies"
# TODO: Change name of the "main" parser to "uk_biobank", then call uk_biobank procedures directly from "interface.py"
# TODO:


################################################################################
# Activate Python Virtual Environment.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_tools=$(<"./waller_tools.txt")
path_environment_main="${path_tools}/python/environments/main"
source "${path_environment_main}/bin/activate"
echo "----------"
echo "confirm activation of Python Virtual Environment..."
which python3
sleep 5s

################################################################################
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"
path_package="${path_process}/psychiatric_metabolism/psychiatric_metabolism"

# Echo each command to console.
set -x

# Determine whether the temporary directory structure already exists.
if [ ! -d $path_dock ]
then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_dock
fi

# Execute procedure(s).

# Routine: main
#python3 $path_package/interface.py main --path_dock $path_dock --scratch

# Routine: uk_biobank
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --assembly # TCW; 07 November 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --importation # TCW; 31 March 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --organization # TCW; 07 June 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --stratification # TCW; 09 April 2022 <-- no longer executable
python3 $path_package/interface.py uk_biobank --path_dock $path_dock --genotype # TCW; 10 November 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --description # TCW; 01 June 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --description # TCW; 1 November 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --regression # TCW; 1 November 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --collection # TCW; 20 April 2022

# Routine: stragglers

#python3 $path_package/interface.py stragglers --path_dock $path_dock --scrap # TCW; 05 October 2022
#python3 $path_package/interface.py stragglers --path_dock $path_dock --mbpdb_prioritize_supplement # TCW; 12 September 2022
#python3 $path_package/interface.py stragglers --path_dock $path_dock --mbpdb_assembly # TCW; 14 September 2022
#python3 $path_package/interface.py stragglers --path_dock $path_dock --mbpdb_organization # TCW; 14 September 2022
#python3 $path_package/interface.py stragglers --path_dock $path_dock --mbpdb_regression # TCW; 22 September 2022
#python3 $path_package/interface.py stragglers --path_dock $path_dock --mcita_assembly # TCW; 06 July 2022

################################################################################
# Deactivate Python Virtual Environment.
deactivate
echo "----------"
echo "confirm deactivation of Python Virtual Environment..."
which python3
