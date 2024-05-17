#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

################################################################################
# Author: T. Cameron Waller, PhD
# Date, first execution: 1 December 2020
# Date, last execution: 14 March 2024
# Review: TCW; 14 March 2024
################################################################################
# Note

# This Bash script calls Python scripts or interfaces to more extensive
# packages and modules. Paths in this script are specific to execution on a
# specific server.


################################################################################
# Organize paths.

echo "read private file path variables and organize paths..."
cd ~/paths
path_tools=$(<"./waller_tools.txt")
path_environment_main="${path_tools}/python/environments/main"
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"
path_package="${path_process}/psychiatry_biomarkers/psychiatry_biomarkers"


################################################################################
# Execute procedure(s).

# Echo each command to console.
set -x

# Determine whether the temporary directory structure already exists.
if [ ! -d $path_dock ]
then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_dock
fi

##########
# Activate Python Virtual Environment.
source "${path_environment_main}/bin/activate"
echo "----------"
echo "...Confirm Python Virtual Environment path..."
which python3
sleep 1s
echo "----------"

##########
# Call Python procedures.

# Routine: main
python3 $path_package/interface.py main --path_dock $path_dock --extraction_ldsc
#python3 $path_package/interface.py main --path_dock $path_dock --assembly_rg
#python3 $path_package/interface.py main --path_dock $path_dock --assembly_prs
#python3 $path_package/interface.py main --path_dock $path_dock --plot

# Routine: sexy_age


# Notes below this point are about execution of old and obsolete procedures (TCW; 12 March 2024).

# Routine: uk_biobank
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --assembly # TCW; 07 November 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --importation # TCW; 31 March 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --organization # TCW; 07 June 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --stratification # TCW; 09 April 2022 <-- no longer executable
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --genotype # TCW; 10 November 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --description # TCW; 12 December 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --regression # TCW; 1 November 2022
#python3 $path_package/interface.py uk_biobank --path_dock $path_dock --collection # TCW; 20 April 2022

# Routine: stragglers

#python3 $path_package/interface.py stragglers --path_dock $path_dock --scrap # TCW; 05 October 2022
#python3 $path_package/interface.py stragglers --path_dock $path_dock --mbpdb_prioritize_supplement # TCW; 12 September 2022
#python3 $path_package/interface.py stragglers --path_dock $path_dock --mbpdb_assembly # TCW; 16 May 2023
#python3 $path_package/interface.py stragglers --path_dock $path_dock --mbpdb_organization # TCW; 24 April 2023
#python3 $path_package/interface.py stragglers --path_dock $path_dock --mbpdb_description # TCW; 16 December 2022
#python3 $path_package/interface.py stragglers --path_dock $path_dock --mbpdb_regression # TCW; 24 April 2023
#python3 $path_package/interface.py stragglers --path_dock $path_dock --mcita_assembly # TCW; 06 July 2022


##########
# Deactivate Python Virtual Environment.
deactivate
echo "----------"
echo "confirm deactivation of Python Virtual Environment..."
which python3



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "execute_psychiatry_biomarkers.sh"
  echo "----------"
fi



#
