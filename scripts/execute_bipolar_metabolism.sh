#!/bin/bash

#chmod u+x script.sh

# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_python_library=$(<"./project_tools_python_library.txt")
path_temporary=$(<"./processing_bipolar_metabolism.txt")
path_waller="$path_temporary/waller"
path_dock="$path_waller/dock"
path_bipolar_metabolism="$path_waller/bipolar_metabolism"
path_package="$path_bipolar_metabolism/bipolar_metabolism"

# Organize paths to custom package installations.
PYTHONPATH=$path_python_library:$PYTHONPATH
export PYTHONPATH

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

# Collect and organize heritability estimations for metabolites from GWAS
# summary statistics of multiple studies on the human metabolome.
python3 $path_package/interface.py main --path_dock $path_dock --heritability

# Collect and aggregate genetic scores for metabolites across the UK Biobank.
#python3 $path_package/interface.py main --path_dock $path_dock --aggregation

# Assemble phenotypic information from UK Biobank.
#python3 $path_package/interface.py main --path_dock $path_dock --assembly

# Organize phenotypic information from UK Biobank.
#python3 $path_package/interface.py main --path_dock $path_dock --organization

# Regress phenotypes against metabolites' genetic scores across the UK Biobank.
#python3 $path_package/interface.py main --path_dock $path_dock --association


# TODO: new module...
# TODO: start with function for 1 metabolite...
# TODO: combine metabolite genetic score(s) with UK Biobank phenotypes
# TODO: adapt regression code from "bimodality"...

# TODO: organize UK Biobank phenotypes for Bipolar Disorder, Depression, etc
