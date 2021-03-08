#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

###########################################################################
# Organize paths.
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_24816252_shin_2014="$path_gwas_summaries/24816252_shin_2014"
path_31959995_schlosser_2021="$path_gwas_summaries/31959995_schlosser_2021"
path_33437055_panyard_2021="$path_gwas_summaries/33437055_panyard_2021"

path_temporary=$(<"./processing_bipolar_metabolism.txt")
path_waller="$path_temporary/waller"
path_bipolar_metabolism="$path_waller/bipolar_metabolism"
path_scripts="$path_waller/bipolar_metabolism/scripts/record/2021-03-09"
path_promiscuity_scripts="$path_waller/promiscuity/scripts"

path_dock="$path_waller/dock"
path_genetic_reference="$path_dock/access/genetic_reference"
path_heritability="$path_dock/heritability"
path_heritability_shin_2014="$path_heritability/24816252_shin_2014"
path_heritability_schlosser_2021="$path_heritability/31959995_schlosser_2021"
path_heritability_panyard_2021="$path_heritability/33437055_panyard_2021"

###########################################################################
# Organize variables.

path_file="${path_33437055_panyard_2021}/metabolite_X57547_meta_analysis_gwas.csv.gz"
path_destination_parent=${path_heritability_panyard_2021}
name_prefix="metabolite_" # file name prefix before metabolite identifier or empty string
name_suffix="_meta_analysis_gwas.csv.gz" # file name suffix after metabolite identifier or empty string
path_script_gwas_organization="${path_scripts}/6_organize_gwas_ldsc_33437055_panyard_2021.sh"

# Initialize directories.
#rm -r $path_heritability
if [ ! -d $path_heritability ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_heritability
fi
#rm -r $path_destination_parent
if [ ! -d $path_destination_parent ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_destination_parent
fi

###########################################################################
# Execute procedure.

$path_scripts/5_execute_procedure_metabolite.sh \
$path_file \
$path_destination_parent \
$path_genetic_reference \
$name_prefix \
$name_suffix \
$path_script_gwas_organization \
$path_scripts \
$path_promiscuity_scripts
