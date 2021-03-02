#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "Access summary statistics from GWAS on human metabolome."
echo "Organize files in format suitable for LDSC."
echo "version 1"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo ""
echo ""
echo ""

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
path_scripts="$path_waller/bipolar_metabolism/scripts/record/2021-03-02"
path_promiscuity_scripts="$path_waller/promiscuity/scripts"
path_calculate_z_score_column_5_of_6="$path_promiscuity_scripts/calculate_z_score_column_5_of_6.sh"
path_calculate_z_score_column_4_of_5="$path_promiscuity_scripts/calculate_z_score_column_4_of_5.sh"

path_dock="$path_waller/dock"
path_heritability="$path_dock/heritability"
path_heritability_shin_2014="$path_heritability/24816252_shin_2014"
path_heritability_schlosser_2021="$path_heritability/31959995_schlosser_2021"
path_heritability_panyard_2021="$path_heritability/33437055_panyard_2021"

# Initialize directories.
#rm -r $path_heritability
if [ ! -d $path_heritability ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_heritability
fi

################################################################################
# Computational time
# Running on head node
# Iteration for each metabolite
# format and standardization of GWAS summary statistics: ~ 45 seconds
# LDSC munge of GWAS summary statistics: ???
# LDSC heritability: ???

# TODO: I DEFINITELY need to parallelize this on the cluster... organize batch jobs...
# TODO: print the format and z-score reports to log file for each metabolite...

################################################################################
# PubMed: 33437055; Author: Panyard; Year: 2021.
echo "--------------------------------------------------"
echo "--------------------------------------------------"
echo "--------------------------------------------------"
echo "PubMed: 33437055; Author: Panyard; Year: 2021"
echo "Human genome version: GRCh37, hg19"
echo "destination path: " $path_heritability_panyard_2021
echo "--------------------------------------------------"
echo "--------------------------------------------------"
echo "--------------------------------------------------"
# Parameters.
path_source=$path_33437055_panyard_2021
file_pattern="metabolite_*_meta_analysis_gwas.csv.gz"
path_destination_parent=${path_heritability_panyard_2021}
path_script_gwas_organization="${path_scripts}/5_organize_gwas_ldsc_33437055_panyard_2021.sh"
# Prepare and submit batch.
/usr/bin/bash "$path_scripts/3_prepare_submit_batch_organize_gwas_heritability.sh" \
$path_source \
$file_pattern \
$path_destination_parent \
$path_script_gwas_organization \
$path_scripts \
$path_promiscuity_scripts
