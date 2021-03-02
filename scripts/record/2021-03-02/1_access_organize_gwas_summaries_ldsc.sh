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
rm -r $path_heritability
if [ ! -d $path_heritability ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_heritability
    mkdir -p $path_heritability_shin_2014
    mkdir -p $path_heritability_schlosser_2021
    mkdir -p $path_heritability_panyard_2021
fi

# Format of GWAS summary statistics for LDSC.
# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics
# description: ............................ LDSC column
# variant identifier: ........................ "SNP"
# alternate allele (effect allele): .......... "A1"
# reference allele (non-effect allele): ...... "A2"
# sample size: ............................... "N"
# effect (coefficient or odds ratio): ........ "BETA" or "OR"
# probability (p-value): ..................... "P"

# Format of GWAS summary statistics for PRS-CS.
# https://github.com/getian107/PRScs
# description: ............................ LDSC column
# variant identifier: ........................ "SNP"
# alternate allele (effect allele): .......... "A1"
# reference allele (non-effect allele): ...... "A2"
# effect (coefficient or odds ratio): ........ "BETA" or "OR"
# probability (p-value): ..................... "P"

################################################################################

# PubMed: 33437055; Author: Panyard; Year: 2021.
echo "----------------------------------------------------------------------"
echo "PubMed: 33437055; Author: Panyard; Year: 2021"
echo "Human genome version: GRCh37, hg19"
echo "----------------------------------------------------------------------"
cd $path_33437055_panyard_2021
metabolite_files=(metabolite_*_meta_analysis_gwas.csv.gz)
count=${#metabolite_files[@]}
echo "count of file paths: " $count
#rm $path_metabolites/metabolite_files.txt
#for path_file in "${metabolite_files[@]}"; do
#    echo $path_file >> $path_metabolites/metabolite_files.txt
#done
# Define paths to temporary files for each iteration.
path_parent="$path_heritability_panyard_2021"
path_temporary_gwas_format="$path_heritability_panyard_2021/temporary_gwas_format.txt"
path_temporary_gwas_format_zip="$path_heritability_panyard_2021/temporary_gwas_format.txt.gz"
# Define glob pattern to recognize relevant files.
pattern="${path_33437055_panyard_2021}/metabolite_*_meta_analysis_gwas.csv.gz"
# Iterate on all files and directories in parent directory.
for path_file in $path_33437055_panyard_2021/*; do
  if [ -f "$path_file" ]; then
    # Current content item is a file.
    #echo $file
    if [[ "$path_file" == ${pattern} ]]; then
      # File name matches glob pattern.
      file_name="$(basename -- $path_file)"
      # Organize information in format for LDSC.
      # Parameters.
      report="true"
      /usr/bin/bash "$path_scripts/2_organize_gwas_ldsc_33437055_panyard_2021.sh" \
      $file_name \
      $path_file \
      $path_temporary_gwas_format \
      $path_temporary_gwas_format_zip \
      $path_parent \
      $path_calculate_z_score_column_5_of_6 \
      $report
      # Munge.
      # Heritability.
    fi
  fi
done
