#!/bin/bash

#chmod u+x script.sh

###########################################################################
###########################################################################
###########################################################################
# This script accesses GWAS summary statistics from repositories for
# multiple studies.
# This script will not actually work on the server behind firewall.
# The paths are only for reference.
###########################################################################
###########################################################################
###########################################################################

###########################################################################
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths

path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

# Yengo et al, Human Molecular Genetics, 2018 (PubMed:30124842)
# phenotype: body mass index
if true; then
  # Define destination path.
  path_destination="$path_gwas_summaries/30124842_yengo_2018"
  #rm -r $path_destination
  mkdir -p $path_destination
  # Download files from server.
  cd $path_destination
  wget "https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz"
  wget "https://portals.broadinstitute.org/collaboration/giant/images/0/01/README_summary_statistics_Yengo_et_al_2018.txt"
fi

# Pulit et al, Human Molecular Genetics, 2018 (PubMed:30239722)
# phenotype: waist-to-hip ratio, adjusted for BMI
if true; then
  # Define destination path.
  path_destination="$path_gwas_summaries/30239722_pulit_2018"
  #rm -r $path_destination
  mkdir -p $path_destination
  # Download files from server.
  cd $path_destination
  wget "https://zenodo.org/record/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz"
  wget "https://portals.broadinstitute.org/collaboration/giant/images/a/a4/README_summary_statistics_Pulit_et_al_2018.txt"
fi
