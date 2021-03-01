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
path_temporary=$(<"./processing_bipolar_metabolism.txt")
path_dock="$path_temporary/waller/dock"
path_access="$path_dock/access/gwas_summaries"

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

# Determine whether the temporary directory structure already exists.
if [ ! -d $path_access ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_access
fi

# Panyard et al, Communications Biology, 2021 (PubMed:33437055)
# ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/MWAS/
# Metabolite information is in Supplemental Table 2.
if true; then
  # Define file destination.
  path_access_study="$path_dock/access/gwas_summaries/33437055_panyard_2021"
  mkdir -p $path_access_study
  # Download files from server.
  # curl -u user:password "path-to-original-file-online" -o "path-to-new-local-copy"
  cd $path_access
  wget -r -l 0 "ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/MWAS/Panyard_GWAS_sumstats_CSF_metabolites/"

  #cp -r "${path_access}/ftp.biostat.wisc.edu/pub/lu_group/Projects/MWAS/Panyard_GWAS_sumstats_CSF_metabolites" \
  #$path_access

fi

# Schlosser et al, Nature Genetics, 2020 (PubMed:31959995)
# https://www.ebi.ac.uk/gwas/publications/31959995
# Metabolite information is in Supplemental Table 2.
if true; then
  # Define file destination.
  path_access_study="$path_dock/access/gwas_summaries/31959995_schlosser_2021"
  mkdir -p $path_access_study

  # Download files from server.
  # curl -u user:password "path-to-original-file-online" -o "path-to-new-local-copy"
  cd $path_access
  wget -r -l 0 "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/SchlosserP_31959995_GCST009733/"

fi

# Shin et al, Nature Genetics, 2014 (PubMed:24816252)
# https://www.ebi.ac.uk/gwas/publications/24816252
# Metabolite information is in metabolite map file???
if true; then
  # Define file destination.
  path_access_study="$path_dock/access/gwas_summaries/24816252_shin_2014"
  mkdir -p $path_access_study

  # Download files from server.
  # curl -u user:password "path-to-original-file-online" -o "path-to-new-local-copy"
  cd $path_access
  #curl "sftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/ShinSY_24816252_GCST002443/"
  wget --no-parent --no-passive -r -l 10 "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/ShinSY_24816252_GCST002443/"

fi
