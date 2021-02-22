#!/bin/bash

#chmod u+x script.sh

###########################################################################
###########################################################################
###########################################################################
# This script accesses GWAS summary statistics from repositories for
# multiple studies.
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
if true: then
  # Define file destination.
  path_access_study="$path_dock/access/gwas_summaries/33437055_panyard_2021"
  mkdir -p $path_access_study
  # Download files from server.
  # curl -u user:password "path-to-original-file-online" -o "path-to-new-local-copy"
  cd $path_access
  wget -r -l 0 "ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/MWAS/Panyard_GWAS_sumstats_CSF_metabolites/"

  #cp -r "${path_access}/ftp.biostat.wisc.edu/pub/lu_group/Projects/MWAS/Panyard_GWAS_sumstats_CSF_metabolites" \
  #$path_access

  # Metabolite information is in Supplemental Table 2.
fi
