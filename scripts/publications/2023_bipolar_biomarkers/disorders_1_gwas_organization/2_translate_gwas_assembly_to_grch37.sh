#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 2 August 2023
# Date, last execution: 31 August 2023
# Date, review: 31 August 2023
################################################################################
# Note

# The purpose of this script is to translate to human genome assembly GRCh37
# (team standard) any GWAS summary statistics that use coordinates from any
# other human genome assembly (NCBI36-hg18, GRCh38, etc).

# On 24 May 2023, TCW chose to prioritize the chain map coordinates from
# Ensembl.

# After running this script, check the product directory to make sure that this
# procedure wrote the appropriate files again, later than all others.



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_source="${path_directory_dock}/gwas_disorders_tcw_2023-08-31/1_gwas_format_standard"
path_directory_product="${path_directory_dock}/gwas_disorders_tcw_2023-08-31/2_gwas_assembly_grch37"
# Files.
#path_file_chain_ncbi36_to_grch37="${path_directory_reference}/crossmap/ucsc/hg18ToHg19.over.chain.gz"
#path_file_chain_grch37_to_grch38="${path_directory_reference}/crossmap/ucsc/hg19ToHg38.over.chain.gz"
#path_file_chain_grch38_to_grch37="${path_directory_reference}/crossmap/ucsc/hg38ToHg19.over.chain.gz"
path_file_chain_ncbi36_to_grch37="${path_directory_reference}/crossmap/ensembl/NCBI36_to_GRCh37.chain.gz"
path_file_chain_grch37_to_grch38="${path_directory_reference}/crossmap/ensembl/GRCh37_to_GRCh38.chain.gz"
path_file_chain_grch38_to_grch37="${path_directory_reference}/crossmap/ensembl/GRCh38_to_GRCh37.chain.gz"
# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_script_map_assembly="${path_directory_partner_scripts}/crossmap/map_gwas_standard_format_bed.sh"

# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
cd $path_directory_product

###########################################################################
# Organize parameters.

threads=1
report="true"

################################################################################
# Execute procedure.

##########
# Copy the GWAS summary statistics from the previous process.
# Most sets of GWAS summary statistics do not need extra processing.
# Subsequent processes on a few studies will replace the appropriate files.
cp $path_directory_source/*.txt.gz $path_directory_product

##########
# Translate genomic assemblies to GRCh37.

##########
# 36477530_saunders_2022 (GRCh38 to GRCh37)

# UCSC:
# Ensembl: 13,268,541 lines to 13,254,829 lines (TCW; 3 August 2023)
/usr/bin/bash $path_script_map_assembly \
"${path_directory_source}/36477530_saunders_2022_alcohol_all.txt.gz" \
"${path_directory_product}/36477530_saunders_2022_alcohol_all.txt.gz" \
$path_file_chain_grch38_to_grch37 \
$threads \
$report

# UCSC:
# Ensembl: ___ lines to ___ lines (TCW; __ August 2023)
/usr/bin/bash $path_script_map_assembly \
"${path_directory_source}/36477530_saunders_2022_alcohol_no_ukb.txt.gz" \
"${path_directory_product}/36477530_saunders_2022_alcohol_no_ukb.txt.gz" \
$path_file_chain_grch38_to_grch37 \
$threads \
$report

# UCSC:
# Ensembl: ___ lines to ___ lines (TCW; __ August 2023)
/usr/bin/bash $path_script_map_assembly \
"${path_directory_source}/36477530_saunders_2022_tobacco_all.txt.gz" \
"${path_directory_product}/36477530_saunders_2022_tobacco_all.txt.gz" \
$path_file_chain_grch38_to_grch37 \
$threads \
$report

# UCSC:
# Ensembl: ___ lines to ___ lines (TCW; __ August 2023)
/usr/bin/bash $path_script_map_assembly \
"${path_directory_source}/36477530_saunders_2022_tobacco_no_ukb.txt.gz" \
"${path_directory_product}/36477530_saunders_2022_tobacco_no_ukb.txt.gz" \
$path_file_chain_grch38_to_grch37 \
$threads \
$report



#
