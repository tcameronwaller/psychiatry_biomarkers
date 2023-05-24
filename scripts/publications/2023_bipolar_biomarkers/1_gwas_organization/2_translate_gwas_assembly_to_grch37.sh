#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 23 December 2022
# Date, last execution: 24 May 2023
################################################################################
# Note

# The purpose of this script is to translate to human genome assembly GRCh37
# (team standard) any GWAS summary statistics that are use coordinates from any
# other human genome assembly (NCBI36-hg18, GRCh38, etc).



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_source="${path_directory_dock}/test_crossmap_gwas_summary/gwas_source"
path_directory_product="${path_directory_dock}/test_crossmap_gwas_summary/gwas_product"
# Files.
#path_file_chain_ncbi36_to_grch37="${path_directory_reference}/crossmap/ucsc/hg18ToHg19.over.chain.gz"
#path_file_chain_grch37_to_grch38="${path_directory_reference}/crossmap/ucsc/hg19ToHg38.over.chain.gz"
#path_file_chain_grch38_to_grch37="${path_directory_reference}/crossmap/ucsc/hg38ToHg19.over.chain.gz"
path_file_chain_ncbi36_to_grch37="${path_directory_reference}/crossmap/ensembl/NCBI36_to_GRCh37.chain.gz"
path_file_chain_grch37_to_grch38="${path_directory_reference}/crossmap/ensembl/GRCh37_to_GRCh38.chain.gz"
path_file_chain_grch38_to_grch37="${path_directory_reference}/crossmap/ensembl/GRCh38_to_GRCh37.chain.gz"
# Scripts.
path_directory_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
path_script_map_assembly="${path_directory_promiscuity_scripts}/crossmap/map_gwas_standard_format_bed.sh"


###########################################################################
# Organize parameters.

threads=1
report="true"

################################################################################
# Execute procedure.


##########
# 36635386_chen_2023 (GRCh38 to GRCh37)

# UCSC: 15,428,167 lines to 15,396,790 lines (TCW; 24 May 2023)
# Ensembl: 15,428,167 lines to 15,384,987 lines
/usr/bin/bash $path_script_map_assembly \
"${path_directory_source}/36635386_chen_2023_cortisol.txt.gz" \
"${path_directory_product}/36635386_chen_2023_cortisol.txt.gz" \
$path_file_chain_grch38_to_grch37 \
$threads \
$report

#/usr/bin/bash $path_script_map_assembly \
#"${path_directory_source}/36635386_chen_2023_thyroxine_total.txt.gz" \
#"${path_directory_product}/36635386_chen_2023_thyroxine_total.txt.gz" \
#$path_file_chain_grch38_to_grch37 \
#$threads \
#$report

##########
# 34662886_backman_2021 (GRCh38 to GRCh37)

#/usr/bin/bash $path_script_map_assembly \
#"${path_directory_source}/34662886_backman_2021_albumin.txt.gz" \
#"${path_directory_product}/34662886_backman_2021_albumin.txt.gz" \
#$path_file_chain_grch38_to_grch37 \
#$threads \
#$report

##########
# 34017140_mbatchou_2021 (GRCh38 to GRCh37)

#/usr/bin/bash $path_script_map_assembly \
#"${path_directory_source}/34017140_mbatchou_2021_albumin.txt.gz" \
#"${path_directory_product}/34017140_mbatchou_2021_albumin.txt.gz" \
#$path_file_chain_grch38_to_grch37 \
#$threads \
#$report

##########
# 24586183_medici_2014 (NCBI36-hg18 to GRCh37)

# UCSC: 2,425,175 lines to 2,424,714 lines (TCW; 24 May 2023)
# Ensembl: 2,425,175 lines to 2,424,987 lines
/usr/bin/bash $path_script_map_assembly \
"${path_directory_source}/24586183_medici_2014_thyroid_peroxidase_antibody.txt.gz" \
"${path_directory_product}/24586183_medici_2014_thyroid_peroxidase_antibody.txt.gz" \
$path_file_chain_ncbi36_to_grch37 \
$threads \
$report

#/usr/bin/bash $path_script_map_assembly \
#"${path_directory_source}/24586183_medici_2014_thyroid_peroxidase_reactivity.txt.gz" \
#"${path_directory_product}/24586183_medici_2014_thyroid_peroxidase_reactivity.txt.gz" \
#$path_file_chain_ncbi36_to_grch37 \
#$threads \
#$report



#
