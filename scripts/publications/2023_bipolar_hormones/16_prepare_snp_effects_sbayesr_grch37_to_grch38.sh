#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 18 April 2022
# Date, last execution: 18 April 2023
# Review: TCW; 18 April 2023
################################################################################
# Note



################################################################################
# Organize arguments.

path_directory_source=${1} # full path to parent directory for source files
name_base_file_source=${2} # base name of file for source effect weights from SBayesR
suffix_grch37_sbayesr=${3} # suffix in name of file in assembly GRCh37 in format of SBayesR
suffix_grch37_bed=${4} # suffix in name of file in assembly GRCh37 in format of UCSC BED
suffix_grch38_bed=${5} # suffix in name of file in assembly GRCh38 in format of UCSC BED
suffix_grch38_standard=${6} # suffix in name of file in assembly GRCh38 in format of team standard
path_directory_product=${7} # full path to parent directory for product files in human genome assembly GRCh38
report=${8} # whether to print reports



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes

# Files.
path_file_grch37_sbayesr="${path_directory_source}/${name_base_file_source}${suffix_grch37_sbayesr}"
path_file_grch37_bed="${path_directory_product}/${name_base_file_source}${suffix_grch37_bed}"
path_file_grch38_bed="${path_directory_product}/${name_base_file_source}${suffix_grch38_bed}"
path_file_grch38_standard="${path_directory_product}/${name_base_file_source}${suffix_grch38_standard}"
#path_file_chain_grch37_to_grch38="${path_directory_reference}/assembly_chains/ucsc/hg19ToHg38.over.chain.gz"
path_file_chain_grch37_to_grch38="${path_directory_reference}/crossmap/ensembl/GRCh37_to_GRCh38.chain.gz"

# Scripts.
path_script_ucsc_bed="${path_directory_process}/promiscuity/scripts/gctb/translate_snp_effects_sbayesr_to_ucsc_bed.sh"
path_script_map="${path_directory_process}/promiscuity/scripts/crossmap/map_genomic_feature_bed.sh"
#path_script_standard="${path_directory_process}/promiscuity/scripts/gctb/translate_snp_effects_ucsc_bed_to_standard.sh"
path_script_standard="${path_directory_process}/promiscuity/scripts/gctb/translate_snp_effects_ucsc_bed_to_standard_identifier.sh"

# Initialize directories.
cd $path_directory_product



###########################################################################
# Organize parameters.

threads=1
report="true"

###########################################################################
# Execute procedure.

##########
# Translate SBayesR SNP effect weights format from SBayesR to CrossMap UCSC BED.

if true; then
  # 1.
  /usr/bin/bash $path_script_ucsc_bed \
  $path_file_grch37_sbayesr \
  $path_file_grch37_bed \
  $report
fi

##########
# Translate SBayesR SNP effect weights in CrossMap from GRCh37 to GRCh38.

if true; then
  # 1.
  /usr/bin/bash $path_script_map \
  $path_file_grch37_bed \
  $path_file_grch38_bed \
  $path_file_chain_grch37_to_grch38 \
  $threads \
  $report
fi

##########
# Translate SBayesR SNP effect weights format from CrossMap UCSC BED to team standard with special identifiers.
# The format of variant (SNP) identifiers must match the target genotypes.

if true; then
  # 1.
  /usr/bin/bash $path_script_standard \
  $path_file_grch38_bed \
  $path_file_grch38_standard \
  $report
fi



#
