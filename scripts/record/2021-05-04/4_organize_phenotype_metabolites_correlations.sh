#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize argument variables.

phenotype_study=${1} # identifier of GWAS study for phenotype
path_phenotype_gwas=${2} # full path to parent directory for formatted GWAS summary statistics for phenotype
path_phenotype_gwas_munge_suffix=${3} # full path to file for formatted and munged GWAS summary statistics for phenotype
path_gwas=${4} # full path to parent directory for formatted GWAS summary statistics for phenotype
path_heritability=${5} # full path to parent directory for LDSC heritability estimation
path_genetic_correlation=${6} # full path to parent directory for LDSC genetic correlation estimation
path_genetic_reference=${7} # full path to parent directory with genetic reference files for LDSC
path_promiscuity_scripts=${8} # complete path to directory of general scripts
path_scripts_record=${9} # full path to pipeline scripts
path_ldsc=${10}
report=${11} # whether to print reports

################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_scripts_format="${path_promiscuity_scripts}/format_gwas_ldsc"

###########################################################################
# Execute procedure.

metabolite_study="27005778_kettunen_2016"
path_metabolite_gwas_source_directory="${path_gwas_summaries}/${metabolite_study}"
metabolite_file_pattern="Summary_statistics_MAGNETIC_*.txt.gz" # do not expand with full path yet
metabolite_file_prefix="Summary_statistics_MAGNETIC_" # file name prefix before metabolite identifier or "null"
metabolite_file_suffix=".txt.gz" # file name suffix after metabolite identifier or "null"

path_study_gwas="${path_gwas}/${metabolite_study}"
path_study_heritability="${path_heritability}/${metabolite_study}"
path_study_genetic_correlation="${path_genetic_correlation}/${phenotype_study}/${metabolite_study}" # notice the directory structure for phenotype and metabolite studies

path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${metabolite_study}.sh"

# Prepare and submit batch.
/usr/bin/bash "${path_scripts_record}/5_prepare_submit_batch_metabolites_correlation.sh" \
$phenotype_study \
$metabolite_study \
$path_metabolite_gwas_source_directory \
$metabolite_file_pattern \
$metabolite_file_prefix \
$metabolite_file_suffix \
$path_genetic_reference \
$path_phenotype_gwas \
$path_phenotype_gwas_munge_suffix \
$path_study_gwas \
$path_study_heritability \
$path_study_genetic_correlation \
$path_scripts_record \
$path_promiscuity_scripts \
$path_script_gwas_format \
$path_ldsc \
$report
