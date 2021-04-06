#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize paths.
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")

path_temporary=$(<"./processing_bipolar_metabolism.txt")
path_waller="$path_temporary/waller"
path_bipolar_metabolism="$path_waller/bipolar_metabolism"
path_scripts_record="$path_waller/bipolar_metabolism/scripts/record/2021-04-06"
path_promiscuity_scripts="$path_waller/promiscuity/scripts"
path_scripts_format="$path_waller/promiscuity/scripts/format_gwas_ldsc"

path_dock="$path_waller/dock"
path_genetic_reference="$path_dock/access/genetic_reference"
path_gwas="$path_dock/gwas"
path_heritability="$path_dock/heritability"
path_genetic_correlation="$path_dock/genetic_correlation"

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

#path_source_directory="${path_gwas_summaries}/${metabolite_study}" # path unique to 24816252_shin_2014
#name_prefix="metabolite_" # file name prefix before metabolite identifier or "null"
#name_suffix="_meta_analysis_gwas.csv.gz" # file name suffix after metabolite identifier or "null"
#file_pattern="metabolite_*_meta_analysis_gwas.csv.gz" # do not expand with full path yet



################################################################################
# PubMed: 24816252; Author: Shin; Year: 2014.
echo "--------------------------------------------------"
echo "--------------------------------------------------"
echo "--------------------------------------------------"
echo "PubMed: 24816252; Author: Shin; Year: 2014"
echo "--------------------------------------------------"
echo "--------------------------------------------------"
echo "--------------------------------------------------"
# Parameters.
# "30124842_yengo_2018", "30239722_pulit_2018",
# "31043756_stahl_2019", "30718901_howard_2019", "30482948_walters_2018",
# "00000000_pgc3_2021_all", "00000000_pgc3_2021_bd1", "00000000_pgc3_2021_bd2",
# "29906448_ruderfer_2018_scz_vs_bpd", "29906448_ruderfer_2018_scz_bpd_vs_ctl",
phenotype_study="29906448_ruderfer_2018_scz_bpd_vs_ctl"
metabolite_study="24816252_shin_2014"
path_source_directory="${path_gwas_summaries}/${metabolite_study}/metabolites_meta" # path unique to 24816252_shin_2014
name_prefix="null" # file name prefix before metabolite identifier or "null"
name_suffix=".metal.pos.txt.gz" # file name suffix after metabolite identifier or "null"
source_file_pattern="*.metal.pos.txt.gz" # do not expand with full path yet
path_phenotype_gwas="${path_gwas}/${phenotype_study}"
path_study_gwas="${path_gwas}/${metabolite_study}"
path_study_heritability="${path_heritability}/${metabolite_study}"
path_study_genetic_correlation="${path_genetic_correlation}/${phenotype_study}/${metabolite_study}" # notice the directory structure for phenotype and metabolite studies
path_script_gwas_format="${path_scripts_format}/format_gwas_ldsc_${metabolite_study}.sh"
report="false" # "true" or "false"
# Prepare and submit batch.
/usr/bin/bash "${path_scripts_record}/5_prepare_submit_batch_organize_gwas_heritability.sh" \
$phenotype_study \
$metabolite_study \
$source_file_pattern \
$path_source_directory \
$name_prefix \
$name_suffix \
$path_genetic_reference \
$path_phenotype_gwas \
$path_study_gwas \
$path_study_heritability \
$path_study_genetic_correlation \
$path_scripts_record \
$path_script_gwas_format \
$path_promiscuity_scripts \
$report
