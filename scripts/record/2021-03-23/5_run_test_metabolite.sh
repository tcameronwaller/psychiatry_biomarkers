#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

###########################################################################
# Organize paths.
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")

path_temporary=$(<"./processing_bipolar_metabolism.txt")
path_waller="$path_temporary/waller"
path_bipolar_metabolism="$path_waller/bipolar_metabolism"
path_scripts_record="$path_waller/bipolar_metabolism/scripts/record/2021-03-23"
path_promiscuity_scripts="$path_waller/promiscuity/scripts"

path_dock="$path_waller/dock"
path_genetic_reference="$path_dock/access/genetic_reference"
path_gwas="$path_dock/gwas"
path_heritability="$path_dock/heritability"
path_genetic_correlation="$path_dock/genetic_correlation"

###########################################################################
# Organize variables.

phenotype_study="30124842_yengo_2018" # "30239722_pulit_2018"
metabolite_study="33437055_panyard_2021" # "24816252_shin_2014", "31959995_schlosser_2021", "33437055_panyard_2021"
source_file="metabolite_X57547_meta_analysis_gwas.csv.gz"
name_prefix="metabolite_" # file name prefix before metabolite identifier or empty string
name_suffix="_meta_analysis_gwas.csv.gz" # file name suffix after metabolite identifier or empty string

path_source_directory="${path_gwas_summaries}/${metabolite_study}"
path_source_file="${path_source_directory}/${source_file}"
path_script_gwas_organization="${path_scripts_record}/7_organize_gwas_ldsc_${metabolite_study}.sh"
path_phenotype_gwas="${path_gwas}/${phenotype_study}"
path_study_gwas="${path_gwas}/${metabolite_study}"
path_study_heritability="${path_heritability}/${metabolite_study}"
path_study_genetic_correlation="${path_genetic_correlation}/${phenotype_study}/${metabolite_study}" # notice the directory structure for phenotype and metabolite studies

rm -r $path_study_gwas
rm -r $path_study_heritability
rm -r $path_study_genetic_correlation
mkdir -p $path_study_gwas
mkdir -p $path_study_heritability
mkdir -p $path_study_genetic_correlation

###########################################################################
# Execute procedure.

/usr/bin/bash "$path_scripts_record/6_execute_procedure_metabolite.sh" \
$phenotype_study \
$metabolite_study \
$source_file \
$path_source_file \
$name_prefix \
$name_suffix \
$path_genetic_reference \
$path_phenotype_gwas \
$path_study_gwas \
$path_study_heritability \
$path_study_genetic_correlation \
$path_script_gwas_organization \
$path_promiscuity_scripts \
$report
