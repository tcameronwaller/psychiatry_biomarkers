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
path_file=$1 # full path to file with GWAS summary statistics for a metabolite
path_destination_parent=$2 # full path to destination directory
path_genetic_reference=$3 # full path to genetic reference access
name_prefix=$4 # file name prefix before metabolite identifier or empty string
name_suffix=$5 # file name suffix after metabolite identifier or empty string
path_script_gwas_organization=$6 # full path to script to use for format organization
path_scripts=$7 # full path to scripts for current implementation pipeline
path_promiscuity_scripts=$8 # full path to scripts from promiscuity package

################################################################################
# Derive variables.

# Determine file name.
file_name="$(basename -- $path_file)"
# Determine metabolite identifier.
# Refer to documnetation for test: https://www.freebsd.org/cgi/man.cgi?test
metabolite=${file_name}
if [[ ! -z "$name_prefix" ]]; then
  metabolite=${metabolite/$name_prefix/""}
fi
if [[ ! -z "$name_suffix" ]]; then
  metabolite=${metabolite/$name_suffix/""}
fi
# Report for log.
echo "--------------------------------------------------"
echo "--------------------------------------------------"
echo "--------------------------------------------------"
echo "file: " $file_name
echo "metabolite: " $metabolite
echo "----------"

###########################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_ldsc=$(<"./tools_ldsc.txt")

path_alleles="$path_genetic_reference/alleles"
path_disequilibrium="$path_genetic_reference/disequilibrium"
path_baseline="$path_genetic_reference/baseline"
path_weights="$path_genetic_reference/weights"
path_frequencies="$path_genetic_reference/frequencies"

# Even temporary files need to have names specific to each metabolite.
# During parallel processing, multiple temporary files will exist
# simultaneously.
path_temporary_collection="${path_destination_parent}/temporary_gwas_collection_${metabolite}.txt"
path_temporary_format="${path_destination_parent}/temporary_gwas_format_${metabolite}.txt"
path_temporary_munge="${path_destination_parent}/temporary_munge_${metabolite}"
path_temporary_munge_gwas="${path_temporary_munge}.sumstats.gz"
path_temporary_munge_log="${path_temporary_munge}.log"
path_metabolite_heritability="${path_destination_parent}/heritability_${metabolite}"
path_metabolite_heritability_log="${path_destination_parent}/heritability_${metabolite}"

###########################################################################
# Execute procedure.

# Organize information in format for LDSC.
# Parameters.
report="false" # "true" or "false"
/usr/bin/bash "$path_script_gwas_organization" \
$metabolite \
$file_name \
$path_file \
$path_temporary_collection \
$path_temporary_format \
$path_destination_parent \
$path_promiscuity_scripts \
$report

# Munge metabolite GWAS.
$path_ldsc/munge_sumstats.py \
--sumstats $path_temporary_format \
--out $path_temporary_munge \
--merge-alleles $path_alleles/w_hm3.snplist \
#--a1-inc

# Heritability.
$path_ldsc/ldsc.py \
--h2 $path_temporary_munge_gwas \
--ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--out $path_metabolite_heritability

# Munge phenotype GWAS.

# Genetic correlation.



###########################################################################
# Remove temporary files.
rm $path_temporary_collection
rm $path_temporary_format
rm $path_temporary_munge_gwas
#rm $path_temporary_munge_log
