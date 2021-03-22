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
metabolite_study=${2} # identifier of GWAS study for metabolites
source_file=${3} # name of source file with GWAS summary statistics for a single metabolite
path_source_file=${4} # full path to source file with GWAS summary statistics for a single metabolite
name_prefix=${5} # file name prefix before metabolite identifier or empty string
name_suffix=${6} # file name suffix after metabolite identifier or empty string
path_genetic_reference=${7} # full path to parent directory with genetic reference files for LDSC
path_phenotype_gwas=${8} # full path to parent directory for formatted GWAS summary statistics for phenotype
path_study_gwas=${9} # full path to parent directory for formatted GWAS summary statistics for metabolites in study
path_study_heritability=${10} # full path to parent directory for LDSC heritability estimation for metabolites in study
path_study_genetic_correlation=${11} # full path to parent directory for LDSC genetic correlation for metabolites in study
path_script_gwas_organization=${12} # full path to script to use to organize format of GWAS summary statistics for metabolites in study
path_promiscuity_scripts=${13} # complete path to directory of scripts for z-score standardization
report=${14} # whether to print reports

################################################################################
# Derive variables.

# Determine file name.
#file_name=$source_file
file_name="$(basename -- $path_source_file)"
# Determine metabolite identifier.
# Refer to documnetation for test: https://www.freebsd.org/cgi/man.cgi?test
metabolite=${file_name}
if [[ ! -z "$name_prefix" ]]; then
  metabolite=${metabolite/$name_prefix/""}
fi
if [[ ! -z "$name_suffix" ]]; then
  metabolite=${metabolite/$name_suffix/""}
fi
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo "phenotype study: " $phenotype_study
  echo "metabolite study: " $metabolite_study
  echo "path to metabolite file: " $path_source_file
  echo "file: " $file_name
  echo "metabolite: " $metabolite
  echo "----------"
fi

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

path_phenotype_gwas_munge_suffix="${path_phenotype_gwas}/gwas_munge.sumstats.gz"

# Even temporary files need to have names specific to each metabolite.
# During parallel processing, multiple temporary files will exist
# simultaneously.
path_gwas_collection="${path_study_gwas}/gwas_collection_${metabolite}.txt"
path_gwas_format="${path_study_gwas}/gwas_format_${metabolite}.txt"
path_gwas_format_compress="${path_gwas_format}.gz"
path_gwas_munge="${path_study_gwas}/gwas_munge_${metabolite}"
path_gwas_munge_suffix="${path_gwas_munge}.sumstats.gz"
path_gwas_munge_log="${path_gwas_munge}.log"

path_heritability_report="${path_study_heritability}/heritability_${metabolite}"
path_heritability_report_suffix="${path_heritability_report}.log"

path_genetic_correlation_report="${path_study_genetic_correlation}/correlation_${metabolite}"
path_genetic_correlation_report_suffix="${path_genetic_correlation_report}.log"

###########################################################################
# Execute procedure.

# Organize information in format for LDSC.
# Parameters.
report="false" # "true" or "false"
/usr/bin/bash "$path_script_gwas_organization" \
$metabolite \
$path_source_file \
$path_gwas_collection \
$path_gwas_format \
$path_gwas_format_compress \
$path_promiscuity_scripts \
$report

# Munge metabolite GWAS.
$path_ldsc/munge_sumstats.py \
--sumstats $path_gwas_format \
--out $path_gwas_munge \
--merge-alleles $path_alleles/w_hm3.snplist \
#--a1-inc

# Heritability.
$path_ldsc/ldsc.py \
--h2 $path_gwas_munge_suffix \
--ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--out $path_heritability_report

# Genetic correlation between metabolite and phenotype.
$path_ldsc/ldsc.py \
--rg $path_phenotype_gwas_munge_suffix,$path_gwas_munge_suffix \
--ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--out $path_genetic_correlation_report

###########################################################################
# Remove temporary files.
rm $path_gwas_collection
rm $path_gwas_format
rm $path_gwas_munge_suffix
rm $path_gwas_munge_log
