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

phenotype_study=${1} # identifier of GWAS study for main phenotype
path_gwas_phenotype=${2} # full path to parent directory for formatted GWAS summary statistics for phenotype
path_gwas_phenotype_munge_suffix=${3} # full path to file for formatted and munged GWAS summary statistics for phenotype
path_gwas_cohorts_models=${4} # full path to parent directory for formatted GWAS summary statistics for cohorts and hormones
file_gwas_cohorts_models_munge_suffix=${5} # name of file for formatted and munged GWAS summary statistics for cohorts and hormones
path_genetic_correlation=${6} # full path to parent directory for LDSC genetic correlation estimation
path_genetic_reference=${7} # full path to parent directory with genetic reference files for LDSC
path_promiscuity_scripts=${8} # complete path to directory of general scripts
path_scripts_record=${9} # full path to pipeline scripts
path_ldsc=${10} # full path to parent directory of LDSC executable files
report=${11} # whether to print reports

###########################################################################
# Organize paths.

path_alleles="$path_genetic_reference/alleles"
path_disequilibrium="$path_genetic_reference/disequilibrium"
path_baseline="$path_genetic_reference/baseline"
path_weights="$path_genetic_reference/weights"
path_frequencies="$path_genetic_reference/frequencies"

###########################################################################
# Execute procedure.

# Iterate on directories for GWAS on cohorts and hormones.
cd $path_gwas_cohorts_models
for path_directory in `find . -maxdepth 1 -mindepth 1 -type d -not -name .`; do
  if [ -d "$path_directory" ]; then
    # Current content item is a directory.
    directory="$(basename -- $path_directory)"

    echo $directory

    # Determine whether directory contains file for GWAS summary statistics
    # after concatenation, format, and munge.
    path_gwas_cohorts_models_munge_suffix="${path_gwas_cohorts_models}/${directory}/${file_gwas_cohorts_models_munge_suffix}"
    if [[ -f "$path_gwas_cohorts_models_munge_suffix" ]]; then

      echo $path_gwas_cohorts_models_munge_suffix

      # Genetic correlation.
      if true; then

        # Organize paths.
        path_genetic_correlation_comparison="${path_genetic_correlation}/${phenotype_study}/cohorts_hormones/${directory}"
        path_genetic_correlation_report="${path_genetic_correlation_comparison}/correlation"
        path_genetic_correlation_report_suffix="${path_genetic_correlation_report}.log"
        # Initialize directories.
        rm -r $path_genetic_correlation_comparison
        mkdir -p $path_genetic_correlation_comparison
        # Estimate genetic correlation in LDSC.
        $path_ldsc/ldsc.py \
        --rg $path_gwas_phenotype_munge_suffix,$path_gwas_cohorts_models_munge_suffix \
        --ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
        --w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
        --out $path_genetic_correlation_report
        # Report.
        if [[ "$report" == "true" ]]; then
          echo "----------------------------------------------------------------------"
          echo "----------------------------------------------------------------------"
          echo "----------------------------------------------------------------------"
          echo "phenotype study: " $phenotype_study
          echo "cohort-hormone study: " $directory
          echo "path to cohort-hormone file: " $path_gwas_cohorts_models_munge_suffix
          echo "----------"
        fi
      fi
    fi
  fi
done
