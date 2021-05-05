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

metabolite_study=${1} # identifier of GWAS study for metabolites
phenotype_study=${2} # identifier of GWAS study for main phenotype
path_gwas_metabolite=${3} # full path to parent directory for formatted GWAS summary statistics for metabolite
gwas_metabolite_munge_suffix=${4} # file name suffix after metabolite identifier
pattern_gwas_metabolite_munge_suffix=${5} # glob pattern by which to recognize relevant files in source directory
path_gwas_phenotype=${6} # full path to parent directory for formatted GWAS summary statistics for phenotype
path_gwas_phenotype_munge_suffix=${7} # full path to file for formatted and munged GWAS summary statistics for phenotype
path_genetic_correlation_comparison=${8} # full path to parent directory for LDSC genetic correlation estimation
path_genetic_reference=${9} # full path to parent directory with genetic reference files for LDSC
path_promiscuity_scripts=${10} # complete path to directory of general scripts
path_scripts_record=${11} # full path to pipeline scripts
path_ldsc=${12} # full path to parent directory of LDSC executable files
report=${13} # whether to print reports

###########################################################################
# Organize paths.

path_alleles="$path_genetic_reference/alleles"
path_disequilibrium="$path_genetic_reference/disequilibrium"
path_baseline="$path_genetic_reference/baseline"
path_weights="$path_genetic_reference/weights"
path_frequencies="$path_genetic_reference/frequencies"

path_batch_instances="${path_gwas_metabolite}/batch_instances.txt"

###########################################################################
# Execute procedure.

# Organize instances for iteration.
echo "----------------------------------------------------------------------"
echo "Organize array of batch instances."
echo "----------------------------------------------------------------------"
# Collect batch instances.
#cd $path_source
#metabolite_files=(metabolite_*_meta_analysis_gwas.csv.gz)
#for path_file in "${metabolite_files[@]}"; do
#    echo $path_file >> $path_metabolites/metabolite_files.txt
#done
# Define glob pattern for file paths.
# This definition expands to an array of all files in the path directory that
# match the pattern.
path_pattern="${path_gwas_metabolite}/${pattern_gwas_metabolite_munge_suffix}"
# Iterate on all files and directories in parent directory.
rm $path_batch_instances
for path_file in ${path_gwas_metabolite}/*; do
  if [ -f "$path_file" ]; then
    # Current content item is a file.
    # Compare to glob pattern to recognize relevant files.
    #echo $path_file
    #echo $path_file >> $path_batch_instances
    if [[ "$path_file" == ${path_pattern} ]]; then
      # File name matches glob pattern.
      # Include full path to file in batch instances.
      #echo $path_file
      echo $path_file >> $path_batch_instances
    fi
  fi
done

# Read batch instances.
readarray -t batch_instances < $path_batch_instances
batch_instances_count=${#batch_instances[@]}
echo "----------"
echo "count of batch instances: " $batch_instances_count
echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
echo "last batch instance: " ${batch_instances[batch_instances_count - 1]}

# Execute batch iteratively without grid scheduler.
if true; then
  for path_gwas_metabolite_munge_suffix in "${batch_instances[@]}"; do
    # Determine file name.
    #file_name=$path_source_file
    file_name="$(basename -- $path_gwas_metabolite_munge_suffix)"
    # Determine metabolite identifier.
    # Refer to documnetation for test: https://www.freebsd.org/cgi/man.cgi?test
    # Bash script more or less ignores empty string argument, so it does not
    # work well to pass an empty string as an argument.
    # Instead use a non-empty string, such as "null".
    # if [[ ! -z "$name_prefix" ]]; then
    metabolite=${file_name}
    metabolite=${metabolite/$gwas_metabolite_munge_suffix/""}
    # Report.
    if [[ "$report" == "true" ]]; then
      echo "----------------------------------------------------------------------"
      echo "----------------------------------------------------------------------"
      echo "----------------------------------------------------------------------"
      echo "phenotype study: " $phenotype_study
      echo "metabolite study: " $metabolite_study
      echo "path to metabolite file: " $path_gwas_metabolite_munge_suffix
      echo "file: " $file_name
      echo "metabolite: " $metabolite
      echo "----------"
    fi
    # Organize paths.
    path_genetic_correlation_report="${path_genetic_correlation_comparison}/${metabolite}_correlation"
    path_genetic_correlation_report_suffix="${path_genetic_correlation_report}.log"
    # Estimate genetic correlation in LDSC.
    $path_ldsc/ldsc.py \
    --rg $path_gwas_phenotype_munge_suffix,$path_gwas_metabolite_munge_suffix \
    --ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
    --w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
    --out $path_genetic_correlation_report
  done
fi

# Remove temporary files.
rm $path_batch_instances
