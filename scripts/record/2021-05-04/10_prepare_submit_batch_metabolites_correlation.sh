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
path_metabolite_gwas_source_directory=${4} # full path to parent directory for raw GWAS summary statistics for metabolites in study
metabolite_file_pattern=${5} # glob pattern by which to recognize relevant files in source directory
metabolite_file_prefix=${6} # file name prefix before metabolite identifier or "null"
metabolite_file_suffix=${7} # file name suffix after metabolite identifier or "null"
path_genetic_reference=${8} # full path to parent directory with genetic reference files for LDSC
path_phenotype_gwas=${9} # full path to parent directory for formatted GWAS summary statistics for phenotype
path_phenotype_gwas_munge_suffix=${10} # full path to file for formatted and munged GWAS summary statistics for phenotype
path_study_gwas=${11} # full path to parent directory for formatted GWAS summary statistics for metabolites in study
path_study_heritability=${12} # full path to parent directory for LDSC heritability estimation for metabolites in study
path_study_genetic_correlation=${13} # full path to parent directory for LDSC genetic correlation for metabolites in study
path_scripts_record=${14} # full path to pipeline scripts
path_promiscuity_scripts=${15} # complete path to directory of general scripts
path_script_gwas_format=${16} # full path to script to use to organize format of GWAS summary statistics for metabolites in study
path_ldsc=${16}
report=${17} # whether to print reports

################################################################################
# Organize variables.

path_batch_instances="${path_study_gwas}/batch_instances.txt"
# Define glob pattern for file paths.
# This definition expands to an array of all files in the path directory that
# match the pattern.
path_pattern="${path_metabolite_gwas_source_directory}/${metabolite_file_pattern}"

# Initialize directories.
#rm -r $path_study_gwas
#rm -r $path_study_heritability
rm -r $path_study_genetic_correlation
mkdir -p $path_study_gwas
mkdir -p $path_study_heritability
mkdir -p $path_study_genetic_correlation

# Report.
echo "----------------------------------------------------------------------"
echo "Report."
echo "----------------------------------------------------------------------"
echo "----------"
echo "phenotype study: " $phenotype_study
echo "metabolite study: " $metabolite_study
echo "source path directory: " $path_metabolite_gwas_source_directory
echo "path to batch instances: " $path_batch_instances
echo "----------"

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
# Iterate on all files and directories in parent directory.
rm $path_batch_instances
for path_file in $path_metabolite_gwas_source_directory/*; do
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

if false; then
  # Submit array batch to Sun Grid Engine.
  # Array batch indices must start at one (not zero).
  echo "----------------------------------------------------------------------"
  echo "Submit array of batches to Sun Grid Engine."
  echo "----------------------------------------------------------------------"
  qsub -t 1-${batch_instances_count}:1 -o \
  "${path_study_gwas}/out.txt" -e "${path_study_gwas}/error.txt" \
  "${path_scripts_record}/6_run_batch_organize_gwas_ldsc_heritability.sh" \
  $path_batch_instances \
  $batch_instances_count \
  $phenotype_study \
  $metabolite_study \
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
fi
if false; then
  for path_source_file in "${batch_instances[@]}"; do
    /usr/bin/bash "${path_scripts_record}/7_execute_procedure_metabolite.sh" \
    $phenotype_study \
    $metabolite_study \
    $path_source_file \
    $name_prefix \
    $name_suffix \
    $path_genetic_reference \
    $path_phenotype_gwas \
    $path_study_gwas \
    $path_study_heritability \
    $path_study_genetic_correlation \
    $path_script_gwas_format \
    $path_promiscuity_scripts \
    $report
  done
fi
