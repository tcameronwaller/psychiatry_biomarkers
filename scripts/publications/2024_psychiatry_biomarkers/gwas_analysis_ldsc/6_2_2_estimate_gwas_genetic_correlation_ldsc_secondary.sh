#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 16 May 2024
# Date, last execution: 16 May 2024
# Date, review: ___ 2024
################################################################################
# Note

# SLURM batch job: ___ (group: "secondaries"; instances: 26,896; date: 15 May 2024)
# TODO: set up the SLURM driver scripts to accommodate 26,896 instances in 3 sub-batches (<= 10,000 per batch)


##########
# Note: TCW; 16 May 2024
# The file "table_gwas_translation_tcw_2023-12-30_dbsnp_rsid.tsv" provided
# parameters for the preparation and organization of a collection of 244 files
# of GWAS summary statistics. Of these 244 files of GWAS summary statistics, 80
# 'primaries' were for psychiatric disorders or substance use disorders or their
# related traits, and 164 'secondaries' were for thyroid disorders, thyroid
# biomarkers, or other biomarkers including sex hormones and their binding
# proteins. On 8 January 2024, TCW ran LDSC munge on all 244 files of GWAS
# summary statistics. On 8 and 9 January 2024, TCW ran LDSC SNP heritability
# estimates on all 244 munged files of GWAS summary statistics. On 9 January
# 2024, TCW ran LDSC genetic correlation estimates on pairwise combinations
# between all primary-primary studies and all primary-secondary studies. On 31
# January 2024, TCW ran LDSC genetic correlation estimates on pairwise
# combinations between all secondary-secondary studies relevant to thyroid
# disorders and thyroid biomarkers. On 16 May 2024, TCW ran LDSC genetic
# correlation estimates on pairwise combinations between all secondary-secondary
# studies.

##########
# On 15 May 2024, TCW set up directories and files within processing space.
# 1. Create parent directory within processing space for new genetic correlations.
#cd ${path_dock_directory_in_process_space}
#mkdir ./gwas_2023-12-30_ldsc_2024-01-08_extra_2024-05-15 # New files created here will be copied to storage space.
# 2. Copy LDSC reference files from storage space to processing space.
#cp -r ${path_directory_reference_ldsc} $"{path_directory_dock}/gwas_2023-12-30_ldsc_2024-01-08_extra_2024-05-15"
#mv $"{path_directory_dock}/gwas_2023-12-30_ldsc_2024-01-08_extra_2024-05-15/ldsc" $"{path_directory_dock}/gwas_2023-12-30_ldsc_2024-01-08_extra_2024-05-15/2_reference_ldsc"
# size: 3.5G
# 3. Copy LDSC munge files from storage space to processing space.
#cp -r $"{path_parent_directory_in_storage_space}/4_gwas_munge_ldsc" $"{path_directory_dock}/gwas_2023-12-30_ldsc_2024-01-08_extra_2024-05-15"
# size: 2.1G
# count content files: 487

################################################################################
# Organize paths.

# Identifiers or designators of parameter version, preparation batch, and
# analysis batch.
identifier_analysis="gwas_2023-12-30_ldsc_2024-01-08_extra_2024-05-15"
identifier_parameter="tcw_2023-12-30_dbsnp_rsid"

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"

path_directory_group_parent="${path_directory_dock}/${identifier_analysis}"
path_directory_reference="${path_directory_group_parent}/2_reference_ldsc"
path_directory_source_primary="${path_directory_group_parent}/4_gwas_munge_ldsc"
path_directory_source_secondary="${path_directory_group_parent}/4_gwas_munge_ldsc"
path_directory_product_parent="${path_directory_group_parent}/6_gwas_correlation_ldsc_secondary"
path_directory_product_child="${path_directory_product_parent}/physiology_biomarkers"
path_directory_disequilibrium="${path_directory_reference}/disequilibrium/eur_w_ld_chr"
path_directory_batch="${path_directory_product_parent}/batch"

# Files.
path_file_batch_instances="${path_directory_batch}/batch_instances.txt"
#path_file_batch_out="${path_directory_batch}/batch_out.txt"
#path_file_batch_error="${path_directory_batch}/batch_error.txt"

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_directory_ldsc="${path_directory_partner_scripts}/ldsc"
path_file_script_ldsc_correlation="${path_directory_ldsc}/estimate_gwas_genetic_correlation_ldsc.sh"
path_file_script_ldsc_correlation_batch_1="${path_directory_ldsc}/ldsc_correlation_batch_1.sh"

# Initialize directories.
rm -r $path_directory_product_parent # caution
rm -r $path_directory_batch # caution
mkdir -p $path_directory_product_parent
mkdir -p $path_directory_product_child
mkdir -p $path_directory_batch

# Initialize files.
rm $path_file_batch_instances

################################################################################
# Organize parameters.



##########
# Common parameters.
threads=2
report="true"



##########
# Secondary studies.

# Define array of secondary studies.
secondaries=()

# All secondary studies, including thyroid disorders, thyroid biomarkers,
# sex hormones and their binding proteins, and other biomarkers.
# Review: TCW; ___ 2024
# Count of secondary studies and versions: __

# TODO: fill in the identifiers of secondary studies here...
# TODO: there ought to be identifiers for 164 studies
# TODO: 164 X 164 = 26,896


##########
# Organize multi-dimensional array of information about comparisons.
# [full path to base name of product file] ; \
# [full path to primary source file of LDSC munge GWAS summary statistics] ; \
# [full path to secondary source file of LDSC munge GWAS summary statistics]

comparisons=()
#comparisons+=(
#  "${path_directory_product}/hypothyroidism_against_hyperthyroidism;\
#  ${path_directory_source_primary}/30367059_teumer_2018_hypothyroidism.sumstats.gz;\
#  ${path_directory_source_secondary}/30367059_teumer_2018_hyperthyroidism.sumstats.gz"
#)

if true; then
  # Assemble array of batch instance details.
  for primary in "${secondaries[@]}"; do
    for secondary in "${secondaries[@]}"; do
      # Organize paths.
      name_comparison="${primary}_-_${secondary}"
      path_file_base_product="${path_directory_product_child}/${name_comparison}"
      path_file_source_primary="${path_directory_source_primary}/${primary}.sumstats.gz"
      path_file_source_secondary="${path_directory_source_secondary}/${secondary}.sumstats.gz"
      # Assemble parameters for comparison.
      comparisons+=("${path_file_base_product};${path_file_source_primary};${path_file_source_secondary}")
    done
  done
fi



################################################################################
# Report.

count_comparisons=${#comparisons[@]}

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Source directory:"
  echo $path_directory_source
  echo "count of comparisons: " $count_comparisons
  echo "first file: " ${comparisons[0]} # notice base-zero indexing
  echo "last file: " ${comparisons[$count_comparisons - 1]}
  echo "----------"
fi

sleep 5s

################################################################################
# Execute procedure.



##########
# Simple iteration.
if false; then
  for comparison in "${comparisons[@]}"; do
    # Separate fields from instance.
    # [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]
    IFS=";" read -r -a array <<< "${comparison}"
    path_file_base_product="${array[0]}"
    path_file_source_primary="${array[1]}"
    path_file_source_secondary="${array[2]}"
    # Estimate Genetic Correlation by LDSC.
    /usr/bin/bash $path_file_script_ldsc_correlation \
    $path_file_source_primary \
    $path_file_source_secondary \
    $path_file_base_product \
    $path_directory_disequilibrium \
    $threads \
    $report
    # Report.
    if [[ "$report" == "true" ]]; then
      echo "----------"
      echo "Script path:"
      echo $path_script
      echo "Product file path:"
      echo $path_file_base_product
      echo "Primary source file path:"
      echo $path_file_source_primary
      echo "Secondary source file path:"
      echo $path_file_source_secondary
      echo "----------"
    fi
  done
fi



##########
# Batch parallelization.
if true; then
  # Organize batch job instances.
  for comparison in "${comparisons[@]}"; do
    # Define parameters in array instance for batch job.
    echo $comparison >> $path_file_batch_instances
  done
  # Call first script in series for batch execution.
  /usr/bin/bash $path_file_script_ldsc_correlation_batch_1 \
  $path_file_batch_instances \
  $path_directory_batch \
  $path_directory_product_parent \
  $path_directory_disequilibrium \
  $path_directory_process \
  $threads \
  $report
fi



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "6_estimate_gwas_genetic_correlation_ldsc_secondary.sh"
  echo "----------"
fi



#
