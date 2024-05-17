#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 17 May 2024
# Date, last execution: 17 May 2024
# Date, review: __ May 2024
################################################################################
# Note

# SLURM batch job: ___ (group: "all"; instances: ___; date: 17 May 2024)

##########
# Note: TCW; 17 May 2024
# The file "table_gwas_translation_tcw_2023-12-30_dbsnp_rsid.tsv" provided
# parameters for the preparation and organization of a collection of 244 files
# of GWAS summary statistics. Of these 244 files of GWAS summary statistics, 80
# 'primaries' were for psychiatric disorders or substance use disorders or their
# related traits, and 164 'secondaries' were for thyroid disorders, thyroid
# biomarkers, or other biomarkers including sex hormones and their binding
# proteins. On 8 January 2024, TCW ran LDSC munge on all 244 files of GWAS
# summary statistics. On 8 and 9 January 2024, TCW ran LDSC SNP heritability
# estimates on all 244 munged files of GWAS summary statistics. On 16 May 2024,
# TCW ran LDSC genetic correlation estimates on pairwise combinations between
# all 244 sets of GWAS summary statistics in the collection.

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
path_directory_source="${path_directory_group_parent}/4_gwas_munge_ldsc"
path_directory_product="${path_directory_group_parent}/6_gwas_correlation_ldsc_all"
path_directory_disequilibrium="${path_directory_reference}/disequilibrium/eur_w_ld_chr"
path_directory_batch="${path_directory_product_parent}/batch"

# Files.
path_file_table_parameter="${path_directory_parameters}/table_gwas_translation_${identifier_parameter}.tsv"

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
rm -r $path_directory_product # caution
rm -r $path_directory_batch # caution
mkdir -p $path_directory_product
mkdir -p $path_directory_batch

# Initialize files.
rm $path_file_batch_instances

################################################################################
# Organize parameters.

##########
# Common parameters.
threads=2
report="false"

##########
# Extract identifiers of studies for which to estimate genetic correlations.
# Review: TCW; __ May 2024
# Count of all studies and versions: 244 (244 X 244 = 59,536)

studies=()

# Read lines from file and split fields within each line by space, tab, or new-line delimiters.
input=$path_file_table_parameter
while IFS=$' \t\n' read -r -a array
do

  # Extract values from individual columns within table's row.
  # These parameters match the format of the table as of 13 December 2023.
  raw_availability="${array[0]}"
  raw_inclusion="${array[1]}"
  raw_directory="${array[2]}"
  raw_name_study="${array[3]}"
  raw_phenotype="${array[4]}"
  raw_sex="${array[5]}"
  raw_name_file_source="${array[6]}"
  raw_suffix_file_source="${array[7]}"
  raw_bgzip="${array[8]}"
  raw_gzip="${array[9]}"
  raw_type="${array[10]}"
  raw_fill_observations="${array[11]}"
  raw_observations_total="${array[12]}"
  raw_fill_case_control="${array[13]}"
  raw_cases="${array[14]}"
  raw_controls="${array[15]}"
  raw_observations_effective="${array[16]}"
  raw_prevalence_sample="${array[17]}"
  raw_prevalence_population="${array[18]}"
  raw_script="${array[19]}"
  raw_note="${array[20]}"

  # Report.
  if [ $raw_inclusion == "1" ] && [ "$report" == "true" ]; then
    echo "----------"
    echo "field 0, availability: ${raw_availability}"
    echo "field 1, inclusion: ${raw_inclusion}"
    echo "field 2, directory: ${raw_directory}"
    echo "field 3, name_study: ${raw_name_study}"
    echo "field 4, phenotype: ${raw_phenotype}"
    echo "field 5, sex: ${raw_sex}"
    echo "field 6, file: ${raw_name_file_source}"
    echo "field 7, suffix: ${raw_suffix_file_source}"
    echo "field 8, bgzip: ${raw_bgzip}"
    echo "field 9, gzip: ${raw_gzip}"
    echo "field 10, type: ${raw_type}"
    echo "field 11, fill_observations: ${raw_fill_observations}"
    echo "field 12, observations_total: ${raw_observations_total}"
    echo "field 13, fill_case_control: ${raw_fill_case_control}"
    echo "field 14, cases: ${raw_cases}"
    echo "field 15, controls: ${raw_controls}"
    echo "field 16, observations_effective: ${raw_observations_effective}"
    echo "field 17, prevalence_sample: ${raw_prevalence_sample}"
    echo "field 18, prevalence_population: ${raw_prevalence_population}"
    echo "field 19, script: ${raw_script}"
    echo "field 20, note: ${raw_note}"
    echo "----------"
  fi
  # Execute procedure for current record's parameters.
  if [ $raw_inclusion == "1" ]; then
    # Collection identifier of study.
    studies+=("${raw_name_study}")
  fi
done < "${input}"

##########
# Common parameters.
threads=2
report="true"

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
  for primary in "${studies[@]}"; do
    for secondary in "${studies[@]}"; do
      # Organize paths.
      name_comparison="${primary}_-_${secondary}"
      path_file_base_product="${path_directory_product}/${name_comparison}"
      path_file_source_primary="${path_directory_source}/${primary}.sumstats.gz"
      path_file_source_secondary="${path_directory_source}/${secondary}.sumstats.gz"
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

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "6_4_estimate_gwas_genetic_correlation_ldsc_all.sh"
  echo "----------"
fi



#
