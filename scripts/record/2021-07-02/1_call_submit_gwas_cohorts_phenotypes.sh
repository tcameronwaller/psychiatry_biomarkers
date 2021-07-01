#!/bin/bash

#chmod u+x script.sh

###########################################################################
###########################################################################
###########################################################################
# This script organizes directories and iteration instances then submits
# script "regress_metabolite_heritability.sh" to the Sun Grid Engine.

# version check: 2

###########################################################################
###########################################################################
###########################################################################

# Organize paths.
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-07-02"
path_dock="$path_process/dock"
path_cohorts_models="${path_dock}/organization/cohorts_models"
path_gwas="${path_dock}/gwas/cohorts_models_maf_0" # run a few GWAS with minimal MAF 0

# Initialize directories.
#rm -r $path_gwas
mkdir -p $path_gwas

# Define covariates common for all cohorts.
covariates_common="genotype_pc_1,genotype_pc_2,genotype_pc_3,genotype_pc_4,genotype_pc_5,genotype_pc_6,genotype_pc_7,genotype_pc_8,genotype_pc_9,genotype_pc_10"

# Define multi-dimensional array of cohorts and covariates.

# Note:
# Each GWAS (30,000 - 200,000 persons; 22 chromosomes) requires about 5-7 hours to run on the grid.

cohorts_models=()

#cohorts_models+=("all_bipolar_disorder_case;table_all_bipolar_disorder_case;sex,age,")
cohorts_models+=("white_bipolar_disorder_case;table_white_bipolar_disorder_case;sex,age,")
#cohorts_models+=("all_bipolar_disorder_control;table_all_bipolar_disorder_control;sex,age,")
cohorts_models+=("white_bipolar_disorder_control;table_white_bipolar_disorder_control;sex,age,")

#cohorts_models+=("all_bipolar_disorder_case_simple;table_all_bipolar_disorder_case;sex,")
cohorts_models+=("white_bipolar_disorder_case_simple;table_white_bipolar_disorder_case;sex,")
#cohorts_models+=("all_bipolar_disorder_control_simple;table_all_bipolar_disorder_control;sex,")
cohorts_models+=("white_bipolar_disorder_control_simple;table_white_bipolar_disorder_control;sex,")

#cohorts_models+=("all_bipolar_disorder_case_unadjust;table_all_bipolar_disorder_case;")
cohorts_models+=("white_bipolar_disorder_case_unadjust;table_white_bipolar_disorder_case;")
#cohorts_models+=("all_bipolar_disorder_control_unadjust;table_all_bipolar_disorder_control;")
cohorts_models+=("white_bipolar_disorder_control_unadjust;table_white_bipolar_disorder_control;")


# Define array of hormones.
phenotypes=()
phenotypes+=("body_mass_index") # <- start on 1 July 2021
#phenotypes+=("body_mass_index_log")

# Assemble array of batch instance details.
path_batch_instances="${path_gwas}/batch_instances.txt"
rm $path_batch_instances
for cohort_model in "${cohorts_models[@]}"; do
  for phenotype in "${phenotypes[@]}"; do
    instance="${phenotype};${cohort_model}${covariates_common}"
    echo $instance
    echo $instance >> $path_batch_instances
  done
done

# Read batch instances.
readarray -t batch_instances < $path_batch_instances
batch_instances_count=${#batch_instances[@]}
echo "----------"
echo "count of batch instances: " $batch_instances_count
echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
echo "last batch instance: " ${batch_instances[batch_instances_count - 1]}

if true; then
  # Submit array batch to Sun Grid Engine.
  # Array batch indices must start at one (not zero).
  echo "----------------------------------------------------------------------"
  echo "Submit array of batches to Sun Grid Engine."
  echo "----------------------------------------------------------------------"
  qsub -t 1-${batch_instances_count}:1 \
  -o "${path_gwas}/out.txt" -e "${path_gwas}/error.txt" \
  "${path_scripts_record}/2_organize_call_run_chromosomes_plink_gwas.sh" \
  $path_batch_instances \
  $batch_instances_count \
  $path_cohorts_models \
  $path_gwas \
  $path_scripts_record
fi
