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

# TODO: I need to access the "cohorts_models" tables from the "stratification" directory...


# Organize paths.
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2021-10-25"
path_dock="$path_process/dock"
path_cohorts_models="${path_dock}/stratification_2021-10-21/cohorts_models_linear"

#path_gwas="${path_dock}/gwas/body_white_bipolar_strict"          # 12 GWAS; TCW started at 13:41 on 22 October 2021
path_gwas="${path_dock}/gwas/body_white_bipolar_loose"          # __ GWAS; TCW started at ___ on ____ 2021

# Initialize directories.
#rm -r $path_gwas
mkdir -p $path_gwas

##########
##########
##########
# Assemble a list of analysis instances with common patterns.

# Assemble array of batch instance details.
path_batch_instances="${path_gwas}/batch_instances.txt"
rm $path_batch_instances

##########
# General models.

# Define covariates common for all cohorts.
covariates_common="genotype_pc_1,genotype_pc_2,genotype_pc_3,genotype_pc_4,genotype_pc_5,genotype_pc_6,genotype_pc_7,genotype_pc_8,genotype_pc_9,genotype_pc_10"
# Define multi-dimensional array of cohorts and model covariates.
cohorts_models=()

### body_white_bipolar_strict
#cohorts_models+=("white_bipolar_strict_control_unadjust;table_white_bipolar_strict_control;")
#cohorts_models+=("white_bipolar_strict_control_sex;table_white_bipolar_strict_control;sex,")
#cohorts_models+=("white_bipolar_strict_control_sex_age;table_white_bipolar_strict_control;sex,age,")
#cohorts_models+=("white_bipolar_strict_case_unadjust;table_white_bipolar_strict_case;")
#cohorts_models+=("white_bipolar_strict_case_sex;table_white_bipolar_strict_case;sex,")
#cohorts_models+=("white_bipolar_strict_case_sex_age;table_white_bipolar_strict_case;sex,age,")

### body_white_bipolar_loose
cohorts_models+=("white_bipolar_loose_control_unadjust;table_white_bipolar_loose_control;")
cohorts_models+=("white_bipolar_loose_control_sex;table_white_bipolar_loose_control;sex,")
cohorts_models+=("white_bipolar_loose_control_sex_age;table_white_bipolar_loose_control;sex,age,")
cohorts_models+=("white_bipolar_loose_case_unadjust;table_white_bipolar_loose_case;")
cohorts_models+=("white_bipolar_loose_case_sex;table_white_bipolar_loose_case;sex,")
cohorts_models+=("white_bipolar_loose_case_sex_age;table_white_bipolar_loose_case;sex,age,")

# Define array of phenotypes.
phenotypes=()

phenotypes+=("body")
phenotypes+=("body_log")

for cohort_model in "${cohorts_models[@]}"; do
  for phenotype in "${phenotypes[@]}"; do
    instance="${phenotype};${cohort_model}${covariates_common}"
    echo $instance
    echo $instance >> $path_batch_instances
  done
done

##########
# Summary of batch instances.

# Array pattern.
#phenotype="${array[0]}"
#cohort_model="${array[1]}"
#table_cohort_model="${array[2]}"
#covariates="${array[3]}"
#path_table_phenotypes_covariates="${path_cohorts_models}/${table_cohort_model}_${phenotype}.tsv"

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
  "${path_scripts_record}/4-1_organize_call_run_chromosomes_plink_linear_gwas.sh" \
  $path_batch_instances \
  $batch_instances_count \
  $path_cohorts_models \
  $path_gwas \
  $path_scripts_record
fi
