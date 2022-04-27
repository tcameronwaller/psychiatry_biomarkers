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
path_scripts_record="$path_process/psychiatric_metabolism/scripts/record/2022-05-04/gwas_association"
path_dock="$path_process/dock"

# Note: TCW, 02 March 2022
# Logistic GWAS require a long time to run.
# Logistic GWAS for "oestradiol_detection" on the "female",
# "female_postmenopause", and "male" cohorts have very large sample size and
# require much longer than 30 hours to run.
# The "1-day" queue on NCSA server aborts batch jobs after about 30 hours.

################
# TODO: TCW, 08 April 2022
# TODO: submit "oestradiol_detection" to the 4-day queue at least... maybe 7-day queue...
####################


# Cohort stratification phenotype tables.

path_stratification_tables="${path_dock}/stratification_2022-04-09/oestradiol" # 09 April 2022, 10 April 2022, 11 April 2022
#path_stratification_tables="${path_dock}/stratification_2022-04-09/testosterone" # 09 April 2022, 10 April 2022, 11 April 2022
#path_stratification_tables="${path_dock}/stratification_2022-04-09/steroid_globulin" # 11 April 2022
#path_stratification_tables="${path_dock}/stratification_2022-04-09/albumin" # 11 April 2022

# Container directory for raw GWAS summary statistics.

#path_gwas_container="${path_dock}/gwas_raw/oestradiol_logistic"              # 18 GWAS; 11 April 2022
path_gwas_container="${path_dock}/gwas_raw/oestradiol_logistic_long_rescue"       # __ GWAS; 27 April 2022

#path_gwas_container="${path_dock}/gwas_raw/oestradiol_linear_1"              # 18 GWAS; 09 April 2022
#path_gwas_container="${path_dock}/gwas_raw/oestradiol_linear_2"              # 12 GWAS; 11 April 2022
#path_gwas_container="${path_dock}/gwas_raw/oestradiol_bioavailable_linear"   # 18 GWAS; 11 April 2022
#path_gwas_container="${path_dock}/gwas_raw/oestradiol_free_linear"           # 18 GWAS; 11 April 2022

#path_gwas_container="${path_dock}/gwas_raw/testosterone_logistic"            # 18 GWAS; 11 April 2022
#path_gwas_container="${path_dock}/gwas_raw/testosterone_logistic_long_rescue"     # __ GWAS; 27 April 2022

#path_gwas_container="${path_dock}/gwas_raw/testosterone_linear"              # 54 GWAS; 09 April 2022
#path_gwas_container="${path_dock}/gwas_raw/testosterone_bioavailable_linear" # 18 GWAS; 10 April 2022
#path_gwas_container="${path_dock}/gwas_raw/testosterone_free_linear"         # 18 GWAS; 10 April 2022

#path_gwas_container="${path_dock}/gwas_raw/steroid_globulin_linear_1"        # 18 GWAS; 11 April 2022
#path_gwas_container="${path_dock}/gwas_raw/steroid_globulin_linear_2"        # 2 GWAS; 11 April 2022
#path_gwas_container="${path_dock}/gwas_raw/albumin_linear_1"                 # 18 GWAS; 11 April 2022
#path_gwas_container="${path_dock}/gwas_raw/albumin_linear_2"                 # 2 GWAS; 11 April 2022

# Initialize directories.
rm -r $path_gwas_container
mkdir -p $path_gwas_container

##########
##########
##########
# Assemble a list of analysis instances with common patterns.

# Assemble array of batch instance details.
path_batch_instances="${path_gwas_container}/batch_instances.txt"
rm $path_batch_instances

##########
# General models.

##########
# Define covariates common for all cohorts.
covariates_genotype="genotype_array_axiom,genotype_pc_1,genotype_pc_2,genotype_pc_3,genotype_pc_4,genotype_pc_5,genotype_pc_6,genotype_pc_7,genotype_pc_8,genotype_pc_9,genotype_pc_10"

##########
# Basis covariates that are the same across all cohorts

# model: 'unadjust'
covariates_unadjust=""

# model: 'joint_1' # Model 'joint_1' is the same for all dependent variables for Sex Hormones and their Binding Proteins (not for Vitamin D).
covariates_joint_1="age,body,region,season,alteration_sex_hormone"

##########
# Define multi-dimensional array of cohorts and model covariates.
# Use comments to exclude cohorts from GWAS for a specific dependent variable.
# [name of cohort and model for analysis description];[table name cohort-model prefix];[independent variable columns in table, beginning with cohort-specific variables]
cohorts_models_instances=()

# cohort: female_male_priority_male: only for steroid_globulin_imputation, albumin_imputation
#cohorts_models_instances+=("female_male_priority_male_unadjust;table_female_male_priority_male_unadjust;sex_y,") # include 'sex_y' at minimum for X chromosome
#cohorts_models_instances+=("female_male_priority_male_joint_1;table_female_male_priority_male_joint_1;sex_y,${covariates_joint_1},")

# cohort: female
#cohorts_models_instances+=("female_unadjust;table_female_unadjust;")
cohorts_models_instances+=("female_joint_1;table_female_joint_1;births,menopause_ordinal,oophorectomy,hysterectomy,${covariates_joint_1},")

# cohort: female_menstruation_regular
#cohorts_models_instances+=("female_menstruation_regular_unadjust;table_female_menstruation_regular_unadjust;")
#cohorts_models_instances+=("female_menstruation_regular_joint_1;table_female_menstruation_regular_joint_1;menstruation_phase_cycle,birth_live_recent,births,menopause_ordinal,oophorectomy,hysterectomy,${covariates_joint_1},")

# cohort: female_premenopause
#cohorts_models_instances+=("female_premenopause_unadjust;table_female_premenopause_unadjust;")
#cohorts_models_instances+=("female_premenopause_joint_1;table_female_premenopause_joint_1;menstruation_phase_cycle,birth_live_recent,births,hysterectomy,${covariates_joint_1},")

# cohort: female_perimenopause
#cohorts_models_instances+=("female_perimenopause_unadjust;table_female_perimenopause_unadjust;")
#cohorts_models_instances+=("female_perimenopause_joint_1;table_female_perimenopause_joint_1;births,hysterectomy,${covariates_joint_1},")

# cohort: female_postmenopause
#cohorts_models_instances+=("female_postmenopause_unadjust;table_female_postmenopause_unadjust;")
#cohorts_models_instances+=("female_postmenopause_joint_1;table_female_postmenopause_joint_1;births,oophorectomy,hysterectomy,${covariates_joint_1},")

# cohort: male
#cohorts_models_instances+=("male_unadjust;table_male_unadjust;")
cohorts_models_instances+=("male_joint_1;table_male_joint_1;${covariates_joint_1},")

# cohort: male_age_low
#cohorts_models_instances+=("male_age_low_unadjust;table_male_age_low_unadjust;")
#cohorts_models_instances+=("male_age_low_joint_1;table_male_age_low_joint_1;${covariates_joint_1},")

# cohort: male_age_middle
#cohorts_models_instances+=("male_age_middle_unadjust;table_male_age_middle_unadjust;")
#cohorts_models_instances+=("male_age_middle_joint_1;table_male_age_middle_joint_1;${covariates_joint_1},")

# cohort: male_age_high
#cohorts_models_instances+=("male_age_high_unadjust;table_male_age_high_unadjust;")
#cohorts_models_instances+=("male_age_high_joint_1;table_male_age_high_joint_1;${covariates_joint_1},")

##########
# Define array of phenotypes.
# [name of phenotype for analysis description];[table name phenotype suffix];[dependent variable column in table]
phenotypes_instances=()

phenotypes_instances+=("oestradiol_detection;_oestradiol_detection;oestradiol_detection_plink")                                                         # 'oestradiol_logistic'; 11 April 2022
#phenotypes_instances+=("oestradiol_log;_oestradiol_log;oestradiol_log")                                                                                 # 'oestradiol_linear_1'; 09 April 2022
#phenotypes_instances+=("oestradiol_imputation;_oestradiol_imputation;oestradiol_imputation")                                                            # 'oestradiol_linear_2'; 11 April 2022
#phenotypes_instances+=("oestradiol_imputation_log;_oestradiol_imputation_log;oestradiol_imputation_log")                                                # 'oestradiol_linear_2'; 11 April 2022
#phenotypes_instances+=("oestradiol_bioavailable_imputation_log;_oestradiol_bioavailable_imputation_log;oestradiol_bioavailable_imputation_log")         # 'oestradiol_bioavailable_linear'; 11 April 2022
#phenotypes_instances+=("oestradiol_free_imputation_log;_oestradiol_free_imputation_log;oestradiol_free_imputation_log")                                 # 'oestradiol_free_linear'; 11 April 2022

#phenotypes_instances+=("testosterone_detection;_testosterone_detection;testosterone_detection_plink")                                                   # 'testosterone_logistic'; 11 April 2022
#phenotypes_instances+=("testosterone_log;_testosterone_log;testosterone_log")                                                                           # 'testosterone_linear'; 09 April 2022
#phenotypes_instances+=("testosterone_imputation;_testosterone_imputation;testosterone_imputation")                                                      # 'testosterone_linear'; 09 April 2022
#phenotypes_instances+=("testosterone_imputation_log;_testosterone_imputation_log;testosterone_imputation_log")                                          # 'testosterone_linear'; 09 April 2022
#phenotypes_instances+=("testosterone_bioavailable_imputation_log;_testosterone_bioavailable_imputation_log;testosterone_bioavailable_imputation_log")   # 'testosterone_bioavailable_linear'; 10 April 2022
#phenotypes_instances+=("testosterone_free_imputation_log;_testosterone_free_imputation_log;testosterone_free_imputation_log")                           # 'testosterone_free_linear'; 10 April 2022

#phenotypes_instances+=("steroid_globulin_imputation_log;_steroid_globulin_imputation_log;steroid_globulin_imputation_log")                              # 'steroid_globulin_linear_1', 'steroid_globulin_linear_2'; 11 April 2022
#phenotypes_instances+=("albumin_imputation;_albumin_imputation;albumin_imputation")                                                                     # 'albumin_linear_1', 'albumin_linear_2'; 11 April 2022



##########
# Assemble details for batch instances.
for cohort_model_instance in "${cohorts_models_instances[@]}"; do
  # Separate fields from instance.
  IFS=";" read -r -a array <<< "${cohort_model_instance}"
  cohort_model_name="${array[0]}" # name of cohort and model for analysis description
  table_cohort_prefix="${array[1]}" # table name prefix specific to cohort and model
  covariates_specific="${array[2]}" # names of columns in table for covariate variables specific to analysis
  for phenotype_instance in "${phenotypes_instances[@]}"; do
    # Separate fields from instance.
    IFS=";" read -r -a array <<< "${phenotype_instance}"
    phenotype_name="${array[0]}" # name of phenotype for unique analysis description
    table_phenotype_suffix="${array[1]}" # phenotype specific suffix to identify table
    phenotypes="${array[2]}" # names of columns in table for phenotype variables
    # Organize variables.
    name_study="${cohort_model_name}_${phenotype_name}" # unique name for study analysis for parent directory
    table_phenotypes_covariates="${table_cohort_prefix}${table_phenotype_suffix}.tsv" # name of table of phenotypes and covariates
    covariates="${covariates_specific}${covariates_genotype}" # names of columns in table for covariate variables
    # Assemble details for batch instance.
    instance="${name_study};${table_phenotypes_covariates};${phenotypes};${covariates}"
    echo $instance
    echo $instance >> $path_batch_instances
  done
done

##########
# Summary of batch instances.

# Array pattern.
#name_study="${array[0]}"
#table_phenotypes_covariates="${array[1]}"
#phenotypes="${array[2]}"
#covariates="${array[3]}"

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
  -o "${path_gwas_container}/out.txt" -e "${path_gwas_container}/error.txt" \
  "${path_scripts_record}/2-1_organize_call_run_gwas_chromosomes_plink_association.sh" \
  $path_batch_instances \
  $batch_instances_count \
  $path_stratification_tables \
  $path_gwas_container \
  $path_scripts_record \
  $path_process
fi
