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
path_temporary=$(<"./processing_sexy_alcohol.txt")
path_waller="$path_temporary/waller"
path_scripts_record="$path_waller/sexy_alcohol/scripts/record/2021-04-01"
path_dock="$path_temporary/waller/dock"
path_cohorts="${path_dock}/organization/cohorts"
path_gwas="${path_dock}/gwas"

# Initialize directories.
#rm -r $path_gwas
mkdir -p $path_gwas

##########
# Cohorts: covariates
# female_male: "sex", "age", "body_mass_index_log", "genotype_pc_[1-10]"
# female: "age", "body_mass_index_log", "menopause_hormone_category_[???]", "genotype_pc_[1-10]"
# female_premenopause: "age", "body_mass_index_log", "menstruation_day", "hormone_alteration", "genotype_pc_[1-10]"
# female_postmenopause: "age", "body_mass_index_log", "hormone_alteration", "genotype_pc_[1-10]"
# male: "age", "body_mass_index_log", "genotype_pc_[1-10]"

##########
# Phenotypes
# "oestradiol_log" (batch almost complete)
# "oestradiol_free_log" (batch submitted at 22:20)
# "testosterone_log"
# "testosterone_free_log"

if true; then

  phenotype="testosterone_free_log" # "oestradiol_log", "oestradiol_free_log", "testosterone_log", "testosterone_free_log"

  # Jobs: ___, ___
  # Parameters.
  file_name="table_female_male_${phenotype}.tsv"
  cohort_phenotype="female_male_${phenotype}"
  #phenotype="oestradiol_free_log"
  covariates="sex,age,body_mass_index_log,genotype_pc_1,genotype_pc_2,genotype_pc_3,genotype_pc_4,genotype_pc_5,genotype_pc_6,genotype_pc_7,genotype_pc_8,genotype_pc_9,genotype_pc_10"
  /usr/bin/bash "$path_scripts_record/2_submit_gwas_cohorts_phenotypes.sh" \
  $file_name \
  $cohort_phenotype \
  $phenotype \
  $covariates \
  $path_scripts_record \
  $path_cohorts \
  $path_gwas

  # Jobs: ___, ___
  # Parameters.
  file_name="table_female_${phenotype}.tsv"
  cohort_phenotype="female_${phenotype}"
  covariates="age,body_mass_index_log,menopause_hormone_category_1,menopause_hormone_category_3,menopause_hormone_category_4,genotype_pc_1,genotype_pc_2,genotype_pc_3,genotype_pc_4,genotype_pc_5,genotype_pc_6,genotype_pc_7,genotype_pc_8,genotype_pc_9,genotype_pc_10"
  /usr/bin/bash "$path_scripts_record/2_submit_gwas_cohorts_phenotypes.sh" \
  $file_name \
  $cohort_phenotype \
  $phenotype \
  $covariates \
  $path_scripts_record \
  $path_cohorts \
  $path_gwas

  # Jobs: ___, ___
  # Parameters.
  file_name="table_female_premenopause_${phenotype}.tsv"
  cohort_phenotype="female_premenopause_${phenotype}"
  covariates="age,body_mass_index_log,menstruation_day,hormone_alteration,genotype_pc_1,genotype_pc_2,genotype_pc_3,genotype_pc_4,genotype_pc_5,genotype_pc_6,genotype_pc_7,genotype_pc_8,genotype_pc_9,genotype_pc_10"
  /usr/bin/bash "$path_scripts_record/2_submit_gwas_cohorts_phenotypes.sh" \
  $file_name \
  $cohort_phenotype \
  $phenotype \
  $covariates \
  $path_scripts_record \
  $path_cohorts \
  $path_gwas

  # Jobs: ___, ___
  # Parameters.
  file_name="table_female_postmenopause_${phenotype}.tsv"
  cohort_phenotype="female_postmenopause_${phenotype}"
  covariates="age,body_mass_index_log,hormone_alteration,genotype_pc_1,genotype_pc_2,genotype_pc_3,genotype_pc_4,genotype_pc_5,genotype_pc_6,genotype_pc_7,genotype_pc_8,genotype_pc_9,genotype_pc_10"
  /usr/bin/bash "$path_scripts_record/2_submit_gwas_cohorts_phenotypes.sh" \
  $file_name \
  $cohort_phenotype \
  $phenotype \
  $covariates \
  $path_scripts_record \
  $path_cohorts \
  $path_gwas

  # Jobs: ___, ___
  # Parameters.
  file_name="table_male_${phenotype}.tsv"
  cohort_phenotype="male_${phenotype}"
  covariates="age,body_mass_index_log,genotype_pc_1,genotype_pc_2,genotype_pc_3,genotype_pc_4,genotype_pc_5,genotype_pc_6,genotype_pc_7,genotype_pc_8,genotype_pc_9,genotype_pc_10"
  /usr/bin/bash "$path_scripts_record/2_submit_gwas_cohorts_phenotypes.sh" \
  $file_name \
  $cohort_phenotype \
  $phenotype \
  $covariates \
  $path_scripts_record \
  $path_cohorts \
  $path_gwas

fi
