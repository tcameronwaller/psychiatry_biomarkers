#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 14 December 2023
# Date, last execution: 22 December 2023
# Date, review: 22 December 2023
################################################################################
# Note

# Note: TCW; 23 December 2023
# For the 40 studies below I submitted batch job array 3178183 to the
# computational grid cluster via the Slurm Workload Manager on 22 December 2023.
# On 23 December 2023, the batch job array completed successfully (exit code 0),
# using up to 277.26 Gigabytes of memory for each job.

# Note: TCW; 29 December 2023
# For example, here are the results of the "fill_dbsnp_rsid.sh" script procedure
# for the study "37872160_williams_2023" as recorded in the file below.
# "/.../gwas_2023-12-22_dbsnp_rsid/5_fill_dbsnp_rs_identifiers/batch/..."
# ".../fill_rsid_3178183.mforgehi2.3178207.0.stdout"
# Count of lines in original source file: 57,523,911
# Count of lines in novel product file: 57,523,911
# Notice that the count of lines was the same since parameter "strict" was set
# to "false".
# With parameter "strict" set to false, the procedure ought to be loss-less,
# preserving all information from the original GWAS summary statistics.
# Count of lines that matched and merged with dbSNP: 36,107,225
# Percentage of lines that matched and merged with dbSNP: 62.769%

# Note:


################################################################################
# Organize paths.

# Identifiers or designators of parameter version and preparation batch.
identifier_preparation="gwas_2023-12-22_dbsnp_rsid"
identifier_parameter="tcw_2023-12-22_dbsnp_rsid"

# Directories.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes

path_directory_source="${path_directory_dock}/${identifier_preparation}/4_filter_constrain_gwas_values"
path_directory_product="${path_directory_dock}/${identifier_preparation}/5_fill_dbsnp_rs_identifiers"
path_directory_batch="${path_directory_product}/batch"

# Files.
path_file_batch_instances="${path_directory_batch}/batch_instances.txt"

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_file_script_slurm_job="${path_directory_partner_scripts}/gwas_clean/slurm_job_fill_dbsnp_rsid.sh"
path_file_script_fill_dbsnp_rsid="${path_directory_partner_scripts}/gwas_clean/fill_dbsnp_rsid.sh"


# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
#cd $path_directory_product
rm -r $path_directory_batch # caution
mkdir -p $path_directory_batch
cd $path_directory_batch # execute batch from within this directory

# Initialize files.
rm $path_file_batch_instances



###########################################################################
# Organize parameters.

strict="false"
report="true"

################################################################################
# Execute procedure.



##########
# Copy the GWAS summary statistics from the previous process.
# Most sets of GWAS summary statistics do not need extra processing.
# Subsequent processes on a few studies will replace the appropriate files.
#cp $path_directory_source/*.txt.gz $path_directory_product

##########
# Fill missing SNP rsIDs from dbSNP.

# Define array of studies.
studies=()

studies+=("37872160_williams_2023")
studies+=("36635386_chen_2023_cortisol") # Unnecessary. Supplement.
studies+=("36635386_chen_2023_thyroxine_total") # Unnecessary. Supplement.
studies+=("34662886_backman_2021_albumin")

studies+=("34594039_sakaue_2021_multi_hypothyroidism")
studies+=("34594039_sakaue_2021_multi_hyperthyroidism")
studies+=("34594039_sakaue_2021_multi_hashimoto")
studies+=("34594039_sakaue_2021_multi_graves")
studies+=("34594039_sakaue_2021_eur_hypothyroidism")
studies+=("34594039_sakaue_2021_eur_hyperthyroidism")
studies+=("34594039_sakaue_2021_eur_hashimoto")
studies+=("34594039_sakaue_2021_eur_graves")
studies+=("34594039_sakaue_2021_eur_rheumatoid_arthritis")
studies+=("34594039_sakaue_2021_gc_hypothyroidism")
studies+=("34594039_sakaue_2021_gc_hyperthyroidism")
studies+=("34594039_sakaue_2021_gc_hashimoto")
studies+=("34594039_sakaue_2021_gc_graves")

studies+=("34226706_barton_2021_albumin")
studies+=("34017140_mbatchou_2021_albumin")

studies+=("33587031_sinnott-armstrong_2021_testosterone_primary_female")
studies+=("33587031_sinnott-armstrong_2021_testosterone_primary_male")
studies+=("33587031_sinnott-armstrong_2021_lutropin")

studies+=("32769997_zhou_2020_thyroid_hormone")
studies+=("32581359_saevarsdottir_2020_thyroid_autoimmunity") # Unnecessary. Supplement.
studies+=("32059762_manousaki_2020_vitamin_d")

studies+=("30367059_teumer_2018_thyroid_hormone_all")
studies+=("30367059_teumer_2018_thyroid_hormone_female")
studies+=("30367059_teumer_2018_thyroid_hormone_male")
studies+=("30367059_teumer_2018_thyroxine_free_all")
studies+=("30367059_teumer_2018_thyroxine_free_female")
studies+=("30367059_teumer_2018_thyroxine_free_male")
studies+=("30367059_teumer_2018_hypothyroidism")
studies+=("30367059_teumer_2018_hyperthyroidism")

studies+=("29875488_sun_2018_complement_c4")
studies+=("29875488_sun_2018_follitropin")
studies+=("29875488_sun_2018_follistatin")
studies+=("29875488_sun_2018_lutropin")
studies+=("29875488_sun_2018_lutropin_beta")
studies+=("29875488_sun_2018_parathyrin")
studies+=("29875488_sun_2018_thyroid_peroxidase")


##########
# Execute process iteratively for each study.
if false; then
  for study in "${studies[@]}"; do
    path_file_source="${path_directory_source}/${study}.txt.gz"
    path_file_product="${path_directory_product}/${study}.txt.gz"
    # Call script.
    /usr/bin/bash $path_file_script_fill_dbsnp_rsid \
    $path_file_source \
    $path_file_product \
    $strict \
    $report
  done
fi


##########
# Organize array of instance parameters.

instances=()
# instances+=("${part_one};${part_two};${part_three};${part_four}")

# Assemble array of batch instance details.
for study in "${studies[@]}"; do
  path_file_source="${path_directory_source}/${study}.txt.gz"
  path_file_product="${path_directory_product}/${study}.txt.gz"
  # Assemble parameters for batch instance.
  instances+=("${path_file_source};${path_file_product}")
done



##########
# Organize batch instances.

count_instances=${#instances[@]}

for instance in "${instances[@]}"; do
  # Define parameters in array instance for batch job.
  echo $instance >> $path_file_batch_instances
done

# Read batch instances.
readarray -t batch_instances < $path_file_batch_instances
batch_instances_count=${#batch_instances[@]}
index_array_maximum=$(($batch_instances_count - 1))



################################################################################
# Report.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Source directory:"
  echo $path_directory_source
  echo "Product directory:"
  echo $path_directory_product
  echo "count of instances: " $count_instances
  echo "first instance: " ${instances[0]} # notice base-zero indexing
  echo "last instance: " ${instances[$count_instances - 1]}
  echo "----------"
  echo "----------"
  echo "count of batch instances: " $batch_instances_count
  echo "maximum array index: " $index_array_maximum
  echo "----------"
  echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
  echo "last batch instance: " ${batch_instances[$index_array_maximum]}
  echo "----------"

fi

sleep 5s



################################################################################
# Submit batch of jobs to grid cluster scheduler for processing.
# Submit to Slurm Scheduler.
# Indices in array of batch jobs start at zero.

if true; then
  sbatch --array 0-${index_array_maximum}:1 --chdir $path_directory_batch \
  $path_file_script_slurm_job \
  $path_file_batch_instances \
  $batch_instances_count \
  $strict \
  $report \
  $path_file_script_fill_dbsnp_rsid
fi



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "5_fill_dbsnp_rs_identifiers.sh"
  echo "----------"
fi



#
