#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 2 October 2023
# Date, last execution: 19 December 2023
# Date, review: 19 December 2023
################################################################################
# Note

# ???
# Study: 32581359_saevarsdottir_2020_thyroid_autoimmunity
# TODO: TCW; 22 November 2023 --> Does this study still fail GWAS2VCF? I don't think so.
# This study failed GWAS2VCF procedure for unknown reasons.
# Rescue this study by copying set of GWAS summary statistics from step 4 of
# preparation procedure before GWAS2VCF.

# Note: TCW; 4 December 2023
# After translation from GRCh38 to GRCh37, the study "37872160_williams_2023"
# completed the GWAS2VCF procedure with 10,152,437 SNPs.

# Note: TCW; 28 November 2023
# The following studies failed the GWAS2VCF procedure.
# Study: 37872160_williams_2023 (failed in batch submitted on 27 November 2023)
# Study: 32581359_saevarsdottir_2020_thyroid_autoimmunity (failed in batch submitted on 27 November 2023)
# Study: 32242144_revez_2020_vitamin_d
# Study: 32099098_polimanti_2020_eur_opioid_dep_exposed
# Study: 32099098_polimanti_2020_eur_opioid_dep_unexposed
# Study: 32099098_polimanti_2020_eur_opioid_exposure
# Study: 30718901_howard_2019_pgc_ukb

# The following studies failed GWAS2VCF procedure because the original GWAS
# summary statistics did not include information about chromosome and base-pair
# coordinates. This information is not necessary for analysis in LDSC.
# Rescue these studies by copying sets of GWAS summary statistics from step 4 of
# preparation procedure before GWAS2VCF.
# Study: 32242144_revez_2020_vitamin_d
# Study: 32099098_polimanti_2020_eur_opioid_dep_exposed
# Study: 32099098_polimanti_2020_eur_opioid_dep_unexposed
# Study: 32099098_polimanti_2020_eur_opioid_exposure
# Study: 30718901_howard_2019_pgc_ukb

# Note: TCW; 28 November 2023
# The following studies had problems previously, but the filters before GWAS2VCF
# might have resolved these problems adequately.
# 36376304_koskeridis_2022_c_reactive_protein
# 34662886_backman_2021_albumin
# 34226706_barton_2021_albumin
# 32769997_zhou_2020_thyroid_hormone
# 32581359_saevarsdottir_2020_thyroid_autoimmunity

# Note: TCW; 28 November 2023
# The conditionals in the "check" script could be more strict than the
# conditionals in the "filter" script for the sake of describing potential
# problems within the sets of GWAS summary statistics.

################################################################################
# Organize paths.

# Identifiers or designators of parameter version and preparation batch.
identifier_preparation="gwas_2023-12-19_alcohol_sex_test_ldsc_2023-12-19"
identifier_parameter="tcw_2023-12-19_alcohol_sex_test"

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"

path_directory_source_rescue="${path_directory_dock}/gwas_preparation_${identifier_preparation}/5_fill_dbsnp_rs_identifiers"
path_directory_source="${path_directory_dock}/gwas_preparation_${identifier_preparation}/6_gwas_clean_gwas2vcf"
path_directory_product="${path_directory_dock}/gwas_preparation_${identifier_preparation}/7_filter_constrain_gwas_values"
path_directory_parent_temporary="${path_directory_process}/temporary_check_9"

# Files.
# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_script_check="${path_directory_partner_scripts}/gwas_clean/check_gwas_summary_values.sh"

# Initialize directories.
#rm -r $path_directory_parent_temporary # caution
#mkdir -p $path_directory_parent_temporary
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.



##########
# Check sets of GWAS summary statistics to describe potential problems.

# $ ls ./*.txt.gz | wc -l
# $ du -hs ./*.txt.gz

# records that do not pass checks: 1 of 21,795,381 (header line)
#/usr/bin/bash $path_script_check \
#"${path_directory_source}/32769997_zhou_2020_thyroid_hormone.txt.gz" \
#$path_directory_parent_temporary \
#$report

##########
# Copy the GWAS summary statistics from the previous process.
# Most sets of GWAS summary statistics do not need extra processing.
# Subsequent processes on a few studies will replace the appropriate files.
cp $path_directory_source/*.txt.gz $path_directory_product



##########
# Remove files of GWAS summary statistics that failed GWAS2VCF procedure.

#rm $path_directory_product/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz
#rm $path_directory_product/32242144_revez_2020_vitamin_d.txt.gz
#rm $path_directory_product/32099098_polimanti_2020_eur_opioid_dep_exposed.txt.gz
#rm $path_directory_product/32099098_polimanti_2020_eur_opioid_dep_unexposed.txt.gz
#rm $path_directory_product/32099098_polimanti_2020_eur_opioid_exposure.txt.gz
#rm $path_directory_product/30718901_howard_2019_pgc_ukb.txt.gz

##########
# Copy GWAS summary statistics from before GWAS2VCF where appropriate.
# Some studies will be useful for analysis in LDSC even if they lack sufficient
# information to pass filters in GWAS2VCF.

#cp $path_directory_source_rescue/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz $path_directory_product
#cp $path_directory_source_rescue/32242144_revez_2020_vitamin_d.txt.gz $path_directory_product
#cp $path_directory_source_rescue/32099098_polimanti_2020_eur_opioid_dep_exposed.txt.gz $path_directory_product
#cp $path_directory_source_rescue/32099098_polimanti_2020_eur_opioid_dep_unexposed.txt.gz $path_directory_product
#cp $path_directory_source_rescue/32099098_polimanti_2020_eur_opioid_exposure.txt.gz $path_directory_product
#cp $path_directory_source_rescue/30718901_howard_2019_pgc_ukb.txt.gz $path_directory_product

# Remove temporary, intermediate files.
#rm -r $path_directory_parent_temporary



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "7_check_filter_constrain_gwas_values.sh"
  echo "----------"
fi



#
