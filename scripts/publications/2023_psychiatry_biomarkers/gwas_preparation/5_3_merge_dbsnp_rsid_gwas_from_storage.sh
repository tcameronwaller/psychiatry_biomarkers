#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 30 December 2023
# Date, last execution: 30 December 2023
# Date, review: 30 December 2023
################################################################################
# Note

# This process merges in the files from the previous batch that included the
# "fill_dbsnp_rsid.sh" script process. This batch ran on 22 December 2023.
# Count of dbSNP rsID files: 40

################################################################################
# Organize paths.

# Identifiers or designators of parameter version and preparation batch.
identifier_preparation="gwas_preparation_2023-12-30"
identifier_parameter="tcw_2023-12-30"
identifier_preparation_dbsnp_rsid="gwas_2023-12-22_dbsnp_rsid"
identifier_parameter_dbsnp_rsid="tcw_2023-12-30_dbsnp_rsid"
# This parameter table includes copies of studies with suffix "_dbsnp_rsid" for
# those studies copied and merged from the previous batch procedure that
# included the "fill_dbsnp_rsid.sh" script process.

# Directories.
cd ~/paths
path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"

path_directory_source_dbsnp_rsid="${path_directory_gwas_summaries}/organization/${identifier_preparation_dbsnp_rsid}/5_fill_dbsnp_rs_identifiers"
path_directory_source="${path_directory_dock}/${identifier_preparation}/4_filter_constrain_gwas_values"
path_directory_product="${path_directory_dock}/${identifier_preparation}/5_fill_dbsnp_rs_identifiers"
path_directory_product_temporary="${path_directory_product}/temporary_dbsnp_rsid_raw"

# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
mkdir -p $path_directory_product_temporary
cd $path_directory_product

################################################################################
# Organize parameters.

report="true"
set -x

################################################################################
# Execute procedure.



##########
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "Count of files in current main batch:"
  ls $path_directory_source | wc -l
  echo "----------"
  echo "Count of files in previous batch including fill dbSNP rsIDs:"
  ls $path_directory_source_dbsnp_rsid | wc -l
  echo "----------"
fi

##########
# Copy the GWAS summary statistics from the previous process.
# Most sets of GWAS summary statistics do not need extra processing.
# Subsequent processes on a few studies will replace the appropriate files.
cp $path_directory_source/*.txt.gz $path_directory_product

##########
# Copy the GWAS summary statistics from the previous batch procedure that
# included the "fill_dbsnp_rsid.sh" script process.
cp $path_directory_source_dbsnp_rsid/*.txt.gz $path_directory_product_temporary

##########
# Append suffix to the names of all files from the previous batch procedure that
# included the "fill_dbsnp_rsid.sh" script process.
cd $path_directory_product_temporary
for name_file in *.txt.gz
do
  name_base_file="$(basename $name_file .txt.gz)"
  mv "$name_file" "${name_base_file}_dbsnp_rsid.txt.gz"
done

##########
# Merge the files from the dbSNP rsID files with those from the main batch.
for path_file in "$path_directory_product_temporary/*.txt.gz"
do
  echo $path_file
  cp $path_file $path_directory_product
done

##########
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "Count of files in product directory after merge:"
  ls $path_directory_product | wc -l
  echo "----------"
fi

##########
# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
