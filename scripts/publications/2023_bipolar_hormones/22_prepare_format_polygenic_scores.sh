#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 21 March 2023
# Date, last execution: 19 April 2023
# Review: TCW; 19 April 2023
################################################################################
# Note

# Source Format
# Description: Format polygenic scores from SBayesR after combination across chromosomes.
# File suffix: ".txt"
# File type: text
# File compression: none
# Delimiter: Tab
# Columns: identifier count_allele_total count_allele_dosage score_sum score_mean (TCW; 2023-03-21)
#          1          2                  3                   4         5

# Source Format
# Description: Format polygenic scores from LDpred2.
# File suffix: ".txt"
# File type: text
# File compression: none
# Delimiter: White space
# Columns: FID   IID   [score name]_auto (TCW; 2023-03-21)
#          1     2     3

# Source Format
# Description: Format polygenic scores from PRSice2.
# File suffix: ".txt"
# File type: text
# File compression: none
# Delimiter: White space
# Columns: FID   IID   In_Regression   PRS (TCW; 2023-03-21)
#          1     2     3               4

# Product Format
# Description: Format polygenic scores standard for collection and standardization.
# File suffix: ".txt"
# File type: text
# File compression: none
# Delimiter: Tab
# Columns: identifier   score (TCW; 2023-03-21)
#          1            2

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes

path_directory_source="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/sbayesr_plink_polygenic_scores_combination"
path_directory_product="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/sbayesr_plink_polygenic_scores_preparation"

# Files.
path_file_source_1="${path_directory_source}/score_32769997_zhou_2020_tsh.txt"
path_file_product_1="${path_directory_product}/score_32769997_zhou_2020_tsh.txt"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

###########################################################################
# Organize parameters.



###########################################################################
# Execute procedure.

# Echo each command to console.
set -x

##########
# Translate polygenic scores to a standard format.
# AWK interprets a single space delimiter (FS=" ") as any white space.

# SBayesR

echo -e "identifier\tscore" > $path_file_product_1
cat $path_file_source_1 | awk 'BEGIN { FS="\t"; OFS="\t" } NR > 1 {
  print $1, $5
}' >> $path_file_product_1

# LDpred2

# PRSice2
# Use column "Pt_0.05" for score calculated on probability threshold 0.05.



#
