#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 21 September 2023
# Date, last execution: ___ 2023
# Date, review: ___ 2023
################################################################################
# Note

# Perform this procedure judiciously.
# The purpose of this procedure is to run a few final checks to prepare GWAS
# summary statistics for analysis in LDSC, LDpred2, or SBayesR.
# 1. Remove any records (rows) with empty values (cells).
# 2. Constrain probabilities (p-values) to double-float precision
# (1E-308 to __).
# 3. Remove records for any SNPs with missing counts of observations, as these
# SNPs raise an error in LDpred2.


################################################################################
# Organize paths.


################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.


################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "6_check_constrain_gwas_values.sh"
  echo "----------"
fi



#
