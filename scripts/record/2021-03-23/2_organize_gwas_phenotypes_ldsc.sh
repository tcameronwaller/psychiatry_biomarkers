#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "Concatenate GWAS summary statistics across chromosomes."
echo "Psychiatric Genetics Consortium GWAS data on alcoholism."
echo "Human genome version: GRCh37, hg19"
echo "version 1"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo ""
echo ""
echo ""

# Organize paths.
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths

# Phenotype GWAS summary statistics.
path_gwas_team_directory=$(<"./gwas_summaries_team.txt")


path_gwas_alcohol_raw="$path_gwas_team_directory/AUD/pgc_alcdep.discovery.aug2018_release.txt.gz"
path_gwas_bipolar_1_raw="$path_gwas_team_directory/BD/daner_bip_pgc3_nm.gz"
path_gwas_bipolar_2_raw="$path_gwas_team_directory/BD/BD2/daner_bip_pgc3_nm.gz"

path_scripts="$path_waller/sexy_alcohol/scripts/record/2021-02-25"
path_calculate_z_score="$path_scripts/calculate_z_score.sh"
path_temporary=$(<"./processing_sexy_alcohol.txt")
path_waller="$path_temporary/waller"
path_dock="$path_temporary/waller/dock"
path_gwas="$path_dock/gwas"
path_gwas_alcohol_format_directory="$path_temporary/waller/dock/gwas/female_male_alcoholism_pgc"
path_gwas_alcohol_format="$path_gwas_alcohol_format_directory/gwas_format_pgc_alcoholism.txt"

# Initialize directories.
#rm -r $path_gwas
if [ ! -d $path_gwas_alcohol_format_directory ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_gwas_alcohol_format_directory
fi

# Organize information from linear GWAS.

echo "SNP A1 A2 N BETA P" > $path_gwas_alcohol_format
# Format of GWAS summary for LDSC.
# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics
# description: ............................ PGC column ... LDSC column
# variant identifier: ....................... "SNP" ........ "SNP"
# alternate allele (effect allele): ......... "A1" ......... "A1"
# reference allele (non-effect allele): ..... "A2" ......... "A2"
# sample size: .............................. 52,848 ....... "N"
# effect (beta): ............................ "Z" .......... "BETA"
# probability (p-value): .................... "P" .......... "P"

# SNP: split($2,a,":"); print a[1]
# A1: toupper($4)
# A2: toupper($5)
# N: (cases: 14,904; controls: 37,944; total: 52,848)
# BETA: $6
# P: $7
zcat $path_gwas_alcohol_raw | awk 'NR > 1 {split($2,a,":"); print a[1], toupper($4), toupper($5), (52848), $6, $7}' >> $path_gwas_alcohol_format
# Calculate Z-score standardization of Beta coefficients.
#/usr/bin/bash $path_calculate_z_score 5 $path_gwas_alcohol_format $path_gwas_alcohol_format
gzip $path_gwas_alcohol_format
echo "after format..."
head -30 "$path_gwas_alcohol_format.gz"

echo "----------"
echo "----------"
echo "----------"
