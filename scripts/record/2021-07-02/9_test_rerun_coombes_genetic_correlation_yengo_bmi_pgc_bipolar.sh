
################################################################################
# This script exists to test method for genetic correlation analysis by comparison to analysis by Brandon J. Coombes, Ph.D.
################################################################################

#!/bin/bash

################################################################################
# Organize paths.

# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_ldsc=$(<"./tools_ldsc.txt")
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_gwas_summaries_team=$(<"./gwas_summaries_team.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")

path_gwas_30124842_yengo_2018="${path_gwas_summaries}/30124842_yengo_2018/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz"
path_gwas_30124842_yengo_2018_coombes="${path_gwas_summaries_team}/REFORMATTED/BMI_GIANTUKBB.txt.gz"
path_gwas_coombes=$(<"./gwas_pgc_bipolar_bmi_coombes.txt")
path_gwas_pgc_bipolar_bmi_coombes="${path_gwas_coombes}/ALL_BMI.pcAdj.assoc.linear_P_MAresultsFE.FUMA.txt.gz"

path_dock="${path_process}/dock"
path_genetic_reference="${path_dock}/access/genetic_reference_test_coombes"
path_alleles="${path_genetic_reference}/alleles"
path_disequilibrium="${path_genetic_reference}/disequilibrium"

path_genetic_correlation="${path_dock}/genetic_correlation_test_coombes"

################################################################################
# Organize directories.

rm -r $path_genetic_reference
rm -r $path_genetic_correlation
# Determine whether the temporary directory structure already exists.
if [ ! -d $path_genetic_reference ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_genetic_reference
    mkdir -p $path_alleles
    mkdir -p $path_disequilibrium
    mkdir -p $path_genetic_correlation
fi

################################################################################
# Access references for LDSC.

cd $path_genetic_reference

# Definitions of Simple Nucleotide Variant alleles.
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 "$path_genetic_reference/w_hm3.snplist.bz2"
mv "$path_genetic_reference/w_hm3.snplist" "$path_alleles/w_hm3.snplist"
# w_hm3.snplist

# Linkage disequilibrium scores for European population.
# For simple heritability estimation.
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -xjvf eur_w_ld_chr.tar.bz2 -C $path_disequilibrium
# dock/access/disequilibrium/eur_w_ld_chr/*

################################################################################
# Format GWAS summary statistics for analysis in LDSC.

echo "test test test"

cd $path_genetic_correlation

# Waller accession of BMI GWAS summary statistics.
echo "SNP A1 A2 N BETA P" > gwas_30124842_yengo_2018_format
zcat $path_gwas_30124842_yengo_2018 | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $3, toupper($4), toupper($5), $10, $7, $9}' >> gwas_30124842_yengo_2018_format
head gwas_30124842_yengo_2018_format

# Coombes accession of BMI GWAS summary statistics.
echo "SNP A1 A2 N BETA P" > gwas_30124842_yengo_2018_coombes_format
zcat $path_gwas_30124842_yengo_2018_coombes | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $1, toupper($2), toupper($3), $10, $5, $4}' >> gwas_30124842_yengo_2018_coombes_format
head gwas_30124842_yengo_2018_coombes_format

# Coombes BMI in Bipolar Disorder (PGC cohort) GWAS summary statistics.
echo "SNP A1 A2 N BETA P" > gwas_pgc_bipolar_bmi_format
zcat $path_gwas_pgc_bipolar_bmi_coombes | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $1, toupper($4), toupper($5), 4332, $7, $9}' >> gwas_pgc_bipolar_bmi_format
head gwas_pgc_bipolar_bmi_format

################################################################################
# Munge GWAS summary statistics for analysis in LDSC.

cd $path_genetic_correlation

# Waller accession of BMI GWAS summary statistics.
$path_ldsc/munge_sumstats.py \
--sumstats gwas_30124842_yengo_2018_format \
--out gwas_bmi_waller \
--merge-alleles $path_alleles/w_hm3.snplist \

# Coombes accession of BMI GWAS summary statistics.
$path_ldsc/munge_sumstats.py \
--sumstats gwas_30124842_yengo_2018_coombes_format \
--out gwas_bmi_coombes \
--merge-alleles $path_alleles/w_hm3.snplist \

# Coombes BMI in Bipolar Disorder (PGC cohort) GWAS summary statistics.
$path_ldsc/munge_sumstats.py \
--sumstats gwas_pgc_bipolar_bmi_format \
--signed-sumstats BETA,0 \
--out gwas_bipolar_bmi_coombes \
--merge-alleles $path_alleles/w_hm3.snplist \

if false; then
  ################################################################################
  # Genetic correlation in LDSC.

  $path_ldsc/ldsc.py \
  --rg gwas_bmi_waller.sumstats.gz,gwas_pgc_bipolar_bmi_format.sumstats.gz \
  --ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
  --w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
  --out bmi_control_waller_versus_bmi_bipolar_coombes

  $path_ldsc/ldsc.py \
  --rg gwas_bmi_coombes.sumstats.gz,gwas_pgc_bipolar_bmi_format.sumstats.gz \
  --ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
  --w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
  --out bmi_control_coombes_versus_bmi_bipolar_coombes

fi
