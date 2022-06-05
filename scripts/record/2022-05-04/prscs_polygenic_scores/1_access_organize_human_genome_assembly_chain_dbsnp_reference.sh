#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_bcftools=$(<"./tools_bcftools.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="${path_process}/dock"
path_parameters="${path_dock}/parameters"
path_translations_chromosomes_dbsnp_grch37="${path_parameters}/promiscuity/translations_chromosomes_refseq_grch37p13.txt"
#path_translations_chromosomes_dbsnp_grch38="${path_parameters}/promiscuity/translations_chromosomes_refseq_grch38p14.txt"
path_human_genome_sequence_parent_container="${path_dock}/access/human_genome_sequence"
path_human_genome_assembly_chain_parent_container="${path_dock}/access/human_genome_assembly_chain"
path_dbsnp_parent_container="${path_dock}/access/dbsnp"
path_dbsnp_grch37_raw_container="${path_dock}/access/dbsnp/grch37_raw"
path_dbsnp_grch37_format_container="${path_dock}/access/dbsnp/grch37_format"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"

path_script_access_human_genome="${path_promiscuity_scripts}/utility/access_genetic_reference/access_genome_sequence_human_grch37_grch38_fasta.sh"
path_script_access_genome_assembly_chain="${path_promiscuity_scripts}/utility/access_genetic_reference/access_genome_human_assembly_map_chain_ucsc_ensembl.sh"

path_script_access_dbsnp="${path_promiscuity_scripts}/utility/access_genetic_reference/access_dbsnp_human_grch37_grch38_vcf.sh"
path_script_chromosome_in_vcf="${path_promiscuity_scripts}/utility/bcftools/translate_chromosomes_in_vcf.sh"

###########################################################################
# Execute procedure.
###########################################################################

################################################################################
# Access and organize information for human genome sequence.

# Access.
if false; then
  /usr/bin/bash "${path_script_access_human_genome}" \
  "grch37" \
  $path_human_genome_sequence_parent_container
fi

if true; then
  /usr/bin/bash "${path_script_access_human_genome}" \
  "grch38" \
  $path_human_genome_sequence_parent_container
fi


################################################################################
# Access and chain files to map between assemblies of human genome.

# Access.
if false; then
  /usr/bin/bash "${path_script_access_genome_assembly_chain}" \
  $path_human_genome_assembly_chain_parent_container
fi

################################################################################
# Access and organize information for annotations on SNPs.

# Access.
if false; then
  /usr/bin/bash "${path_script_access_dbsnp}" \
  "grch37" \
  $path_dbsnp_parent_container
fi

# Format genotype information in VCF.
if false; then
  # Initialize directory.
  rm -r $path_dbsnp_grch37_format_container
  mkdir -p $path_dbsnp_grch37_format_container
  # Organize specific paths and parameters.
  path_vcf_source="${path_dbsnp_grch37_raw_container}/GCF_000001405.25.gz"
  path_vcf_product="${path_dbsnp_grch37_format_container}/GCF_000001405.25.gz" # determine suffix from BCFTools argument
  threads=16
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_chromosome_in_vcf}" \
  $path_translations_chromosomes_dbsnp_grch37 \
  $path_vcf_source \
  $path_vcf_product \
  $threads \
  $path_bcftools \
  $report
fi
