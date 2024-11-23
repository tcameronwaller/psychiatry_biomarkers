"""
Supply functionality for process and analysis of data from genetic correlation
analyses between sets of GWAS summary statistics.

This module 'thyroid_organization' is part of the 'genetic_correlation' package
within the 'psychiatry_biomarkers' package.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Rochester, Minnesota 55902
    United States of America

License:

    This file is part of the project package directory 'psychiatry_biomarkers'
    (https://github.com/tcameronwaller/psychiatry_biomarkers/).

    Project 'psychiatry_biomarkers' supports data analysis with team in
    psychiatry and pharmacogenomics.
    Copyright (C) 2024 Thomas Cameron Waller

    The code within project 'psychiatry_biomarkers' is free software: you can
    redistribute it and/or modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation, either version 3 of
    the GNU General Public License, or (at your option) any later version.

    The code within project 'psychiatry_biomarkers' is distributed in the hope
    that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with project 'psychiatry_biomarkers'. If not, see
    <http://www.gnu.org/licenses/>.
"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import sys
#print(sys.path)
import os
import math
import statistics
import pickle
import copy
import random
import itertools
import time

# Relevant

import numpy
import scipy.stats
import pandas
pandas.options.mode.chained_assignment = None # default = "warn"

# Custom
import partner.utility as putly
import partner.extraction as pextr
import partner.organization as porg
import partner.description as pdesc
#import partner.regression as preg
import partner.plot as pplot
import partner.parallelization as prall

###############################################################################
# Functionality



##########
# 1. Initialize directories for read of source and write of product files.


def initialize_directories(
    project=None,
    routine=None,
    procedure=None,
    path_directory_dock=None,
    restore=None,
    report=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        project (str): name of project
        routine (str): name of routine, either 'transcriptomics' or
            'proteomics'
        procedure (str): name of procedure, a set or step in the routine
            process
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        restore (bool): whether to remove previous versions of data
        report (bool): whether to print reports

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    paths["dock"] = path_directory_dock
    paths["in_data"] = os.path.join(
        paths["dock"], "in_data", str(project), str(routine),
    )
    paths["in_parameters"] = os.path.join(
        paths["dock"], "in_parameters", str(project), str(routine),
    )
    paths["in_parameters_private"] = os.path.join(
        paths["dock"], "in_parameters_private", str(project), str(routine),
    )
    paths["out_project"] = os.path.join(
        paths["dock"], str("out_" + project),
    )
    paths["out_routine"] = os.path.join(
        paths["out_project"], str(routine),
    )
    paths["out_procedure"] = os.path.join(
        paths["out_routine"], str(procedure),
    )
    paths["out_test"] = os.path.join(
        paths["out_procedure"], "test",
    )
    paths["out_data"] = os.path.join(
        paths["out_procedure"], "data",
    )
    paths["out_plot"] = os.path.join(
        paths["out_procedure"], "plot",
    )
    # Initialize directories in main branch.
    paths_initialization = [
        #paths["out_project"],
        #paths["out_routine"],
        paths["out_procedure"], # omit to avoid conflict in parallel branches
        paths["out_test"],
        paths["out_data"],
        paths["out_plot"],
    ]
    # Remove previous directories and files to avoid version or batch
    # confusion.
    if restore:
        for path in paths_initialization:
            putly.remove_directory(path=path) # caution
            pass
    # Create directories.
    for path in paths_initialization:
        putly.create_directories(
            path=path,
        )
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print(
            "module: psychiatry_biomarkers.genetic_correlation." +
            "thyroid_organization.py"
        )
        print("function: initialize_directories()")
        putly.print_terminal_partition(level=5)
        print("path to dock directory for procedure's files: ")
        print(path_directory_dock)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return paths


def initialize_directory_group_analysis(
    group_analysis=None,
    paths=None,
    restore=None,
):
    """
    Initialize child directory for analysis group.

    arguments:
        group_analysis (str): name for analysis group of comparisons
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        restore (bool): whether to remove previous versions of data

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Copy information.
    paths = copy.deepcopy(paths)
    # Define paths to directories.
    paths["out_data_group_analysis"] = os.path.join(
        paths["out_data"], group_analysis,
    )
    paths["out_plot_group_analysis"] = os.path.join(
        paths["out_plot"], group_analysis,
    )
    # Remove previous files to avoid version or batch confusion.
    if restore:
        putly.remove_directory(path=paths["out_data_group_analysis"])
        putly.remove_directory(path=paths["out_plot_group_analysis"])
    # Initialize directories.
    putly.create_directories(
        path=paths["out_data_group_analysis"]
    )
    putly.create_directories(
        path=paths["out_plot_group_analysis"]
    )
    # Return information.
    return paths


##########
# 2.1. Read parameters about studies


def read_source_parameter_studies_process(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Primary studies: 80 studies; Disorders of Neurology, Psychiatry, and
    # Substance Use.
    # Secondary studies: 164 studies; Disorders of Thyroid Physiology,
    # Biomarkers of Thyroid Physiology, Biomarkers of Sex Hormone Physiology,
    # Biomarkers of other physiology and metabolism.

    # Define path to parent directory.
    path_directory_parent = os.path.join(
        paths["in_parameters_private"],
    )
    # Define paths to files.
    path_file_table_studies = os.path.join(
        path_directory_parent,
        "table_gwas_translation_tcw_2023-12-30.tsv",
    )
    # Read information from file.
    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["inclusion"] = "int32"
    types_columns["study"] = "string"
    types_columns["phenotype"] = "string"
    types_columns["sex"] = "string"
    types_columns["observations_total"] = "float32"
    types_columns["cases"] = "float32"
    types_columns["controls"] = "float32"
    types_columns["observations_effective"] = "float32"
    types_columns["prevalence_sample"] = "float32"
    types_columns["prevalence_population"] = "float32"
    table = pandas.read_csv(
        path_file_table_studies,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )
    # Return information.
    return table


def read_source_parameter_studies_polish(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Primary studies: 80 studies; Disorders of Neurology, Psychiatry, and
    # Substance Use.
    # Secondary studies: 164 studies; Disorders of Thyroid Physiology,
    # Biomarkers of Thyroid Physiology, Biomarkers of Sex Hormone Physiology,
    # Biomarkers of other physiology and metabolism.

    # Define path to parent directory.
    path_directory_parent = os.path.join(
        paths["in_parameters_private"],
    )
    # Define paths to files.
    path_file_table_studies = os.path.join(
        path_directory_parent,
        "table_studies_attributes.tsv",
    )
    # Read information from file.
    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["inclusion_thyroid_table"] = "float"
    types_columns["inclusion_thyroid_figure"] = "float"
    types_columns["inclusion_thyroid_figure_1"] = "float"
    types_columns["inclusion_thyroid_figure_2"] = "float"
    types_columns["inclusion_sex_table"] = "float"
    types_columns["inclusion_sex_alcohol_tobacco_table"] = "float"
    types_columns["sort"] = "float"
    types_columns["group"] = "string"
    types_columns["sort_group"] = "float"
    types_columns["identifier"] = "string"
    types_columns["abbreviation"] = "string"
    types_columns["description"] = "string"
    types_columns["sex"] = "string"
    types_columns["ancestry"] = "string"
    types_columns["author"] = "string"
    types_columns["year"] = "string"
    types_columns["pubmed"] = "string"
    table = pandas.read_csv(
        path_file_table_studies,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )
    # Return information.
    return table


##########
# 2.2. Read and organize SNP heritabilities in table for supplement


def read_source_data_snp_heritability(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Primary studies: 80 studies; Disorders of Neurology, Psychiatry, and
    # Substance Use.
    # Secondary studies: 164 studies; Disorders of Thyroid Physiology,
    # Biomarkers of Thyroid Physiology, Biomarkers of Sex Hormone Physiology,
    # Biomarkers of other physiology and metabolism.

    # Define path to parent directory.
    path_directory_parent = os.path.join(
        paths["in_data"],
        "gwas_2023-12-30_ldsc_2024-01-08_extraction_2024-05-22",
        "extraction_2024-05-22",
    )
    # Define path to child directories.
    path_directory_heritability = os.path.join(
        path_directory_parent,
        "5_gwas_heritability_ldsc",
    )
    # Define paths to files.
    path_file_table_heritability = os.path.join(
        path_directory_heritability,
        "table_heritability.tsv",
    )
    # Collect information.
    pail = dict()
    # Read and organize information from file.
    types_columns = pextr.define_snp_heritability_table_column_types()
    pail["table_heritability"] = pandas.read_csv(
        path_file_table_heritability,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Return information.
    return pail


def read_organize_source_parameter_snp_heritability(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Primary studies: 80 studies; Disorders of Neurology, Psychiatry, and
    # Substance Use.
    # Secondary studies: 164 studies; Disorders of Thyroid Physiology,
    # Biomarkers of Thyroid Physiology, Biomarkers of Sex Hormone Physiology,
    # Biomarkers of other physiology and metabolism.

    # Read information from file.
    table_process = read_source_parameter_studies_process(
        paths=paths,
        report=report,
    )
    table_polish = read_source_parameter_studies_polish(
        paths=paths,
        report=report,
    )

    # Translate names of columns.
    translations_column = dict()
    translations_column["study"] = "identifier"
    table_process.rename(
        columns=translations_column,
        inplace=True,
    )

    # Specify sequence of columns within table.
    columns_process = [
        "inclusion",
        "identifier",
        "observations_total",
        "cases",
        "controls",
        "observations_effective",
        "prevalence_sample",
        "prevalence_population",
    ]
    # Filter and sort columns within table.
    table_process = porg.filter_sort_table_columns(
        table=table_process,
        columns_sequence=columns_process,
        report=report,
    )
    # Specify sequence of columns within table.
    columns_polish = [
        "inclusion_thyroid_table",
        "inclusion_sex_table",
        "inclusion_sex_alcohol_tobacco_table",
        "inclusion_gonad",
        "sort",
        "group",
        "sort_group",
        "identifier",
        "abbreviation",
        "description",
        "sex",
        "ancestry",
        "author",
        "year",
        "pubmed",
    ]
    # Filter and sort columns within table.
    table_polish = porg.filter_sort_table_columns(
        table=table_polish,
        columns_sequence=columns_polish,
        report=report,
    )

    # Merge tables with information about studies.
    table_study = porg.merge_columns_two_tables(
        identifier_first="identifier",
        identifier_second="identifier",
        table_first=table_polish,
        table_second=table_process,
        preserve_index=False,
        report=report,
    )

    # Specify sequence of columns within table.
    columns_sequence = [
        "inclusion_thyroid_table",
        "inclusion_sex_table",
        "inclusion_sex_alcohol_tobacco_table",
        "inclusion_gonad",
        "sort",
        "group",
        "sort_group",
        "identifier",
        "abbreviation",
        "description",
        "sex",
        "ancestry",
        "observations_total",
        "cases",
        "controls",
        "observations_effective",
        "prevalence_sample",
        "prevalence_population",
        "author",
        "year",
        "pubmed",
    ]
    # Filter and sort columns within table.
    table_study = porg.filter_sort_table_columns(
        table=table_study,
        columns_sequence=columns_sequence,
        report=report,
    )
    # Collect information.
    pail = dict()
    pail["table_study"] = table_study
    # Return information.
    return pail


# TODO: TCW; 21 October 2024
# Temporarily necessary to switch within this function between variables for the
# "inclusion" filter.
# Currently set for the sex-hormones-alcohol-tobacco project... <-- not anymore
# Currently set for the sex-hormones versus psychiatric and substance use disorders


def control_read_organize_snp_heritability_table_supplement(
    paths=None,
    report=None,
):
    """
    Control procedure to read from file the SNP heritabilities from LDSC and
    organize these within a table for the article's supplement.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Read source information from file.
    source_heritability = read_source_data_snp_heritability(
        paths=paths,
        report=report,
    )
    table_heritability_raw = source_heritability["table_heritability"]
    # Organize identifier of study.
    table_heritability_raw["identifier"] = table_heritability_raw.apply(
        lambda row: str(row["name_file"]).strip().replace(".log", ""),
        axis="columns", # apply function to each row
    )
    # Specify sequence of columns within table.
    columns_heritability = (
        pextr.define_snp_heritability_table_column_sequence()
    )
    # Filter and sort columns within table.
    table_process = porg.filter_sort_table_columns(
        table=table_heritability_raw,
        columns_sequence=columns_heritability,
        report=report,
    )

    # Read source information from file.
    source_study = read_organize_source_parameter_snp_heritability(
        paths=paths,
        report=report,
    )
    # Merge tables with information about studies.
    table_heritability = porg.merge_columns_two_tables(
        identifier_first="identifier",
        identifier_second="identifier",
        table_first=source_study["table_study"],
        table_second=table_heritability_raw,
        preserve_index=False,
        report=report,
    )
    # Filter rows in table.
    inclusion = "inclusion_thyroid_table"
    #inclusion = "inclusion_sex_alcohol_tobacco_table"
    #inclusion = "inclusion_gonad"
    table_heritability = table_heritability.loc[
        (
            (table_heritability[inclusion] == 1)
        ), :
    ].copy(deep=True)
    # Sort rows within table.
    table_heritability.sort_values(
        by=[
            "sort",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )

    # Specify sequence of columns within table.
    columns_sequence = [
        #"inclusion_thyroid_table",
        #"path_directory",
        #"name_file",
        #"type_analysis",
        #"sort",
        "group",
        #"sort_group",
        "identifier",
        "abbreviation",
        "description",
        "sex",
        "ancestry",
        "observations_total",
        "cases",
        "controls",
        "observations_effective",
        "prevalence_sample",
        "prevalence_population",
        "author",
        "year",
        "pubmed",
        "variants",
        "heritability",
        "heritability_error",
        "heritability_ci95_low",
        "heritability_ci95_high",
        "heritability_ci99_low",
        "heritability_ci99_high",
        "lambda_gc",
        "chi_square",
        "intercept",
        "intercept_error",
        "ratio",
        "ratio_error",
        #"summary_heritability_error",
        #"summary_heritability_ci95",
        #"summary_heritability_ci99",
    ]
    # Filter and sort columns within table.
    table_heritability = porg.filter_sort_table_columns(
        table=table_heritability,
        columns_sequence=columns_sequence,
        report=report,
    )

    ##########
    # Collect information.
    # Collections of files.
    pail_write_tables = dict()
    pail_write_tables[str("table_heritability")] = table_heritability
    ##########
    # Write product information to file.
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=paths["out_data"],
        reset_index_rows=False,
        write_index_rows=False,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Source table before organization:")
        print(table_heritability)
        print("Column labels:")
        labels_columns = table_heritability.columns.to_list()
        print(labels_columns)
        putly.print_terminal_partition(level=4)
        count_columns = (table_heritability.shape[1])
        print("count of columns in table: " + str(count_columns))
        count_rows = (table_heritability.shape[0])
        print("count of rows in source table: " + str(count_rows))
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return table_heritability


##########
# 2.3. Read and assemble genetic correlations


def read_organize_source_data_genetic_correlations(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables with entry names
            (keys) derived from original names of files

    """

    # Primary studies: 80 studies; Disorders of Neurology, Psychiatry, and
    # Substance Use.
    # Secondary studies: 164 studies; Disorders of Thyroid Physiology,
    # Biomarkers of Thyroid Physiology, Biomarkers of Sex Hormone Physiology,
    # Biomarkers of other physiology and metabolism.

    # Define path to parent directory.
    path_directory_parent = os.path.join(
        paths["in_data"],
        "gwas_2023-12-30_ldsc_2024-01-08_extraction_2024-05-22",
        "extraction_2024-05-22",
    )
    # Define path to child directories.
    path_directory_one_one = os.path.join(
        path_directory_parent,
        "6_gwas_correlation_ldsc_primary",
    )
    #path_directory_two_two = os.path.join(
    #    path_directory_parent,
    #    "6_gwas_correlation_ldsc_secondary_thyroid",
    #)
    path_directory_two_two = os.path.join(
        path_directory_parent,
        "6_gwas_correlation_ldsc_secondary",
    )
    path_directory_one_two = os.path.join(
        path_directory_parent,
        "6_gwas_correlation_ldsc_primary_secondary",
    )
    # Define paths to files.
    path_file_table_one_one = os.path.join(
        path_directory_one_one, "table_neuropsychiatry_substance_disorders.tsv",
    )
    path_file_table_two_two = os.path.join(
        path_directory_two_two, "table_physiology_biomarkers.tsv",
    )
    # Collect information.
    pail = dict()
    # Read and organize information from file.
    # Genetic correlations between primary-primary studies, secondary-secondary
    # studies, and primary-secondary studies.
    pail["table_rg_one_one"] = (
        pextr.read_organize_table_ldsc_correlation_single(
            path_file_table=path_file_table_one_one,
            report=report,
    ))
    pail["table_rg_two_two"] = (
        pextr.read_organize_table_ldsc_correlation_single(
            path_file_table=path_file_table_two_two,
            report=report,
    ))
    # Correlations between primary-secondary studies.
    pail["table_rg_one_two"] = (
        pextr.read_organize_table_ldsc_correlation_multiple(
            path_directory_parent=path_directory_one_two,
            report=report,
    ))
    # Return information.
    return pail


def control_assemble_genetic_correlations(
    paths=None,
    report=None,
):
    """
    Control procedure to assemble within a single table the information about
    genetic correlations from LDSC.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Read source information from file.
    source = read_organize_source_data_genetic_correlations(
        paths=paths,
        report=report,
    )

    # Organize source information within a single table.
    # Concatenate tables.
    table = pandas.concat(
        [
            source["table_rg_one_one"],
            source["table_rg_two_two"],
            source["table_rg_one_two"],
        ],
        axis="index",
        join="outer",
        ignore_index=True,
        copy=True,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Source table before organization:")
        print(table)
        print("Column labels:")
        labels_columns = table.columns.to_list()
        print(labels_columns)
        putly.print_terminal_partition(level=4)
        count_columns = (table.shape[1])
        print("Count of columns in table: " + str(count_columns))

        count_rows = (table.shape[0])
        print("Count of rows in source table: " + str(count_rows))
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return table


##########
# 2.4. Read and organize information to plot charts


def organize_source_plot_table_index(
    table=None,
    report=None,
):
    """
    Organizes index of table from source.

    arguments:
        table (object): Pandas data-frame table of genetic correlations
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    ##########
    # Copy information in table.
    table_copy = table.copy(deep=True)

    ##########
    # Organize information in table.
    table_copy.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_copy.set_index(
        [
            "group_analysis",
            "abbreviation_primary",
        ],
        append=False,
        drop=True,
        inplace=True,
    )

    ##########
    # Return information.
    return table_copy


def read_organize_source_plot(
    group_analysis=None,
    name_table=None,
    path_directory_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        group_analysis (str): name of child directory for analysis group of
            comparisons
        name_table (str): base name for files to which to save tables
        path_directory_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Genetic correlations:
    # Thyroid Disorders and Thyroid Biomarkers
    # against Psychiatric and Substance Use Disorders
    ##########

    # Define path to parent directory.
    path_directory_parent = os.path.join(
        path_directory_dock, "out_psychiatry_biomarkers",
        "genetic_correlation", "thyroid_organization", "data"
    )
    # Define path to child directory.
    path_directory_child = os.path.join(
        path_directory_parent, group_analysis, "pickle",
    )
    # Define path to file.
    path_file_table = os.path.join(
        path_directory_child,
        str(name_table + "_for_plot.pickle"),
    )
    # Read information from file.
    table_raw = pandas.read_pickle(
        path_file_table,
    )
    print("here is the table before organizing the index...")
    print(table_raw)
    # Organize information in table from source.
    # Collect information.
    table = organize_source_plot_table_index(
        table=table_raw,
        report=report,
    )
    print("here is the table after organizing the index...")
    print(table)
    # This functionality is useful for reading a similarly-formatted table from
    # tab-delimited text file.
    if False:
        # Read and organize information from file.
        pail = putly.read_table_multiindex_columns_transpose(
            path_file_table=path_file_table,
            row_index_place_holder="group_secondary",
            name_row_index="group_secondary",
            name_column_index_1="group_primary",
            name_column_index_2="type_value",
            report=True,
        )
    # Return information.
    return table


##########
# 2.5. Read and organize information for network nodes and links


def read_organize_source_supplemental_tables_for_network(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Define path to parent directory.
    path_directory_parent = os.path.join(
        paths["dock"], "out_psychiatry_biomarkers",
        "genetic_correlation", "thyroid_organization", "data"
    )
    # Define instances for iteration.
    instances = [
        {
            "group_analysis": "primary_primary_table_s2_figure_s1", # name of child directory
            "name_table": "table_psychiatry_substance_table_s2_figure_s1", # base name of file
        },
        {
            "group_analysis": "secondary_secondary_table_s3_figure_s2",
            "name_table": "table_thyroid_table_s3_figure_s2",
        },
        {
            "group_analysis": "primary_secondary_table_s4",
            "name_table": "table_psychiatry_substance_thyroid_table_s4",
        },
    ]
    # Collect information.
    tables = list()
    # Iterate on instances.
    for instance in instances:
        # Define path to child directory.
        path_directory_child = os.path.join(
            path_directory_parent, instance["group_analysis"], "pickle",
        )
        # Define path to file.
        path_file_table = os.path.join(
            path_directory_child,
            str(instance["name_table"] + "_for_supplement.pickle"),
        )
        # Read information from file.
        table_instance = pandas.read_pickle(
            path_file_table,
        )
        # Collect table.
        tables.append(table_instance)
        pass
    # Concatenate tables.
    table = pandas.concat(
        tables,
        axis="index",
        join="outer",
        ignore_index=True,
        copy=True,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        putly.print_terminal_partition(level=3)
        print(str(
            "module: psychiatry_biomarkers.genetic_correlation." +
            "thyroid_organization.py"
        ))
        print(
            "function: read_organize_source_supplemental_tables_for_network()"
        )
        print("column labels:")
        labels_columns = table.columns.to_list()
        print(labels_columns)
        putly.print_terminal_partition(level=4)
        count_columns = (table.shape[1])
        print("count of columns in table: " + str(count_columns))
        count_rows = (table.shape[0])
        print("count of rows in source table: " + str(count_rows))
        putly.print_terminal_partition(level=4)
        print("table:")
        print(table)
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return table


##########
# 2.6. Read and organize information queries


def read_organize_source_supplemental_tables_for_query(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables with entry names
            (keys) derived from original names of files

    """

    # Define path to parent directory.
    path_directory_parent = os.path.join(
        paths["dock"], "out_psychiatry_biomarkers",
        "genetic_correlation", "thyroid_organization", "data"
    )
    # Define instances for iteration.
    instances = [
        {
            "handle": "table_rg_one_one",
            "group_analysis": "primary_primary_table_s2_figure_s1", # name of child directory
            "name_table": "table_psychiatry_substance_table_s2_figure_s1", # base name of file
        },
        {
            "handle": "table_rg_two_two",
            "group_analysis": "secondary_secondary_table_s3_figure_s2",
            "name_table": "table_thyroid_table_s3_figure_s2",
        },
        {
            "handle": "table_rg_one_two",
            "group_analysis": "primary_secondary_table_s4",
            "name_table": "table_psychiatry_substance_thyroid_table_s4",
        },
    ]

    # Collect information.
    pail = dict()
    # Iterate on instances.
    for instance in instances:
        # Define path to child directory.
        path_directory_child = os.path.join(
            path_directory_parent, instance["group_analysis"], "pickle",
        )
        # Define path to file.
        path_file_table = os.path.join(
            path_directory_child,
            str(instance["name_table"] + "_for_supplement.pickle"),
        )
        # Read information from file.
        table_instance = pandas.read_pickle(
            path_file_table,
        )
        # Collect table.
        pail[instance["handle"]] = table_instance
        pass

    # Return information.
    return pail


##########
# 3. Organize genetic correlations in tables for supplements and plots


def organize_genetic_correlation_table_general(
    group_analysis=None,
    name_table=None,
    inclusion_table=None,
    inclusion_figure=None,
    groups_primary=None,
    groups_secondary=None,
    symmetry=None,
    table_studies=None,
    table_rg=None,
    paths=None,
    report=None,
):
    """
    Organize within table the information about genetic correlations from LDSC.

    arguments:
        group_analysis (str): name of child directory for analysis group of
            comparisons
        name_table (str): base name for files to which to save tables
        inclusion_table (str): name of column to use for inclusion filter
            of studies for the supplemental tables
        inclusion_figure (str): name of column to use for inclusion filter
            of studies for the figures
        groups_primary (list<str>): identifiers or names of groups of
            primary studies for which to keep information
        groups_secondary (list<str>): identifiers or names of groups of
            secondary studies for which to keep information
        symmetry (bool): whether to filter to half diagonal of a symmetrical
            matrix for plot
        table_studies (object): Pandas data-frame table of information about
            studies
        table_rg (object): Pandas data-frame table of genetic correlations
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information in table.
    table_studies = table_studies.copy(deep=True)
    table_rg = table_rg.copy(deep=True)

    ##########
    # Extract identifiers of primary and secondary studies to keep.
    table_studies_inclusion = table_studies.loc[
        (table_studies[inclusion_table] == 1), :
    ]
    table_studies_primary = table_studies_inclusion.loc[
        (table_studies_inclusion["group"].isin(groups_primary)), :
    ]
    table_studies_secondary = table_studies_inclusion.loc[
        (table_studies_inclusion["group"].isin(groups_secondary)), :
    ]
    studies_primary = copy.deepcopy(
        table_studies_primary["identifier"].to_list()
    )
    studies_secondary = copy.deepcopy(
        table_studies_secondary["identifier"].to_list()
    )

    ##########
    # Filter table's rows by primary and secondary studies in comparisons.
    # This first filter ensures that the table only includes studies for which
    # sort information is available.
    table_filter_rows = pextr.filter_table_rows_ldsc_correlation(
        table=table_rg,
        studies_primary_keep=studies_primary,
        studies_secondary_keep=studies_secondary,
        name_primary="study_primary",
        name_secondary="study_secondary",
        keep_double=False, # whether to keep double redundant pairs for symmetry
        match_redundancy=False, # whether to match redundant pairs of studies
        match_self_pair=False, # whether to match self pairs of studies
        remove_else_null=False, # whether to remove or nullify matches
        report=True,
    )

    ##########
    # Sort table's rows by primary and secondary studies in comparisons.
    # Sort table's rows before filtering table's rows so that any removal of
    # redundant pairs leaves the pair in desirable order.
    table_sort_rows = porg.sort_table_rows_primary_secondary_reference(
        table_main=table_filter_rows,
        column_main_1="study_primary",
        column_main_2="study_secondary",
        table_reference_1=table_studies_inclusion,
        table_reference_2=table_studies_inclusion,
        column_reference_1="identifier",
        column_reference_2="identifier",
        column_reference_sequence="sort", # not column 'sort_group'
        report=True,
    )

    ##########
    # Filter table's rows by primary and secondary studies in comparisons.
    # This second filter nullifies with missing values, but does not remove,
    # any self pairs and redundant pairs of the interchangeable primary and
    # secondary studies. This second filter ensures that there are not self
    # pairs or redundant pairs that would otherwise burden unnecessarily the
    # correction for multiple hypothesis testing. The subsequent operation to
    # calculate Benjamini-Hochberg false discovery rate q-values ignores any
    # comparisons with missing values.
    table_filter_rows = pextr.filter_table_rows_ldsc_correlation(
        table=table_sort_rows,
        studies_primary_keep=studies_primary,
        studies_secondary_keep=studies_secondary,
        name_primary="study_primary",
        name_secondary="study_secondary",
        keep_double=False, # whether to keep double redundant pairs for symmetry
        match_redundancy=True, # whether to match redundant pairs of studies
        match_self_pair=True, # whether to match self pairs of studies
        remove_else_null=False, # whether to remove, otherwise nullify matches
        report=True,
    )

    ##########
    # Transfer extra attributes of primary and secondary studies.
    table_extra = porg.transfer_table_rows_attributes_reference(
        table_main=table_filter_rows,
        column_main_key="study_primary",
        table_reference=table_studies_inclusion,
        column_reference_key="identifier",
        columns_reference_transfer=[
            "group", "abbreviation", "description", "sex",
        ],
        prefix_reference_main="",
        suffix_reference_main="_primary",
        report=report,
    )
    table_extra = porg.transfer_table_rows_attributes_reference(
        table_main=table_extra,
        column_main_key="study_secondary",
        table_reference=table_studies_inclusion,
        column_reference_key="identifier",
        columns_reference_transfer=[
            "group", "abbreviation", "description", "sex",
        ],
        prefix_reference_main="",
        suffix_reference_main="_secondary",
        report=report,
    )

    ##########
    # Calculate Benjamini-Hochberg q-values for False-Discovery Rate (FDR).
    # Calculate q-values across all comparisons in table.
    # FDR 5% (q <= 0.05).
    table_q = pdesc.calculate_table_false_discovery_rate_q_values(
        threshold=0.05, # alpha; family-wise error rate
        name_column_p_value="p_value_ldsc",
        name_column_q_value="q_value_ldsc",
        name_column_significance="q_significance_ldsc",
        table=table_extra,
    )
    table_q = pdesc.calculate_table_false_discovery_rate_q_values(
        threshold=0.05, # alpha; family-wise error rate
        name_column_p_value="p_value_not_zero",
        name_column_q_value="q_value_not_zero",
        name_column_significance="q_significance_not_zero",
        table=table_q,
    )
    table_q = pdesc.calculate_table_false_discovery_rate_q_values(
        threshold=0.05, # alpha; family-wise error rate
        name_column_p_value="p_value_less_one",
        name_column_q_value="q_value_less_one",
        name_column_significance="q_significance_less_one",
        table=table_q,
    )

    ##########
    # Indicate analysis group.
    table_group = table_q.copy(deep=True)
    table_group["group_analysis"] = group_analysis

    ##########
    # Copy information in table.
    table_general = table_group.copy(deep=True)

    ##########
    # Organize information in table.
    table_general.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )

    ##########
    # Return information.
    return table_general


def organize_genetic_correlation_table_supplement(
    group_analysis=None,
    name_table=None,
    inclusion=None,
    groups_primary=None,
    groups_secondary=None,
    symmetry=None,
    table_studies=None,
    table_rg=None,
    paths=None,
    report=None,
):
    """
    Organize within table the information about genetic correlations from LDSC.

    This function applies a filter to the table for specific, relevant studies
    using the same parameters of "inclusion" and "table_studies" as in the
    previous function that filtered "table_general". These filters in this
    function will have the additional effect of removing the rows for self
    pairs and redundant pairs which the previous filter had nullified by
    filling with missing values.

    arguments:
        group_analysis (str): name of child directory for analysis group of
            comparisons
        name_table (str): base name for files to which to save tables
        inclusion (str): name of column to use for inclusion filter
        groups_primary (list<str>): identifiers or names of groups of
            primary studies for which to keep information
        groups_secondary (list<str>): identifiers or names of groups of
            secondary studies for which to keep information
        symmetry (bool): whether to filter to half diagonal of a symmetrical
            matrix for plot
        table_studies (object): Pandas data-frame table of information about
            studies
        table_rg (object): Pandas data-frame table of genetic correlations
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information in table.
    table_studies = table_studies.copy(deep=True)
    table_rg = table_rg.copy(deep=True)

    ##########
    # Extract identifiers of primary and secondary studies to keep.
    table_studies_inclusion = table_studies.loc[
        (table_studies[inclusion] == 1), :
    ]
    table_studies_primary = table_studies_inclusion.loc[
        (table_studies_inclusion["group"].isin(groups_primary)), :
    ]
    table_studies_secondary = table_studies_inclusion.loc[
        (table_studies_inclusion["group"].isin(groups_secondary)), :
    ]
    studies_primary = copy.deepcopy(
        table_studies_primary["identifier"].to_list()
    )
    studies_secondary = copy.deepcopy(
        table_studies_secondary["identifier"].to_list()
    )

    ##########
    # Filter table's rows by primary and secondary studies in comparisons.
    # For text tables:
    # - match_redundancy=True
    # - match_self_pair=True
    # - remove_else_null=True
    table_filter_rows = pextr.filter_table_rows_ldsc_correlation(
        table=table_rg,
        studies_primary_keep=studies_primary,
        studies_secondary_keep=studies_secondary,
        name_primary="study_primary",
        name_secondary="study_secondary",
        keep_double=False, # whether to keep double redundant pairs for symmetry
        match_redundancy=True, # whether to match redundant pairs of studies
        match_self_pair=True, # whether to match self pairs of studies
        remove_else_null=True, # whether to remove or nullify matches
        report=True,
    )

    ##########
    # Sort table's rows by primary and secondary studies in comparisons.
    # Sort table's rows before filtering table's rows so that any removal of
    # redundant pairs leaves the pair in desirable order.
    table_sort_rows = porg.sort_table_rows_primary_secondary_reference(
        table_main=table_filter_rows,
        column_main_1="study_primary",
        column_main_2="study_secondary",
        table_reference_1=table_studies_inclusion,
        table_reference_2=table_studies_inclusion,
        column_reference_1="identifier",
        column_reference_2="identifier",
        column_reference_sequence="sort", # not column 'sort_group'
        report=True,
    )

    ##########
    # Filter and sort table's columns.
    columns_sequence = pextr.define_genetic_correlation_table_column_sequence()
    table_filter_sort_columns = porg.filter_sort_table_columns(
        table=table_sort_rows,
        columns_sequence=columns_sequence,
        report=report,
    )

    ##########
    # Copy information in table.
    table_supplement = table_filter_sort_columns.copy(deep=True)

    ##########
    # Organize information in table.
    table_supplement.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_supplement.set_index(
        [
            "group_analysis",
            "abbreviation_primary",
            "abbreviation_secondary",
        ],
        append=False,
        drop=True,
        inplace=True,
    )

    ##########
    # Return information.
    return table_supplement


def simplify_transform_long_table_plot_symmetrical(
    studies_primary_keep=None,
    studies_secondary_keep=None,
    table_rg=None,
    table_studies_inclusion=None,
    paths=None,
    report=None,
):
    """
    Beginning with a symmetrical table of genetic correlations, this function
    simplifies contents, filters to half diagonal of corresponding matrices of
    values, and transforms to long format.

    To enable the filter to the half diagonal of the otherwise symmetrical
    matrix, it is necessary to accommodate the prior filters in preparation of
    "table_general" that nullified by filling missing values for any self pairs
    or redundant pairs of primary and secondary studies. This function applies
    a clever patch to fill missing values from the reciprocal, interchangeable
    combination of primary and secondary studies.

    arguments:
        studies_primary_keep (list<str>): identifiers or names of primary
            studies for which to keep information in table
        studies_secondary_keep (list<str>): identifiers or names of secondary
            studies for which to keep information in table
        table_rg (object): Pandas data-frame table of genetic correlations
        table_studies_inclusion (object): Pandas data-frame table of
            information about studies
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information in table.
    table_rg = table_rg.copy(deep=True)
    table_studies_inclusion = table_studies_inclusion.copy(deep=True)

    ##########
    # Filter table's rows by primary and secondary studies in comparisons.
    # For symmetrical plot tables:
    # - match_redundancy=False
    # - match_self_pair=True
    # - remove_else_null=False
    # For symmetrical tables with self pairs, this second filter assigns
    # missing values to the self pairs but does not remove redundant pairs.
    # A subsequent filter to half of the symmetrical diagonal will remove
    # the redundant pairs.
    table_filter_rows = pextr.filter_table_rows_ldsc_correlation(
        table=table_rg,
        studies_primary_keep=studies_primary_keep,
        studies_secondary_keep=studies_secondary_keep,
        name_primary="study_primary",
        name_secondary="study_secondary",
        keep_double=False, # whether to keep double redundant pairs for symmetry
        match_redundancy=False, # whether to match redundant pairs of studies
        match_self_pair=False, # whether to match self pairs of studies
        remove_else_null=False, # whether to remove, otherwise nullify matches
        report=True,
    )

    ##########
    # To ensure that non-missing values are available for the subsequent filter
    # to the half diagonal matrix, fill missing values with any non-missing
    # values from the reciprocal, interchangeable pair of primary and secondary
    # studies.
    table_fill = pextr.fill_missing_from_reciprocal_interchangeable_study_pair(
        table=table_filter_rows,
        name_primary="study_primary",
        name_secondary="study_secondary",
        report=False,
    )

    ##########
    # Sort table's rows by primary and secondary studies in comparisons.
    # Sort table's rows before filtering table's rows so that any removal of
    # redundant pairs leaves the pair in desirable order.
    table_sort_rows = porg.sort_table_rows_primary_secondary_reference(
        table_main=table_fill,
        column_main_1="study_primary",
        column_main_2="study_secondary",
        table_reference_1=table_studies_inclusion,
        table_reference_2=table_studies_inclusion,
        column_reference_1="identifier",
        column_reference_2="identifier",
        column_reference_sequence="sort", # not column 'sort_group'
        report=True,
    )

    ##########
    # The next few steps transform the table to partial wide format to filter
    # half diagonal of matrix.

    ##########
    # Simplify content of table.
    # Transform table to long format.
    table_long = pextr.simplify_transform_genetic_correlation_table_long(
        table_rg=table_sort_rows,
        q_values=True,
        report=report,
    )
    # Remove temporarily the index level for analysis group.
    # There ought to be only a single value of the analysis group.
    group_analysis = table_long.index.get_level_values("group_analysis")[0]
    table_long.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_long.drop(
        labels=["group_analysis",],
        axis="columns",
        inplace=True
    )
    table_long.set_index(
        [
            "abbreviation_primary",
            "abbreviation_secondary",
            "type_value",
        ],
        append=False,
        drop=True,
        inplace=True,
    )

    ##########
    # Transform table to partial wide format.
    # This partial wide format corresponds to the symmetrical matrix for
    # graphical representation in a heatmap plot.
    table_wide_partial = (
        porg.transform_table_triple_quadruple_index_long_to_wide_partial(
            table_long=table_long,
            column_index_pivot="abbreviation_primary",
            columns_index_stay=["abbreviation_secondary", "type_value",],
            column_value_long="value",
            report=report,
    ))

    ##########
    # Filter to the diagonal lower half of a table representing a symmetrical
    # matrix.
    table_wide_half = (
        porg.filter_symmetrical_table_half_diagonal_two_value_types(
            table=table_wide_partial,
            name_index_type_value="type_value",
            types_value=["signal", "p_value", "q_value"],
            report=report,
        )
    )

    ##########
    # Transform table from partial wide to full long format.
    # Pandas dataframe methods "stack", "melt", and "wide_to_long", can all be
    # useful in this context.
    # Method "stack" converts to a multi-index series when the column index only
    # has a single level.
    # Method "wide_to_long" assumes that the information about multiple levels
    # in the column index is stored in delimited strings of compound column
    # names.
    if False:
        table_long = table_simple.stack(
            level=name_row_index,
            #future_stack=True,
        )
    table_long = table_wide_half.melt(
        id_vars=None,
        value_vars=None,
        var_name="abbreviation_primary",
        value_name="value",
        ignore_index=False,
    )

    # Introduce analysis group.
    table_long["group_analysis"] = group_analysis
    # Organize information in table.
    table_long.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_long.set_index(
        [
            "group_analysis",
            "abbreviation_primary",
            "abbreviation_secondary",
            "type_value",
        ],
        append=False,
        drop=True,
        inplace=True,
    )

    ##########
    # Return information.
    return table_long


def simplify_transform_long_table_plot_asymmetrical(
    studies_primary_keep=None,
    studies_secondary_keep=None,
    table_rg=None,
    table_studies_inclusion=None,
    paths=None,
    report=None,
):
    """
    Beginning with an asymmetrical table of genetic correlations, this function
    simplifies contents and transforms to long format.

    This function applies a filter to the table for specific, relevant studies
    using the same parameters of "inclusion" and "table_studies" as in the
    previous function that filtered "table_general". These filters in this
    function will have the additional effect of removing the rows for self
    pairs and redundant pairs which the previous filter had nullified by
    filling with missing values.

    arguments:
        studies_primary_keep (list<str>): identifiers or names of primary
            studies for which to keep information in table
        studies_secondary_keep (list<str>): identifiers or names of secondary
            studies for which to keep information in table
        table_rg (object): Pandas data-frame table of genetic correlations
        table_studies_inclusion (object): Pandas data-frame table of
            information about studies
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information in table.
    table_rg = table_rg.copy(deep=True)
    table_studies_inclusion = table_studies_inclusion.copy(deep=True)

    ##########
    # Filter table's rows by primary and secondary studies in comparisons.
    # For asymmetrical plot tables:
    # - match_redundancy=True
    # - match_self_pair=True
    # - remove_else_null=True
    # For asymmetrical tables without self pairs, this second filter removes
    # redundant pairs. The prior sort of primary and secondary studies ought
    # to preserve the desirable designations of primary and secondary
    # studies.
    table_filter_rows = pextr.filter_table_rows_ldsc_correlation(
        table=table_rg,
        studies_primary_keep=studies_primary_keep,
        studies_secondary_keep=studies_secondary_keep,
        name_primary="study_primary",
        name_secondary="study_secondary",
        keep_double=False, # whether to keep double redundant pairs for symmetry
        match_redundancy=True, # whether to match redundant pairs of studies
        match_self_pair=True, # whether to match self pairs of studies
        remove_else_null=True, # whether to remove or nullify matches
        report=True,
    )

    ##########
    # Sort table's rows by primary and secondary studies in comparisons.
    # Sort table's rows before filtering table's rows so that any removal of
    # redundant pairs leaves the pair in desirable order.
    table_sort_rows = porg.sort_table_rows_primary_secondary_reference(
        table_main=table_filter_rows,
        column_main_1="study_primary",
        column_main_2="study_secondary",
        table_reference_1=table_studies_inclusion,
        table_reference_2=table_studies_inclusion,
        column_reference_1="identifier",
        column_reference_2="identifier",
        column_reference_sequence="sort", # not column 'sort_group'
        report=True,
    )

    ##########
    # Simplify content of table.
    # Transform table to long format.
    table_long = pextr.simplify_transform_genetic_correlation_table_long(
        table_rg=table_sort_rows,
        q_values=True,
        report=report,
    )

    ##########
    # Return information.
    return table_long


def organize_genetic_correlation_table_plot(
    group_analysis=None,
    name_table=None,
    inclusion=None,
    groups_primary=None,
    groups_secondary=None,
    symmetry=None,
    table_studies=None,
    table_rg=None,
    paths=None,
    report=None,
):
    """
    Organize within table the information about genetic correlations from LDSC
    for reporting in graphical plots.

    The product table will be in square format with floating-point values of
    signals, p-values, and q-values organized within two levels of categorical
    indices across both columns and rows. There will be a subsequent
    transposition of this table before passing it to the plotting function.

    study_2nd           study_2_a                    study_2_b   ...
    type_value          signal   p_value   q_value   signal   p_value   q_value
    group   study_1st
    a       study_1_a   -0.15    0.001     0.01      -0.15    0.001     0.01
    a       study_1_b   -0.2     0.001     0.01      -0.2     0.001     0.01
    a       study_1_c   -0.25    0.001     0.01      -0.25    0.001     0.01

    arguments:
        group_analysis (str): name of child directory for analysis group of
            comparisons
        name_table (str): base name for files to which to save tables
        inclusion (str): name of column to use for inclusion filter
        groups_primary (list<str>): identifiers or names of groups of
            primary studies for which to keep information
        groups_secondary (list<str>): identifiers or names of groups of
            secondary studies for which to keep information
        symmetry (bool): whether to filter to half diagonal of a symmetrical
            matrix for plot
        table_studies (object): Pandas data-frame table of information about
            studies
        table_rg (object): Pandas data-frame table of genetic correlations
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information in table.
    table_studies = table_studies.copy(deep=True)
    table_rg = table_rg.copy(deep=True)

    ##########
    # Extract identifiers of primary and secondary studies to keep.
    table_studies_inclusion = table_studies.loc[
        (table_studies[inclusion] == 1), :
    ]
    table_studies_primary = table_studies_inclusion.loc[
        (table_studies_inclusion["group"].isin(groups_primary)), :
    ]
    table_studies_secondary = table_studies_inclusion.loc[
        (table_studies_inclusion["group"].isin(groups_secondary)), :
    ]
    studies_primary = copy.deepcopy(
        table_studies_primary["identifier"].to_list()
    )
    studies_secondary = copy.deepcopy(
        table_studies_secondary["identifier"].to_list()
    )

    ##########
    # Simplify content of table and transform to long format.
    if symmetry:
        table_long = simplify_transform_long_table_plot_symmetrical(
            studies_primary_keep=studies_primary,
            studies_secondary_keep=studies_secondary,
            table_rg=table_rg,
            table_studies_inclusion=table_studies_inclusion,
            paths=paths,
            report=report,
        )
        pass
    else:
        table_long = simplify_transform_long_table_plot_asymmetrical(
            studies_primary_keep=studies_primary,
            studies_secondary_keep=studies_secondary,
            table_rg=table_rg,
            table_studies_inclusion=table_studies_inclusion,
            paths=paths,
            report=report,
        )
        pass

    ##########
    # Calculate Benjamini-Hochberg q-values to control the False-Discovery Rate
    # (FDR) for testing multiple hypotheses.
    # Calculate q-values across all comparisons in table.
    # FDR 5% (q <= 0.05).
    # The tables for plots and figures now use the same Benjamini-Hochberg
    # false discovery rate q-values as calculated for "table_general" and used
    # in the supplemental table (TCW; 9 August 2024).
    #table_q = pdesc.calculate_table_long_false_discovery_rate_q_values(
    #    table=table_long,
    #    name_index_type_value="type_value",
    #    type_p_value="p_value",
    #    threshold=0.05,
    #    names_indices_rows=[
    #        "group_analysis",
    #        "abbreviation_primary",
    #        "abbreviation_secondary",
    #        "type_value",
    #    ],
    #    report=report,
    #)

    ##########
    # Transform table to wide format.
    # Across columns: abbreviation_secondary, type_value
    # Across rows: group_analysis, abbreviation_primary
    table_wide = porg.transform_table_quadruple_index_long_to_wide_square(
        table_long=table_long,
        columns_index_pivot=["abbreviation_secondary", "type_value",],
        columns_index_stay=["group_analysis", "abbreviation_primary",],
        column_value_long="value",
        report=report,
    )

    ##########
    # Copy information in table.
    table_plot = table_wide.copy(deep=True)

    ##########
    # Return information.
    return table_plot


def control_prepare_genetic_correlation_table_supplement_plot(
    instance=None,
    parameters=None,
):
    """
    Control procedure to organize within tables the information about genetic
    correlations from LDSC.

    Note: TCW; 9 August 2024

    For organization of "table_general", filter table's rows to include only
    select, relevant studies and to exclude self pairs (same primary and
    secondary studies) and redundant pairs (interchangeable primary and
    secondary studies). These filters will ensure that both the supplemental
    tables and the heatmap plots only represent the relevant genetic
    correlations. Hence, it will only be necessary to calculate
    Benjamini-Hochberg false discovery rate q-values once for "table_general"
    and then use these same q-values in both the supplemental tables and the
    heatmap plots. To allow the symmetrical heatmap plot (same studies across
    horizontal and vertical axes) to have appropriate values when reducing to
    the half diagonal, there is a special function for preparation of this
    table that fills missing values using the reciprocal equivalent of
    interchangeable primary and secondary studies.

    Notice that the "inclusion_figure" studies for the plots and figures are a
    subset of the larger "inclusion_table" studies for the supplemental table.
    The supplemental table will include all of the studies in the narrative,
    whereas the plots and figures will only include some of these. Importantly,
    the Benjamini-Hochberg false discovery rate q-values will represent the
    all inclusive set of studies and comparisons that are in the supplemental
    table. The plots and figures will use these same q-values.

    arguments:
        instance (dict): parameters specific to current instance
            group_analysis (str): name of child directory for analysis group of
                comparisons
            name_table (str): base name for files to which to save tables
            inclusion_table (str): name of column to use for inclusion filter
                of studies for the supplemental tables
            inclusion_figure (str): name of column to use for inclusion filter
                of studies for the figures
            groups_primary (list<str>): identifiers or names of groups of
                primary studies for which to keep information
            groups_secondary (list<str>): identifiers or names of groups of
                secondary studies for which to keep information
            symmetry (bool): whether to filter to half diagonal of a symmetrical
                matrix for plot
        parameters (dict): parameters common to all instances
            table_studies (object): Pandas data-frame table of information about
                studies
            table_rg (object): Pandas data-frame table of genetic correlations
            paths : (dict<str>): collection of paths to directories for procedure's
                files
            report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Extract parameters.
    # Extract parameters specific to each instance.
    group_analysis = instance["group_analysis"]
    name_table = instance["name_table"]
    inclusion_table = instance["inclusion_table"]
    inclusion_figure = instance["inclusion_figure"]
    groups_primary = instance["groups_primary"]
    groups_secondary = instance["groups_secondary"]
    symmetry = instance["symmetry"]
    # Extract parameters common across all instances.
    table_studies = parameters["table_studies"]
    table_rg = parameters["table_rg"]
    paths = parameters["paths"]
    report = parameters["report"]

    ##########
    # Filter and sort information in table for general use.
    table_general = organize_genetic_correlation_table_general(
        group_analysis=group_analysis,
        name_table=name_table,
        inclusion_table=inclusion_table,
        inclusion_figure=inclusion_figure,
        groups_primary=groups_primary,
        groups_secondary=groups_secondary,
        symmetry=symmetry,
        table_studies=table_studies,
        table_rg=table_rg,
        paths=paths,
        report=report,
    )

    ##########
    # Organize and transform information in table to a format for reporting as a
    # text supplemental table.
    table_supplement = organize_genetic_correlation_table_supplement(
        group_analysis=group_analysis,
        name_table=name_table,
        inclusion=inclusion_table,
        groups_primary=groups_primary,
        groups_secondary=groups_secondary,
        symmetry=symmetry,
        table_studies=table_studies,
        table_rg=table_general,
        paths=paths,
        report=report,
    )

    ##########
    # Organize and transform information in table to a format for reporting in a
    # heatmap plot.
    table_plot = organize_genetic_correlation_table_plot(
        group_analysis=group_analysis,
        name_table=name_table,
        inclusion=inclusion_figure,
        groups_primary=groups_primary,
        groups_secondary=groups_secondary,
        symmetry=symmetry,
        table_studies=table_studies,
        table_rg=table_general,
        paths=paths,
        report=report,
    )

    ##########
    # Initialize child directory for analysis group.
    paths = initialize_directory_group_analysis(
        group_analysis=group_analysis,
        paths=paths,
        restore=True,
    )

    ##########
    # Collect information.
    # Collections of files.
    pail_write_files = dict()
    pail_write_files[str(name_table + "_for_supplement")] = table_supplement
    pail_write_files[str(name_table + "_for_plot")] = table_plot
    # Collections of directories.
    pail_write_directories_text = dict()
    pail_write_directories_pickle = dict()
    pail_write_directories_text["text"] = pail_write_files
    pail_write_directories_pickle["pickle"] = pail_write_files

    ##########
    # Write product information to file.
    putly.write_tables_to_file_in_child_directories(
        pail_write=pail_write_directories_text,
        path_directory_parent=paths["out_data_group_analysis"],
        reset_index_rows=False,
        write_index_rows=True,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )
    putly.write_tables_to_file_in_child_directories(
        pail_write=pail_write_directories_pickle,
        path_directory_parent=paths["out_data_group_analysis"],
        reset_index_rows=False,
        write_index_rows=None,
        write_index_columns=None,
        type="pickle",
        delimiter=None,
        suffix=".pickle",
    )
    pass


def control_prepare_genetic_correlation_tables_supplement_plot(
    table_rg=None,
    paths=None,
    report=None,
):
    """
    Control procedure to organize within tables the information about genetic
    correlations from LDSC.

    arguments:
        table_rg (object): Pandas data-frame table of genetic correlations
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports


    raises:

    returns:

    """

    # Read source information from file.
    # Organize source information within tables.
    table_studies = read_source_parameter_studies_polish(
        paths=paths,
        report=report,
    )
    # Collect parameters common across all instances.
    parameters = dict()
    parameters["table_studies"] = table_studies
    parameters["table_rg"] = table_rg
    parameters["paths"] = paths
    parameters["report"] = report

    # Collect parameters specific to each instance.
    instances = [
        {
            "group_analysis": "primary_primary_table_s2_figure_s1", # name of child directory
            "name_table": "table_psychiatry_substance_table_s2_figure_s1", # base name of file
            "inclusion_table": "inclusion_thyroid_table",
            "inclusion_figure": "inclusion_thyroid_figure",
            "groups_primary": ["psychiatry","substance",],
            "groups_secondary": ["psychiatry","substance",],
            "symmetry": True,
        },
        {
            "group_analysis": "secondary_secondary_table_s3_figure_s2",
            "name_table": "table_thyroid_table_s3_figure_s2",
            "inclusion_table": "inclusion_thyroid_table",
            "inclusion_figure": "inclusion_thyroid_figure",
            "groups_primary": ["thyroid",],
            "groups_secondary": ["thyroid",],
            "symmetry": True,
        },
        {
            "group_analysis": "primary_secondary_table_s4",
            "name_table": "table_psychiatry_substance_thyroid_table_s4",
            "inclusion_table": "inclusion_thyroid_table",
            "inclusion_figure": "inclusion_thyroid_figure",
            "groups_primary": ["psychiatry","substance",],
            "groups_secondary": ["thyroid",],
            "symmetry": False,
        },
        {
            "group_analysis": "primary_secondary_figure_1", # must use same parameters for q-values in supplemental table
            "name_table": "table_psychiatry_substance_thyroid_figure_1",
            "inclusion_table": "inclusion_thyroid_table",
            "inclusion_figure": "inclusion_thyroid_figure_1",
            "groups_primary": ["psychiatry","substance",],
            "groups_secondary": ["thyroid",],
            "symmetry": False,
        },
        {
            "group_analysis": "primary_secondary_thyroid_figure_2", # must use same parameters for q-values in supplemental table
            "name_table": "table_psychiatry_substance_thyroid_figure_2",
            "inclusion_table": "inclusion_thyroid_table",
            "inclusion_figure": "inclusion_thyroid_figure_2",
            "groups_primary": ["psychiatry","substance",],
            "groups_secondary": ["thyroid",],
            "symmetry": False,
        },
        {
            "group_analysis": "primary_secondary_table_sex_alcohol_tobacco", # name of child directory
            "name_table": "table_sex_substance", # base name of file
            "inclusion_table": "inclusion_sex_alcohol_tobacco_table",
            "inclusion_figure": "inclusion_sex_alcohol_tobacco_table",
            "groups_primary": ["substance",],
            "groups_secondary": ["hormone_sex",],
            "symmetry": False,
        },
        {
            "group_analysis": "gonad_psychiatry_substance", # name of child directory
            "name_table": "table_gonad_psychiatry_substance", # base name of file
            "inclusion_table": "inclusion_gonad",
            "inclusion_figure": "inclusion_gonad",
            "groups_primary": ["psychiatry", "substance",],
            "groups_secondary": ["hormone_sex",],
            "symmetry": False,
        },
        {
            "group_analysis": "gonad_psychiatry", # name of child directory
            "name_table": "table_gonad_psychiatry", # base name of file
            "inclusion_table": "inclusion_gonad",
            "inclusion_figure": "inclusion_gonad",
            "groups_primary": ["psychiatry",],
            "groups_secondary": ["hormone_sex",],
            "symmetry": False,
        },
        {
            "group_analysis": "gonad_substance", # name of child directory
            "name_table": "table_gonad_substance", # base name of file
            "inclusion_table": "inclusion_gonad",
            "inclusion_figure": "inclusion_gonad",
            "groups_primary": ["substance",],
            "groups_secondary": ["hormone_sex",],
            "symmetry": False,
        },
        #{
        #    "group_analysis": "secondary_secondary_sex",
        #    "name_table": "table_sex_hormone",
        #    "inclusion": "inclusion_sex",
        #    "groups_primary": ["hormone_sex",],
        #    "groups_secondary": ["hormone_sex",],
        #    "symmetry": True,
        #},
        #{
        #    "group_analysis": "primary_secondary_sex",
        #    "name_table": "table_psychiatry_substance_sex_hormone",
        #    "inclusion": "inclusion_sex",
        #    "groups_primary": ["psychiatry","substance",],
        #    "groups_secondary": ["hormone_sex",],
        #    "symmetry": False,
        #},

    ]

    # Execute procedure iteratively with parallelization across instances.
    if True:
        prall.drive_procedure_parallel(
            function_control=(
                control_prepare_genetic_correlation_table_supplement_plot
            ),
            instances=instances,
            parameters=parameters,
            cores=2,
            report=True,
        )
    else:
        # Execute procedure directly for testing.
        control_prepare_genetic_correlation_table_supplement_plot(
            instance=instances[0],
            parameters=parameters,
        )
    pass


def test():

    ##################################################
    # TEST
    ##################################################
    ##########
    table_test = table_general
    # Collect information.
    # Collections of files.
    pail_write_files = dict()
    pail_write_files[str("test_table_check")] = table_test
    ##########
    # Write product information to file.
    putly.write_tables_to_file(
        pail_write=pail_write_files,
        path_directory=paths["out_test"],
        reset_index_rows=False,
        write_index_rows=True,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )

    ##################################################
    pass


##########
# 4. Network
# Genome-Wide Association Studies (GWAS) as Nodes and Genetic Correlations
# as Links between them.


def control_prepare_genetic_correlation_network_nodes_links(
    paths=None,
    report=None,
):
    """
    Control procedure to organize within tables the information about genetic
    correlations from LDSC.

    arguments:
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports


    raises:

    returns:

    """

    # Organize parameters.
    group_analysis = "network_links"
    name_table_links = "table_network_links"
    name_table_nodes = "table_network_nodes"
    groups_primary = ["psychiatry","substance", "thyroid",]
    groups_secondary = ["psychiatry","substance", "thyroid",]
    symmetry = False

    # Read source information from file.
    table_studies = read_source_parameter_studies_polish(
        paths=paths,
        report=report,
    )
    table_rg = read_organize_source_supplemental_tables_for_network(
        paths=paths,
        report=report,
    )

    # Filter rows in table by selection of studies.
    inclusion = "inclusion_thyroid_figure"
    # Extract identifiers of primary and secondary studies to keep.
    table_studies_inclusion = table_studies.loc[
        (table_studies[inclusion] == 1), :
    ]
    table_studies_primary = table_studies_inclusion.loc[
        (table_studies_inclusion["group"].isin(groups_primary)), :
    ]
    table_studies_secondary = table_studies_inclusion.loc[
        (table_studies_inclusion["group"].isin(groups_secondary)), :
    ]
    studies_primary = copy.deepcopy(
        table_studies_primary["identifier"].to_list()
    )
    studies_secondary = copy.deepcopy(
        table_studies_secondary["identifier"].to_list()
    )
    # Filter table's rows by primary and secondary studies in comparisons.
    # For text tables:
    # - match_redundancy=True
    # - match_self_pair=True
    # - remove_else_null=True
    table_rg = pextr.filter_table_rows_ldsc_correlation(
        table=table_rg,
        studies_primary_keep=studies_primary,
        studies_secondary_keep=studies_secondary,
        name_primary="study_primary",
        name_secondary="study_secondary",
        keep_double=False, # whether to keep double redundant pairs for symmetry
        match_redundancy=True, # whether to match redundant pairs of studies
        match_self_pair=True, # whether to match self pairs of studies
        remove_else_null=True, # whether to remove or nullify matches
        report=True,
    )

    # Filter rows in table for non-missing values across relevant columns.
    table_rg.dropna(
        axis="index",
        how="any",
        subset=[
            "correlation",
            "correlation_error",
            "p_value_ldsc",
            "q_value_ldsc",
        ],
        inplace=True,
    )
    # Filter rows in table by applying a threshold on the q-value.
    table_rg = table_rg.loc[
        (
            (table_rg["q_value_ldsc"] < 0.05)
        ), :
    ].copy(deep=True)

    ##########
    # Initialize child directory for analysis group.
    paths = initialize_directory_group_analysis(
        group_analysis=group_analysis,
        paths=paths,
        restore=True,
    )

    ##########
    # Collect information.
    # Collections of files.
    pail_write_files = dict()
    pail_write_files[str(name_table_nodes)] = table_studies
    pail_write_files[str(name_table_links)] = table_rg
    # Collections of directories.
    pail_write_directories_text = dict()
    pail_write_directories_text["text"] = pail_write_files

    ##########
    # Write product information to file.
    putly.write_tables_to_file_in_child_directories(
        pail_write=pail_write_directories_text,
        path_directory_parent=paths["out_data_group_analysis"],
        reset_index_rows=False,
        write_index_rows=True,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )
    pass


# This function will probably become obsolete... TCW; 19 September 2024
def control_prepare_genetic_correlation_network_nodes_links_from_scratch(
    table_rg=None,
    paths=None,
    report=None,
):
    """
    Control procedure to organize within tables the information about genetic
    correlations from LDSC.

    arguments:
        table_rg (object): Pandas data-frame table of genetic correlations
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports


    raises:

    returns:

    """

    # Organize parameters.
    group_analysis = "network_links"
    name_table = "table_network_links"
    groups_primary = ["psychiatry","substance", "thyroid",]
    groups_secondary = ["psychiatry","substance", "thyroid",]
    symmetry = False

    # Read source information from file.
    # Organize source information within tables.
    table_studies = read_source_parameter_studies_polish(
        paths=paths,
        report=report,
    )

    # TODO: organize the table of network links

    # TODO: TCW; 11 June 2024
    # filter to keep all relevant primary + secondary studies
    # filter to remove all self pairs
    # filter to remove any redundant first<-->second, second<-->first pairs
    # calculate q-values
    # filter links by q-values


    ##########
    # Extract identifiers of primary and secondary studies to keep.
    table_studies_inclusion = table_studies.loc[
        (table_studies["inclusion"] == 1), :
    ]
    table_studies_primary = table_studies_inclusion.loc[
        (table_studies_inclusion["group"].isin(groups_primary)), :
    ]
    table_studies_secondary = table_studies_inclusion.loc[
        (table_studies_inclusion["group"].isin(groups_secondary)), :
    ]
    studies_primary = copy.deepcopy(
        table_studies_primary["identifier"].to_list()
    )
    studies_secondary = copy.deepcopy(
        table_studies_secondary["identifier"].to_list()
    )

    ##########
    # Filter table's rows by primary and secondary studies in comparisons.
    # This filter ensures that the table only includes studies from relevant
    # groups, regardless of whether those studies are designated as primary or
    # secondary studies in each comparison.
    # This filter also ensures that there are not redundant pairs or self
    # pairs that would otherwise burden unnecessarily the correction for
    # multiple hypothesis testing.
    # For consideration of redundancy in the pairs of studies, primary and
    # secondary studies are considered to be interchangeable.
    # For text tables:
    # - match_redundancy=True
    # - match_self_pair=True
    # - remove_else_null=True
    table_filter_rows = pextr.filter_table_rows_ldsc_correlation(
        table=table_rg,
        studies_primary_keep=studies_primary,
        studies_secondary_keep=studies_secondary,
        name_primary="study_primary",
        name_secondary="study_secondary",
        match_redundancy=True, # whether to match redundant pairs of studies
        match_self_pair=True, # whether to match self pairs of studies
        remove_else_null=True, # whether to remove, otherwise nullify matches
        report=True,
    )



    ##########
    # Copy information in table.
    #table_links = pandas.DataFrame()
    table_links = table_general.copy(deep=True)

    ##########
    # Initialize child directory for analysis group.
    paths = initialize_directory_group_analysis(
        group_analysis=group_analysis,
        paths=paths,
        restore=True,
    )

    ##########
    # Collect information.
    # Collections of files.
    pail_write_files = dict()
    pail_write_files[str(name_table)] = table_links
    # Collections of directories.
    pail_write_directories_text = dict()
    pail_write_directories_text["text"] = pail_write_files

    ##########
    # Write product information to file.
    putly.write_tables_to_file_in_child_directories(
        pail_write=pail_write_directories_text,
        path_directory_parent=paths["out_data_group_analysis"],
        reset_index_rows=False,
        write_index_rows=True,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )
    pass



##########
# 5. Plot charts


def create_write_figure_heat_map(
    table=None,
    name_figure=None,
    path_directory=None,
    report=None,
):
    """
    Creation dot plot and write to file.

    Assume that rows for records in the source table are already in proper sort
    order corresponding to the categorical labels on the abscissa (horizontal
    axis).

    arguments:
        table (object): Pandas data-frame table in long format with
            floating-point values of signal and p-values.
        name_figure (str): name of figure to use in file name
        path_directory (str): path to parent directory within which to write a
            file for figure
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_heat_map_few_signal_significance_labels(
        table=table,
        transpose_table=True,
        index_group_columns="abbreviation_secondary",
        index_group_rows="abbreviation_primary",
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=-1.0, # -1.0
        value_maximum=1.0, # 1.0
        significance_p=True,
        significance_q=True,
        thresholds_p=[0.001, 0.00001,],
        #thresholds_q=[0.05, 0.01,],
        thresholds_q=[0.05,],
        show_legend=False, # whether to show legend on individual figures
        show_scale_bar=True, # whether to show scale bar on individual figures
        #label_legend=str(
        #    "Labels:\n"
        #    + str("  " + "$\u2020$: p < 0.05" + "  \n") # dagger U+2020
        #    + str("  " + "$\u2021$: p < 0.01" + "  \n") # double dagger U+2021
        #    + str("  " + "  \n") # blank line
        #    + str("  " + "$\u002A$: q < 0.05" + "  \n") # asterisk U+002A
        #    + str("  " + "$\u2051$: q < 0.01") # double asterisk U+2051
        #),
        label_legend=str(
            "Labels:\n"
            + str("  " + "  \n") # blank line
            + str("  " + "$\u2020$: p < 1E-3" + "  \n") # dagger U+2020
            + str("  " + "$\u2021$: p < 1E-5" + "  \n") # double dagger U+2021
            + str("  " + "  \n") # blank line
            + str("  " + "$\u002A$: q < 0.05") # asterisk U+002A
        ),
        labels_ordinate_categories=[""],
        labels_abscissa_categories=[""],
        title_chart_top_right="",
        title_ordinate="",
        title_abscissa="",
        title_bar="Genetic Correlations (rg)", # Genetic Correlations (rg)
        size_label_sig_p="twelve", # multi-panel: twelve; individual: thirteen
        size_label_sig_q="twelve", # multi-panel: twelve; individual: thirteen
        size_title_ordinate="eight", # ten
        size_title_abscissa="eight", # ten
        size_label_ordinate="eleven", # multi-panel: ten; individual: twelve
        size_label_abscissa="eleven", # multi-panel: ten; individual: twelve
        size_label_legend="twelve", # twelve
        size_title_bar="twelve", # twelve
        size_label_bar="thirteen", # thirteen for whole; five for bar itself
        aspect="landscape", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=True,
    )
    # Write figure to file.
    pplot.write_product_plot_figure(
        figure=figure,
        format="jpg", # svg, jpg, png
        resolution=600,
        name_file=name_figure,
        path_directory=path_directory,
    )
    # Return information.
    return figure


def create_heatmap_figure_for_table_groups(
    factors=None,
    table=None,
    path_directory_parent=None,
    report=None,
):
    """
    Splits rows within table by groups of factor columns.
    Applies procedure to to each group of rows.

    arguments:
        factors (list<str>): names of columns in table by which to split groups
        table (object): Pandas data-frame table of columns for feature variables
            across rows for observation records
        path_directory_parent (str): path to parent directory within which to
            write files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Organize table.
    if False:
        table.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table.drop_duplicates(
            subset=None,
            keep="first",
            inplace=True,
        )
        table.set_index(
            factors,
            append=False,
            drop=True,
            inplace=True
        )
    # Split rows within table by factor columns.
    groups = table.groupby(
        level=factors,
    )
    for name, table_group in groups:
        # Copy information in table.
        table_group = table_group.copy(deep=True)
        # Transpose table.
        table_group = table_group.transpose(copy=True)
        # Organize table.
        if False:
            table_group.reset_index(
                level=None,
                inplace=True,
                drop=False, # do not remove index; move to regular columns
            )

        # Report.
        if report:
            putly.print_terminal_partition(level=4)
            print("Name of group:")
            print(name)
            print("Table for group after split:")
            print(table_group)
            putly.print_terminal_partition(level=4)
        # Complete procedures on each group table after split.
        # For example, calculate summary statistics on each group and then
        # collect within a new summary table.
        if True:
            # Create figure and write to file.
            create_write_figure_heat_map(
                table=table_group,
                name_figure=name,
                path_directory=path_directory_parent,
                report=True,
            )
        pass
    # Return information.
    pass


def control_plot_charts(
    paths=None,
    report=None,
):
    """
    Control plotting of charts.

    arguments:
        paths : (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports


    raises:

    returns:

    """

    # Define instances for iteration.
    instances = [
        {
            "group_analysis": "primary_primary_table_s2_figure_s1", # name of child directory
            "name_table": "table_psychiatry_substance_table_s2_figure_s1", # base name of file
        },
        {
            "group_analysis": "secondary_secondary_table_s3_figure_s2",
            "name_table": "table_thyroid_table_s3_figure_s2",
        },
        {
            "group_analysis": "primary_secondary_table_s4",
            "name_table": "table_psychiatry_substance_thyroid_table_s4",
        },
        {
            "group_analysis": "primary_secondary_figure_1", # must use same parameters for q-values in supplemental table
            "name_table": "table_psychiatry_substance_thyroid_figure_1",
        },

        {
            "group_analysis": "primary_secondary_thyroid_figure_2", # must use same parameters for q-values in supplemental table
            "name_table": "table_psychiatry_substance_thyroid_figure_2",
        },
        {
            "group_analysis": "gonad_psychiatry_substance", # name of child directory
            "name_table": "table_gonad_psychiatry_substance", # base name of file
        },
        {
            "group_analysis": "gonad_psychiatry", # name of child directory
            "name_table": "table_gonad_psychiatry", # base name of file
        },
        {
            "group_analysis": "gonad_substance", # name of child directory
            "name_table": "table_gonad_substance", # base name of file
        },
        #{
        #    "group_analysis": "primary_secondary_sex",
        #    "name_table": "table_psychiatry_substance_sex_hormone",
        #},
    ]
    # Iterate on instances.
    for instance in instances:
        # Initialize child directory for analysis group.
        paths = initialize_directory_group_analysis(
            group_analysis=instance["group_analysis"],
            paths=paths,
            restore=False,
        )
        # Read source information from file.
        table_source = read_organize_source_plot(
            group_analysis=instance["group_analysis"],
            name_table=instance["name_table"],
            path_directory_dock=paths["dock"],
            report=report,
        )
        # Prepare product.
        # Create chart and write to file.
        create_heatmap_figure_for_table_groups(
            factors=["group_analysis"],
            table=table_source,
            path_directory_parent=paths["out_plot_group_analysis"],
            report=report,
        )
        pass
    pass


##########
# 6. Query values in tables


def control_query_genetic_correlation_tables(
    paths=None,
    report=None,
):
    """
    Control procedure to query genetic correlations within tables.

    These queries are fairly low-throughput for writing the Results section of
    the article.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Read source information from file.
    # Organize source information within tables.
    source = read_organize_source_supplemental_tables_for_query(
        paths=paths,
        report=report,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=1)
        print("Queries follow on values of genetic correlation.")
        putly.print_terminal_partition(level=2)
        pass

    # Queries.
    # On 13 March 2024, TCW confirmed the primary and secondary studies in
    # Queries 1-12 (primary-primary).
    # On 13 March 2024, TCW confirmed the primary and secondary studies in
    # Queries 13-15 (secondary-secondary).
    # On 14 March 2024, TCW confirmed the primary and secondary studies in
    # Queries 16-19 (primary-secondary).

    variables_query = [
        #"variants_valid",
        "correlation",
        "p_value_ldsc",
    ]

    # Queries on table of primary-primary genetic correlations.
    # Query 1; Table: primary-secondary
    pextr.query_table_correlation_range_variables(
        table=source["table_rg_one_two"],
        label="Query 1; Table: primary-secondary.",
        studies_primary=[
            "34002096_mullins_2021_bd_all", # BD
            "35396580_trubetskoy_2022_all", # SCZ
        ],
        studies_secondary=[
            "36093044_mathieu_2022_hypothyroidism", # Hypothyroidism
            "34594039_sakaue_2021_eur_hashimoto", # Hashimoto
        ],
        name_primary="study_primary",
        name_secondary="study_secondary",
        exclusions=[],
        variables=variables_query,
        remove_self_pair=True, # whether to exclude self pairs from report
        report=True,
    )

    # Queries on table of primary-primary genetic correlations.
    # Query 2; Table: primary-secondary
    pextr.query_table_correlation_range_variables(
        table=source["table_rg_one_two"],
        label="Query 2; Table: primary-secondary.",
        studies_primary=[
            "34002096_mullins_2021_bd_all", # BD
            "35396580_trubetskoy_2022_all", # SCZ
        ],
        studies_secondary=[
            "32581359_saevarsdottir_2020_thyroid_autoimmunity_dbsnp_rsid", # AITD
        ],
        name_primary="study_primary",
        name_secondary="study_secondary",
        exclusions=[],
        variables=variables_query,
        remove_self_pair=True, # whether to exclude self pairs from report
        report=True,
    )

    # Query 3; Table: primary-secondary
    pextr.query_table_correlation_range_variables(
        table=source["table_rg_one_two"],
        label="Query 3; Table: primary-secondary.",
        studies_primary=[
            "30718901_howard_2019_pgc", # MDD
            "36702997_demontis_2023_adhd", # ADHD
            "31748690_purves_2020_meta", # GAD
            "31594949_nievergelt_2019_europe_all", # PTSD
        ],
        studies_secondary=[
            "36093044_mathieu_2022_hypothyroidism", # Hypothyroidism
            "34594039_sakaue_2021_eur_hashimoto", # Hashimoto
            "32581359_saevarsdottir_2020_thyroid_autoimmunity_dbsnp_rsid", # AITD
        ],
        name_primary="study_primary",
        name_secondary="study_secondary",
        exclusions=[],
        variables=variables_query,
        remove_self_pair=True, # whether to exclude self pairs from report
        report=True,
    )

    # Query 4; Table: primary-secondary
    pextr.query_table_correlation_range_variables(
        table=source["table_rg_one_two"],
        label="Query 4; Table: primary-secondary.",
        studies_primary=[
            "36477530_saunders_2022_tobacco_ever_all", # TOB
            "36477530_saunders_2022_tobacco_all", # TOB-Q
        ],
        studies_secondary=[
            "36093044_mathieu_2022_hypothyroidism", # Hypothyroidism
            "34594039_sakaue_2021_eur_hashimoto", # Hashimoto
            "32581359_saevarsdottir_2020_thyroid_autoimmunity_dbsnp_rsid", # AITD
        ],
        name_primary="study_primary",
        name_secondary="study_secondary",
        exclusions=[],
        variables=variables_query,
        remove_self_pair=True, # whether to exclude self pairs from report
        report=True,
    )

    if False:

        # Query 2; Table: primary-primary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_one"],
            label="Query 2; Table: primary-primary.",
            studies_primary=[
                "30718901_howard_2019_pgc", # MDD
                "36702997_demontis_2023_adhd", # ADHD
                "31748690_purves_2020_meta", # GAD
                "31594949_nievergelt_2019_europe_all", # PTSD
            ],
            studies_secondary=[
                "30718901_howard_2019_pgc", # MDD
                "36702997_demontis_2023_adhd", # ADHD
                "31748690_purves_2020_meta", # GAD
                "31594949_nievergelt_2019_europe_all", # PTSD
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 3; Table: primary-primary
        if True:
            exclusions_query_3 = [
                {
                    "primary": "36477530_saunders_2022_tobacco_all", # TOB-Q
                    "secondary": "32099098_polimanti_2020_eur_opioid_dep_unexposed", # OUD
                },
                {
                    "primary": "32099098_polimanti_2020_eur_opioid_dep_unexposed", # OUD
                    "secondary": "36477530_saunders_2022_tobacco_all", # TOB-Q
                },
            ]
        else:
            exclusions_query_3 = []
            pass
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_one"],
            label="Query 3; Table: primary-primary.",
            studies_primary=[
                "30482948_walters_2018_eur_all", # AUD
                "36477530_saunders_2022_alcohol_all", # ALC-Q
                "36477530_saunders_2022_tobacco_ever_all", # TOB
                "36477530_saunders_2022_tobacco_all", # TOB-Q
                "33096046_johnson_2020_eur_all", # CUD
                "32099098_polimanti_2020_eur_opioid_dep_unexposed", # OUD
            ],
            studies_secondary=[
                "30482948_walters_2018_eur_all", # AUD
                "36477530_saunders_2022_alcohol_all", # ALC-Q
                "36477530_saunders_2022_tobacco_ever_all", # TOB
                "36477530_saunders_2022_tobacco_all", # TOB-Q
                "33096046_johnson_2020_eur_all", # CUD
                "32099098_polimanti_2020_eur_opioid_dep_unexposed", # OUD
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=exclusions_query_3,
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 4; Table: primary-primary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_one"],
            label="Query 4; Table: primary-primary.",
            studies_primary=[
                "34002096_mullins_2021_bd_all", # BD
                "35396580_trubetskoy_2022_all", # SCZ
            ],
            studies_secondary=[
                "30482948_walters_2018_eur_all", # AUD
                "36477530_saunders_2022_alcohol_all", # ALC-Q
                "36477530_saunders_2022_tobacco_ever_all", # TOB
                "36477530_saunders_2022_tobacco_all", # TOB-Q
                "33096046_johnson_2020_eur_all", # CUD
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 5; Table: primary-primary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_one"],
            label="Query 5; Table: primary-primary.",
            studies_primary=[
                "30718901_howard_2019_pgc", # MDD
                "36702997_demontis_2023_adhd", # ADHD
                "31748690_purves_2020_meta", # GAD
                "31594949_nievergelt_2019_europe_all", # PTSD
            ],
            studies_secondary=[
                "30482948_walters_2018_eur_all", # AUD
                "36477530_saunders_2022_tobacco_ever_all", # TOB
                "36477530_saunders_2022_tobacco_all", # TOB-Q
                "33096046_johnson_2020_eur_all", # CUD
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 6; Table: primary-primary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_one"],
            label="Query 6; Table: primary-primary.",
            studies_primary=[
                "35396580_trubetskoy_2022_all", # SCZ
                "30718901_howard_2019_pgc", # MDD
                "36702997_demontis_2023_adhd", # ADHD
                "31748690_purves_2020_meta", # GAD
            ],
            studies_secondary=[
                "32099098_polimanti_2020_eur_opioid_dep_unexposed", # OUD
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 7; Table: primary-primary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_one"],
            label="Query 7; Table: primary-primary.",
            studies_primary=[
                "28761083_arnold_2018", # OCD
            ],
            studies_secondary=[
                "31308545_watson_2019", # ANN
                "30818990_yu_2019", # TRS
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 8; Table: primary-primary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_one"],
            label="Query 8; Table: primary-primary.",
            studies_primary=[
                "34002096_mullins_2021_bd_all", # BD
                "35396580_trubetskoy_2022_all", # SCZ
                "30718901_howard_2019_pgc", # MDD
                "31748690_purves_2020_meta", # GAD
            ],
            studies_secondary=[
                "32747698_matoba_2020_europe", # ASD
                "28761083_arnold_2018", # OCD
                "31308545_watson_2019", # ANN
                "30818990_yu_2019", # TRS
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 9; Table: primary-primary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_one"],
            label="Query 9; Table: primary-primary.",
            studies_primary=[
                "36702997_demontis_2023_adhd", # ADHD
            ],
            studies_secondary=[
                "32747698_matoba_2020_europe", # ASD
                "31308545_watson_2019", # ANN
                "30818990_yu_2019", # TRS
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 10; Table: primary-primary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_one"],
            label="Query 10; Table: primary-primary.",
            studies_primary=[
                "31594949_nievergelt_2019_europe_all", # PTSD
            ],
            studies_secondary=[
                "32747698_matoba_2020_europe", # ASD
                "31308545_watson_2019", # ANN
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 11; Table: primary-primary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_one"],
            label="Query 11; Table: primary-primary.",
            studies_primary=[
                "32747698_matoba_2020_europe", # ASD
                "30818990_yu_2019", # TRS
            ],
            studies_secondary=[
                "36477530_saunders_2022_tobacco_all", # TOB-Q
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 12; Table: primary-primary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_one"],
            label="Query 12; Table: primary-primary.",
            studies_primary=[
                "28761083_arnold_2018", # OCD
                "31308545_watson_2019", # ANN
            ],
            studies_secondary=[
                "36477530_saunders_2022_tobacco_all", # TOB-Q
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )

        # Queries on table of secondary-secondary genetic correlations.
        # Query 13; Table: secondary-secondary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_two_two"],
            label="Query 13; Table: secondary-secondary.",
            studies_primary=[
                "32581359_saevarsdottir_2020_thyroid_autoimmunity_dbsnp_rsid", # AITD
                "36093044_mathieu_2022_hypothyroidism", # Hypothyroidism
                "34594039_sakaue_2021_eur_hashimoto", # Hashimoto
            ],
            studies_secondary=[
                "32581359_saevarsdottir_2020_thyroid_autoimmunity_dbsnp_rsid", # AITD
                "36093044_mathieu_2022_hypothyroidism", # Hypothyroidism
                "34594039_sakaue_2021_eur_hashimoto", # Hashimoto
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 14; Table: secondary-secondary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_two_two"],
            label="Query 14; Table: secondary-secondary.",
            studies_primary=[
                "32581359_saevarsdottir_2020_thyroid_autoimmunity_dbsnp_rsid", # AITD
                "36093044_mathieu_2022_hypothyroidism", # Hypothyroidism
                "34594039_sakaue_2021_eur_hashimoto", # Hashimoto
            ],
            studies_secondary=[
                "34594039_sakaue_2021_eur_hyperthyroidism", # Hyperthyroidism
                "34594039_sakaue_2021_eur_graves", # Graves
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 15; Table: secondary-secondary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_two_two"],
            label="Query 15; Table: secondary-secondary.",
            studies_primary=[
                "32581359_saevarsdottir_2020_thyroid_autoimmunity_dbsnp_rsid", # AITD
                "36093044_mathieu_2022_hypothyroidism", # Hypothyroidism
                "34594039_sakaue_2021_eur_hashimoto", # Hashimoto
                "34594039_sakaue_2021_eur_hyperthyroidism", # Hyperthyroidism
                "34594039_sakaue_2021_eur_graves", # Graves
            ],
            studies_secondary=[
                "24586183_medici_2014_thyroid_peroxidase_reactivity", # Anti-TPO-React
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )

        # Queries on table of primary-secondary genetic correlations.
        # Query 16; Table: primary-secondary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_two"],
            label="Query 16; Table: primary-secondary.",
            studies_primary=[
                "34002096_mullins_2021_bd_all", # BD
                "35396580_trubetskoy_2022_all", # SCZ
            ],
            studies_secondary=[
                "36093044_mathieu_2022_hypothyroidism", # Hypothyroidism
                "34594039_sakaue_2021_eur_hashimoto", # Hashimoto
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 17; Table: primary-secondary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_two"],
            label="Query 17; Table: primary-secondary.",
            studies_primary=[
                "34002096_mullins_2021_bd_all", # BD
                "35396580_trubetskoy_2022_all", # SCZ
            ],
            studies_secondary=[
                "32581359_saevarsdottir_2020_thyroid_autoimmunity_dbsnp_rsid", # AITD
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 18; Table: primary-secondary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_two"],
            label="Query 18; Table: primary-secondary.",
            studies_primary=[
                "30718901_howard_2019_pgc", # MDD
                "36702997_demontis_2023_adhd", # ADHD
                "31748690_purves_2020_meta", # GAD
                "31594949_nievergelt_2019_europe_all", # PTSD
            ],
            studies_secondary=[
                "32581359_saevarsdottir_2020_thyroid_autoimmunity_dbsnp_rsid", # AITD
                "36093044_mathieu_2022_hypothyroidism", # Hypothyroidism
                "34594039_sakaue_2021_eur_hashimoto", # Hashimoto
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
        # Query 19; Table: primary-secondary
        pextr.query_table_correlation_range_variables(
            table=source["table_rg_one_two"],
            label="Query 19; Table: primary-secondary.",
            studies_primary=[
                "36477530_saunders_2022_tobacco_ever_all", # TOB
                "36477530_saunders_2022_tobacco_all", # TOB-Q
            ],
            studies_secondary=[
                "32581359_saevarsdottir_2020_thyroid_autoimmunity_dbsnp_rsid", # AITD
                "36093044_mathieu_2022_hypothyroidism", # Hypothyroidism
                "34594039_sakaue_2021_eur_hashimoto", # Hashimoto
            ],
            name_primary="study_primary",
            name_secondary="study_secondary",
            exclusions=[],
            variables=variables_query,
            remove_self_pair=True, # whether to exclude self pairs from report
            report=True,
        )
    pass





################################################
################################################
################################################
# TODO: TCW; 6 June 2024
# Everything below this point (except the "execute_procedure" main function) needs updates.




##########
# Organize tables for node-link diagram

# TODO:

# 1. iterate across files in the child directory for primary-secondary correlations
# for each table, create a new column with the file name (minus ".tsv") as the new column "study_primary".
# Use Pandas concat to combine all of these tables into one.
# reduce the columns to include "study_primary", "study_secondary", "correlation", "error", and "p_value"

# 2. read and organize the primary-primary and secondary-secondary tables
# split the identifiers into "study_primary" and "study_secondary"
# organize the columns to match the primary-secondary table from step 1

# 3. concatenate all "link" tables

# 4. filter the links

# filter node attribute table by inclusion == 1
# extract list of primary studies (category == "psychiatry" or "substance_use") and secondary studies

# collect information about all links
# iterate over the primary-secondary tables
#    primary study = name of file - ".tsv"
#    secondary study = name in the table
# combine the primary-secondary tables with the primary-primary and primary-secondary tables
# filter links in the combined link table by the primary and secondary studies for inclusion from node attribute table



###############################################################################
# Procedure


def execute_procedure(
    path_directory_dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_directory_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:

    """

    ##########
    # Parameters.
    project="psychiatry_biomarkers"
    routine="genetic_correlation"
    procedure="thyroid_organization"
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print(
            "module: psychiatry_biomarkers.genetic_correlation." +
            "thyroid_organization.py"
        )
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("project: " + str(project))
        print("routine: " + str(routine))
        print("procedure: " + str(procedure))
        putly.print_terminal_partition(level=5)
        pass

    if False:

        ##########
        # 1. Initialize directories for read of source and write of product files.
        paths = initialize_directories(
            project=project,
            routine=routine,
            procedure=procedure,
            path_directory_dock=path_directory_dock,
            restore=True,
            report=report,
        )

        ##########
        # 2. Organize SNP heritabilities within tables for supplement.
        table_h2 = control_read_organize_snp_heritability_table_supplement(
            paths=paths,
            report=True,
        )

        ##########
        # 3. Read and assemble all genetic correlations within main table for
        # subsequent organization.
        table_rg_total = control_assemble_genetic_correlations(
            paths=paths,
            report=True,
        )

        ##########
        # 4. Organize genetic correlations within tables for reporting as text or
        # plot.
        control_prepare_genetic_correlation_tables_supplement_plot(
            table_rg=table_rg_total,
            paths=paths,
            report=True,
        )

        ##########
        # 5. Organize tables of nodes and links for network representation of
        # genetic correlations.
        control_prepare_genetic_correlation_network_nodes_links(
            paths=paths,
            report=True,
        )

        pass


    ##########
    # 6. Query genetic correlations for convenient reporting within the text
    # of the article's Results.
    if True:
        paths = initialize_directories(
            project=project,
            routine=routine,
            procedure=procedure,
            path_directory_dock=path_directory_dock,
            restore=False,
            report=report,
        )
        control_query_genetic_correlation_tables(
            paths=paths,
            report=report,
        )

    ##########
    # 7. Plot charts to represent genetic correlations.
    if False:
        control_plot_charts(
            paths=paths,
            report=True,
        )

    pass





###############################################################################
# End
