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
#import partner.plot as pplot
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
    paths["group_analysis"] = os.path.join(
        paths["out_data"], group_analysis,
    )
    # Remove previous files to avoid version or batch confusion.
    if restore:
        putly.remove_directory(path=paths["group_analysis"])
    # Initialize directories.
    putly.create_directories(
        path=paths["group_analysis"]
    )
    # Return information.
    return paths


##########
# 2. Organize SNP heritabilities in table for supplement


##########
# 3. Read and assemble genetic correlations


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
# 3. Read parameters about studies


def read_organize_source_parameter_genetic_correlations(
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
        paths["in_parameters_private"],
    )
    # Define paths to files.
    path_file_table_studies = os.path.join(
        path_directory_parent,
        "table_studies_attributes_filter_sort.tsv",
    )
    # Collect information.
    pail = dict()
    # Read information from file.
    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["inclusion_thyroid_table"] = "float"
    types_columns["inclusion_thyroid_figure"] = "float"
    types_columns["inclusion_sex_table"] = "float"
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
    pail["table_studies"] = pandas.read_csv(
        path_file_table_studies,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )
    # Return information.
    return pail





##########
# Organize genetic correlations in tables for supplement


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


# TODO: TCW; 9 August 2024
# It seems to be necessary to choose between "p_value" or "q_value" when
# organizing the tables for the plots, especially the symmetrical plot.

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
        path_directory_parent=paths["group_analysis"],
        reset_index=False,
        write_index=True,
        type="text",
    )
    putly.write_tables_to_file_in_child_directories(
        pail_write=pail_write_directories_pickle,
        path_directory_parent=paths["group_analysis"],
        reset_index=False,
        write_index=True,
        type="pickle",
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
    source = read_organize_source_parameter_genetic_correlations(
        paths=paths,
        report=report,
    )
    # Collect parameters common across all instances.
    parameters = dict()
    parameters["table_studies"] = source["table_studies"]
    parameters["table_rg"] = table_rg
    parameters["paths"] = paths
    parameters["report"] = report

    # Collect parameters specific to each instance.
    instances = [
        {
            "group_analysis": "primary_primary", # name of child directory
            "name_table": "table_psychiatry_substance", # base name of file
            "inclusion_table": "inclusion_thyroid_table",
            "inclusion_figure": "inclusion_thyroid_figure",
            "groups_primary": ["psychiatry","substance",],
            "groups_secondary": ["psychiatry","substance",],
            "symmetry": True,
        },
        {
            "group_analysis": "secondary_secondary_thyroid",
            "name_table": "table_thyroid",
            "inclusion_table": "inclusion_thyroid_table",
            "inclusion_figure": "inclusion_thyroid_figure",
            "groups_primary": ["thyroid",],
            "groups_secondary": ["thyroid",],
            "symmetry": True,
        },
        {
            "group_analysis": "primary_secondary_thyroid",
            "name_table": "table_psychiatry_substance_thyroid",
            "inclusion_table": "inclusion_thyroid_table",
            "inclusion_figure": "inclusion_thyroid_figure",
            "groups_primary": ["psychiatry","substance",],
            "groups_secondary": ["thyroid",],
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
    if False:
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
    table_test.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Collect information.
    # Collections of files.
    pail_write_files = dict()
    pail_write_files[str("test_table_check")] = table_test
    ##########
    # Write product information to file.
    putly.write_tables_to_file(
        pail_write=pail_write_files,
        path_directory_parent=paths["out_test"],
        reset_index=False,
        write_index=False,
        type="text",
    )
    ##################################################
    pass


##########
# Network


def control_prepare_genetic_correlation_network_links(
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
    source = read_organize_source_parameter_genetic_correlations(
        path_directory_dock=paths["dock"],
        report=report,
    )
    table_studies = source["table_studies"]

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
        path_directory_parent=paths["group_analysis"],
        type="text",
    )
    pass





################################################
################################################
################################################
# TODO: TCW; 6 June 2024
# Everything below this point (except the "execute_procedure" main function) needs updates.


##########
# Query tables of genetic correlations


def read_source_tables_genetic_correlation(
    path_directory_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        path_directory_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables with entry names
            (keys) derived from original names of files

    """

    # Define path to parent directory.
    path_directory_parent = os.path.join(
        path_directory_dock, "parameters", "scratch",
        "tables_genetic_correlation",
    )
    # Define paths to files.
    path_file_table_one_one = os.path.join(
        path_directory_parent, "table_rg_primary_primary.tsv",
    )
    path_file_table_two_two = os.path.join(
        path_directory_parent, "table_rg_secondary_secondary.tsv",
    )
    path_file_table_one_two = os.path.join(
        path_directory_parent, "table_rg_primary_secondary.tsv",
    )

    # Collect tables.
    pail = dict()

    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["study_primary"] = "string"
    types_columns["study_secondary"] = "string"
    types_columns["abbreviation_primary"] = "string"
    types_columns["abbreviation_secondary"] = "string"
    types_columns["trait_primary"] = "string"
    types_columns["trait_secondary"] = "string"
    types_columns["sex_primary"] = "string"
    types_columns["sex_secondary"] = "string"
    types_columns["variants_valid"] = "float"
    types_columns["correlation"] = "float"
    types_columns["correlation_error"] = "float"
    types_columns["z_score"] = "float" # <-- will be obsolete after update extraction
    types_columns["p_not_zero"] = "float" # <-- will be obsolete after update extraction
    types_columns["q_not_zero"] = "float" # <-- will be obsolete after update extraction
    #types_columns["z_statistic_ldsc"] = "float"
    #types_columns["p_value_ldsc"] = "float"
    #types_columns["q_value_ldsc"] = "float"
    #types_columns["z_statistic_not_zero"] = "float"
    #types_columns["p_value_not_zero"] = "float"
    #types_columns["q_value_not_zero"] = "float"
    #types_columns["z_statistic_less_one"] = "float"
    #types_columns["p_value_less_one"] = "float"
    #types_columns["q_value_less_one"] = "float"
    types_columns["correlation_ci95_low"] = "float"
    types_columns["correlation_ci95_high"] = "float"
    types_columns["correlation_ci99_low"] = "float"
    types_columns["correlation_ci99_high"] = "float"
    types_columns["summary_correlation_error"] = "string"
    types_columns["summary_correlation_ci95"] = "string"
    types_columns["summary_correlation_ci99"] = "string"
    # Read information from file.
    pail["table_rg_one_one"] = pandas.read_csv(
        path_file_table_one_one,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )
    pail["table_rg_two_two"] = pandas.read_csv(
        path_file_table_two_two,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )
    pail["table_rg_one_two"] = pandas.read_csv(
        path_file_table_one_two,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
    )
    # Return information.
    return pail


def control_query_genetic_correlation_tables(
    paths=None,
):
    """
    Control procedure to query genetic correlations within tables.

    These queries are fairly low-throughput for writing the Results section of
    the article.

    arguments:
        paths : (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Read source information from file.
    # Organize source information within tables.
    source = read_source_tables_genetic_correlation(
        path_directory_dock=paths["dock"],
        report=True,
    )

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
        "correlation_ci95_low",
        "correlation_ci95_high",
        "p_not_zero", # <-- obsolete after extraction update
        "q_not_zero", # <-- obsolete after extraction update
        #"p_value_not_zero",
        #"q_value_not_zero",
    ]

    # Queries on table of primary-primary genetic correlations.
    # Query 1; Table: primary-primary
    pextr.query_table_correlation_range_variables(
        table=source["table_rg_one_one"],
        label="Query 1; Table: primary-primary.",
        studies_primary=[
            "34002096_mullins_2021_bd_all", # BD
            "35396580_trubetskoy_2022_all", # SCZ
            "30718901_howard_2019_pgc", # MDD
        ],
        studies_secondary=[
            "34002096_mullins_2021_bd_all", # BD
            "35396580_trubetskoy_2022_all", # SCZ
            "30718901_howard_2019_pgc", # MDD
        ],
        name_primary="study_primary",
        name_secondary="study_secondary",
        exclusions=[],
        variables=variables_query,
        remove_self_pair=True, # whether to exclude self pairs from report
        report=True,
    )
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




    if False:
        # Collect tables.
        pail_write = dict()
        #pail_write["tables_genetic_correlation"] = dict()
        pail_write_sub = dict()
        pail_write_sub["table_rg_primary_primary"] = table_rg_one_one_clean
        pail_write_sub["table_rg_secondary_secondary"] = table_rg_two_two_clean
        pail_write_sub["table_rg_primary_secondary"] = table_rg_one_two_clean
        pail_write["tables_genetic_correlation"] = pail_write_sub
        # Write product information to file.
        putly.write_tables_to_file_in_child_directories(
            pail_write=pail_write,
            path_directory_parent=paths["wrangler"],
        )
    pass



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


# TODO: needs update
def read_source_genetic_correlations_node_link_diagram(
    path_directory_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        path_directory_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables with entry names
            (keys) derived from original names of files

    """

    # Define path to parent directory.
    path_dir_parent = os.path.join(
        path_directory_dock, "parameters", "scratch",
        "table_chart_correlation_network_2024-02-22",
    )
    # Define path to child directory.
    path_directory_primary_secondary = os.path.join(
        path_dir_parent, "6_gwas_correlation_ldsc_primary_secondary_extraction",
    )
    # Define paths to files.
    path_file_table_studies = os.path.join(
        path_dir_parent, "table_studies_attributes.tsv",
    )
    path_file_table_primary_primary = os.path.join(
        path_dir_parent, "table_neuropsychiatry_substance_disorders.tsv",
    )
    path_file_table_secondary_secondary = os.path.join(
        path_dir_parent, "table_thyroid_disorders_biomarkers.tsv",
    )

    # Specify variable types in tables.
    types_columns = dict()
    types_columns["path_directory"] = "string"
    types_columns["name_file"] = "string"
    types_columns["summary_correlation_error"] = "string"
    types_columns["summary_correlation_ci95"] = "string"
    types_columns["summary_correlation_ci99"] = "string"
    types_columns["variants"] = "float"
    types_columns["variants_valid"] = "float"
    types_columns["covariance"] = "float"
    types_columns["covariance_error"] = "float"
    types_columns["z_score"] = "float"
    types_columns["correlation"] = "float"
    types_columns["correlation_error"] = "float"
    types_columns["p_not_zero"] = "float"

    # Collect tables.
    pail = dict()

    # Read information from file.
    pail["table_studies"] = pandas.read_csv(
        path_file_table_studies,
        sep="\t",
        header=0,
        na_values=["nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",],
    )

    # Read and organize information from file.

    # Correlations between primary-primary studies and secondary-secondary
    # studies.
    pail["primary_primary"] = (
        pextr.read_organize_table_ldsc_correlation_single(
            path_file_table=path_file_table_primary_primary,
            types_columns=types_columns,
            report=True,
    ))
    pail["secondary_secondary"] = (
        pextr.read_organize_table_ldsc_correlation_single(
            path_file_table=path_file_table_secondary_secondary,
            types_columns=types_columns,
            report=True,
    ))

    # Correlations between primary-secondary studies.
    pail["primary_secondary"] = (
        pextr.read_organize_table_ldsc_correlation_multiple(
            path_directory_parent=path_directory_primary_secondary,
            types_columns=types_columns,
            report=True,
    ))

    # Return information.
    return pail


def combine_tables_correlation(
    table_1=None,
    table_2=None,
    table_3=None,
    report=None,
):
    """
    Combines information from multiple tables.

    arguments:
        table_1 (object): Pandas data-frame table
        table_2 (object): Pandas data-frame table
        table_3 (object): Pandas data-frame table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in first table.
    table = table_1.copy(deep=True)
    # Concatenate tables.
    table = pandas.concat(
        [table, table_2,],
        axis="index",
        join="outer",
        ignore_index=True,
        copy=True,
    )
    table = pandas.concat(
        [table, table_3,],
        axis="index",
        join="outer",
        ignore_index=True,
        copy=True,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Concatenation of tables:")
        print(table)
        print("Column labels:")
        labels_columns = table.columns.to_list()
        print(labels_columns)
        putly.print_terminal_partition(level=4)
        count_table_rows = (table.shape[0])
        count_table_columns = (table.shape[1])
        print("Table rows: " + str(count_table_rows))
        print("Table columns: " + str(count_table_columns))
        putly.print_terminal_partition(level=4)
    # Return information.
    return table


# TODO: needs update
def control_prepare_genetic_correlation_network(
    paths=None,
):
    """
    Control procedure to test extraction of information from LDSC report text
    logs for SNP heritability and genetic correlation.

    arguments:
        paths : (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Genome-Wide Association Studies (GWAS) as Nodes and Genetic Correlations
    # as Links between them.

    # Read source information from file.
    # Organize source information within tables.
    source = read_source_genetic_correlations_node_link_diagram(
        path_directory_dock=paths["dock"],
        report=True,
    )
    # Combine tables of genetic correlations.
    table = combine_tables_correlation(
        table_1=source["table_correlations_1_1"],
        table_2=source["table_correlations_2_2"],
        table_3=source["table_correlations_1_2"],
        report=True,
    )
    # Filter genetic correlations.
    table_studies_inclusion = source["table_studies"].loc[
        (source["table_studies"]["inclusion"] == 1), :
    ]
    studies_inclusion = table_studies_inclusion["identifier"].to_list()
    order_columns = [
        #"index",
        "study_primary",
        "study_secondary",
        "variants",
        "variants_valid",
        "correlation",
        "correlation_absolute",
        "correlation_error",
        "p_not_zero",
        "q_not_zero",
    ]
    table_filter = pextr.filter_table_ldsc_correlation_studies(
        table=table,
        studies_keep=studies_inclusion,
        threshold_q=0.05,
        order_columns=order_columns,
        report=True,
    )
    # Collect tables.
    pail_write = dict()
    pail_write["tables"] = dict()
    pail_write["tables"]["table_study_correlations_q-value_0-005"] = table_filter
    # Write product information to file.
    putly.write_tables_to_file_in_child_directories(
        pail_write=pail_write,
        path_directory_parent=paths["wrangler"],
    )
    pass



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
    # TODO: TCW; 27 May 2024
    # TODO: implement this... only 1 SNP h2 table for supplement.
    # TODO: need to filter and organize

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



    if False:

        ##########
        # Organize genetic correlations within tables for network.
        if False:
            control_prepare_genetic_correlation_network_links(
                table_rg=table_rg_total,
                paths=paths,
                report=True,
            )

        if False:
            control_prepare_genetic_correlation_network(
                paths=paths,
            )

        ##########
        # Query genetic correlations within tables for article's Results.
        if False:
            control_query_genetic_correlation_tables(
                paths=paths,
            )

    pass





#
