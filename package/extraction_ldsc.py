
"""
Organize procedural code for extraction of information from raw text report logs
that the LDSC tool creates for estimates of SNP heritability (h2) and
genetic correlation (rg).

This module makes use of the 'extraction' module within the 'partner' package.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Monroe, North Carolina 28110
    United States of America

License:

    This file is part of project 'psychiatry_biomarkers'
    (https://github.com/tcameronwaller/psychiatry_biomarkers/).

    Project 'psychiatry_biomarkers' supports data analysis for publications
    about biomarkers in psychiatric disorders.
    Copyright (C) 2024 Thomas Cameron Waller

    Project 'psychiatry_biomarkers' is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    Project 'psychiatry_biomarkers' is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with project 'psychiatry_biomarkers'.
    If not, see <http://www.gnu.org/licenses/>.
"""

###############################################################################
# Notes

# Primary studies include a curated collection of summary statistics from prior
# GWAS on psychiatric disorders and substance-use disorders.

# Secondary studies include a curated collection of summary statistics from
# prior GWAS on thyroid disorders, thyroid biomarkers, sex hormones,
# other hormones, and other biomarkers.

###############################################################################
# Installation and importation

# Import modules from specific path without having to install a general package
# I would have to figure out how to pass a path variable...
# https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path


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
import pandas
pandas.options.mode.chained_assignment = None # default = "warn"

import scipy.stats
import scipy.linalg
import statsmodels.multivariate.pca

# Custom
import partner.utility as putly # this import path for subpackage
import partner.extraction as pextr

###############################################################################
# Functionality


##########
# Initialization


def initialize_directories(
    restore=None,
    path_dock=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        restore (bool): whether to remove previous versions of data
        path_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    paths["dock"] = path_dock

    # Source paths.
    paths["parameters"] = os.path.join(
        path_dock, "parameters", "psychiatric_metabolism",
    )
    paths["heritability"] = os.path.join(
        path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "5_gwas_heritability_ldsc",
    )
    paths["heritability_no_liability"] = os.path.join(
        path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "5_gwas_heritability_ldsc_no_liability",
    )


    paths["correlation_primary"] = os.path.join(
        path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "6_gwas_correlation_ldsc_primary",
    )
    paths["correlation_secondary"] = os.path.join(
        path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "6_gwas_correlation_ldsc_secondary",
    )

    paths["correlation_primary_secondary"] = os.path.join(
        path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "6_gwas_correlation_ldsc_primary_secondary",
    )



    paths["extraction_h2"] = os.path.join(
        path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "extraction_5_gwas_heritability_ldsc",
    )
    paths["extraction_rg_one_one"] = os.path.join(
        path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "extraction_6_gwas_correlation_ldsc_primary",
    )
    paths["extraction_rg_two_two_thyroid"] = os.path.join(
        path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "extraction_6_gwas_correlation_ldsc_secondary_thyroid",
    )
    paths["extraction_rg_two_two_all"] = os.path.join(
        path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "extraction_6_gwas_correlation_ldsc_secondary_all",
    )
    paths["extraction_rg_one_two"] = os.path.join(
        path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "extraction_6_gwas_correlation_ldsc_primary_secondary",
    )


    # Remove previous files to avoid version or batch confusion.
    if restore:
        putly.remove_directory(path=paths["extraction_h2"])
        putly.remove_directory(path=paths["extraction_rg_one_one"])
        putly.remove_directory(path=paths["extraction_rg_two_two_thyroid"])
        putly.remove_directory(path=paths["extraction_rg_two_two_all"])
        putly.remove_directory(path=paths["extraction_rg_one_two"])
    # Initialize directories.
    putly.create_directories(
        path=paths["extraction_h2"]
    )
    putly.create_directories(
        path=paths["extraction_rg_one_two"]
    )
    putly.create_directories(
        path=paths["extraction_rg_one_one"]
    )
    putly.create_directories(
        path=paths["extraction_rg_two_two_thyroid"]
    )
    # Return information.
    return paths



###############################################################################
# Procedure



def execute_procedure(
    path_dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Report version.
    putly.print_terminal_partition(level=1)
    print(path_dock)
    print("version check: 1")
    # Pause procedure.
    time.sleep(5.0)

    # Initialize directories.
    paths = initialize_directories(
        restore=True,
        path_dock=path_dock,
    )
    #paths["correlation"]
    #paths["correlation_extraction"]
    #paths["heritability_biomarkers"]
    #paths["heritability_biomarkers_extraction"]
    #paths["heritability_disorders"]
    #paths["heritability_disorders_extraction"]

    ##########
    # Manage extraction of information about SNP heritability.

    if True:
        # Optionally read list of indices by which to sort rows in each
        # extraction table of genetic correlations.
        path_file_list_sort = os.path.join(
            paths["parameters"], "list_sort_heritability_psychiatry_substance_thyroid.txt",
        )
        list_sort = putly.read_file_text_list(
            delimiter="\n",
            path_file=path_file_list_sort,
        )
        #print("Length of list of indices: " + str(len(list_sort)))
        #print(list_sort)
        # Collect information.
        pail_write_heritability = dict()
        # Extract information from reports of analyses in LDSC.
        table_heritability = pextr.read_extract_from_all_ldsc_files_in_directory(
            path_directory=paths["heritability"],
            file_name_pattern=".log",
            file_name_pattern_not=".....",
            analysis="heritability",
            report=True,
        )
        if True:
            table_heritability = (
                putly.sort_table_rows_by_list_indices(
                    table=table_heritability,
                    list_sort=list_sort,
                    name_column="sort_rows_temporary",
                    report=True,
            ))
        pail_write_heritability["table_heritability"] = table_heritability
        # Extract information from reports of analyses in LDSC.
        table_heritability_no_liability = pextr.read_extract_from_all_ldsc_files_in_directory(
            path_directory=paths["heritability_no_liability"],
            file_name_pattern=".log",
            file_name_pattern_not=".....",
            analysis="heritability",
            report=True,
        )
        if True:
            table_heritability_no_liability = (
                putly.sort_table_rows_by_list_indices(
                    table=table_heritability_no_liability,
                    list_sort=list_sort,
                    name_column="sort_rows_temporary",
                    report=True,
            ))
        pail_write_heritability["table_heritability_no_liability"] = table_heritability_no_liability
        # Write product information to file.
        putly.write_product_tables(
            pail_write=pail_write_heritability,
            path_directory=paths["extraction_h2"],
        )



    ##########
    # Manage extraction of information about genetic correlations between
    # pairs of primary and secondary traits.

    if True:
        # Optionally read list of indices by which to sort rows in each
        # extraction table of genetic correlations.
        path_file_list_sort = os.path.join(
            paths["parameters"], "list_sort_correlation_thyroid_sex_biomarkers_full.txt",
        )
        list_sort = putly.read_file_text_list(
            delimiter="\n",
            path_file=path_file_list_sort,
        )
        # Collect information.
        pail_write_correlation = dict()
        # Extract names of child directories within parent directory.
        names_directories = putly.extract_subdirectory_names(
            path=paths["correlation_primary_secondary"]
        )
        names_directories_ldsc = list(filter(
            lambda name: (name != "batch"),
            names_directories
        ))
        print("--------------------")
        print(names_directories_ldsc)
        print("--------------------")
        # Write each table to file.
        for name_directory in names_directories_ldsc:
            path_directory = os.path.join(
                paths["correlation_primary_secondary"], name_directory,
            )
            table_correlation = pextr.read_extract_from_all_ldsc_files_in_directory(
                path_directory=path_directory,
                file_name_pattern=".log",
                file_name_pattern_not=".....",
                analysis="correlation",
                report=True,
            )
            if True:
                table_correlation = (
                    putly.sort_table_rows_by_list_indices(
                        table=table_correlation,
                        list_sort=list_sort,
                        name_column="sort_rows_temporary",
                        report=True,
                ))
            pail_write_correlation[str("table_" + name_directory)] = table_correlation
        # Write product information to file.
        putly.write_product_tables(
            pail_write=pail_write_correlation,
            path_directory=paths["extraction_rg_one_two"],
        )


    ##########
    # Manage extraction of information about genetic correlations between
    # pairs of primary and primary traits.

    if True:
        # Optionally read list of indices by which to sort rows in each
        # extraction table of genetic correlations.
        path_file_list_sort = os.path.join(
            paths["parameters"], "list_sort_correlation_psychiatry_substance.txt",
        )
        list_sort = putly.read_file_text_list(
            delimiter="\n",
            path_file=path_file_list_sort,
        )
        # Collect information.
        pail_write_correlation = dict()
        # Extract names of child directories within parent directory.
        names_directories = putly.extract_subdirectory_names(
            path=paths["correlation_primary"]
        )
        names_directories_ldsc = list(filter(
            lambda name: (name != "batch"),
            names_directories
        ))
        print("--------------------")
        print(names_directories_ldsc)
        print("--------------------")
        # Write each table to file.
        for name_directory in names_directories_ldsc:
            path_directory = os.path.join(
                paths["correlation_primary"], name_directory,
            )
            table_correlation = pextr.read_extract_from_all_ldsc_files_in_directory(
                path_directory=path_directory,
                file_name_pattern=".log",
                file_name_pattern_not=".....",
                analysis="correlation",
                report=True,
            )
            if True:
                table_correlation = (
                    putly.sort_table_rows_by_list_indices(
                        table=table_correlation,
                        list_sort=list_sort,
                        name_column="sort_rows_temporary",
                        report=True,
                ))
            pail_write_correlation[str("table_" + name_directory)] = table_correlation
        # Write product information to file.
        putly.write_product_tables(
            pail_write=pail_write_correlation,
            path_directory=paths["extraction_rg_one_one"],
        )


    ##########
    # Manage extraction of information about genetic correlations between
    # pairs of secondary and secondary traits.

    if True:
        # Optionally read list of indices by which to sort rows in each
        # extraction table of genetic correlations.
        path_file_list_sort = os.path.join(
            paths["parameters"], "list_sort_correlation_thyroid_disorders_biomarkers.txt",
        )
        if True:
            list_sort = putly.read_file_text_list(
                delimiter="\n",
                path_file=path_file_list_sort,
            )
        # Collect information.
        pail_write_correlation = dict()
        # Extract names of child directories within parent directory.
        names_directories = putly.extract_subdirectory_names(
            path=paths["correlation_secondary"]
        )
        names_directories_ldsc = list(filter(
            lambda name: (name != "batch"),
            names_directories
        ))
        print("--------------------")
        print(names_directories_ldsc)
        print("--------------------")
        # Write each table to file.
        for name_directory in names_directories_ldsc:
            path_directory = os.path.join(
                paths["correlation_secondary"], name_directory,
            )
            table_correlation = pextr.read_extract_from_all_ldsc_files_in_directory(
                path_directory=path_directory,
                file_name_pattern=".log",
                file_name_pattern_not=".....",
                analysis="correlation",
                report=True,
            )
            if True:
                table_correlation = (
                    putly.sort_table_rows_by_list_indices(
                        table=table_correlation,
                        list_sort=list_sort,
                        name_column="sort_rows_temporary",
                        report=True,
                ))
            pail_write_correlation[str("table_" + name_directory)] = table_correlation
        # Write product information to file.
        putly.write_product_tables(
            pail_write=pail_write_correlation,
            path_directory=paths["extraction_rg_two_two_thyroid"],
        )



    pass





#
