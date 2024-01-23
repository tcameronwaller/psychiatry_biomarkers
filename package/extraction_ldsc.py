
"""

This module contains functions for general organization of the raw extraction
data from the UK Biobank accessions.

Combine information from multiple UK Biobank accessions.
Integrate identifier matchings.
Exclude individual cases (persons) who subsequently withdrew consent.
It is also necessary to handle the somewhat awkward format of array fields.

"""

###############################################################################
# Notes

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
import partner.utility as putility
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
    if True:

        paths["parameters"] = os.path.join(
            path_dock, "parameters", "psychiatric_metabolism",
        )


        paths["heritability"] = os.path.join(
            path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "5_gwas_heritability_ldsc",
        )
        paths["heritability_no_liability"] = os.path.join(
            path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "5_gwas_heritability_ldsc_no_liability",
        )
        paths["heritability_extraction"] = os.path.join(
            path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "5_gwas_heritability_ldsc_extraction",
        )


        paths["correlation_primary_secondary"] = os.path.join(
            path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "6_gwas_correlation_ldsc_primary_secondary",
        )
        paths["correlation_primary_secondary_extraction"] = os.path.join(
            path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "6_gwas_correlation_ldsc_primary_secondary_extraction",
        )


        paths["correlation_primary"] = os.path.join(
            path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "6_gwas_correlation_ldsc_primary",
        )
        paths["correlation_primary_extraction"] = os.path.join(
            path_dock, "gwas_2023-12-30_ldsc_2024-01-08", "6_gwas_correlation_ldsc_primary_extraction",
        )


    # Remove previous files to avoid version or batch confusion.
    if restore:
        putility.remove_directory(path=paths["heritability_extraction"])
        putility.remove_directory(path=paths["correlation_primary_secondary_extraction"])
        putility.remove_directory(path=paths["correlation_primary_extraction"])
    # Initialize directories.
    putility.create_directories(
        path=paths["heritability_extraction"]
    )
    putility.create_directories(
        path=paths["correlation_primary_secondary_extraction"]
    )
    putility.create_directories(
        path=paths["correlation_primary_extraction"]
    )
    # Return information.
    return paths



##########
# Write


def write_product_table(
    name=None,
    table=None,
    path_directory=None,
):
    """
    Writes product information in a Pandas data frame table to file.

    arguments:
        name (str): base name for file
        table (object): Pandas data frame table
        path_directory (str): path to directory in which to write file

    raises:

    returns:

    """

    # Specify directories and files.
    path_table_text = os.path.join(
        path_directory, str(name + ".tsv")
    )
    # Write information to file.
    table.to_csv(
        path_or_buf=path_table_text,
        sep="\t",
        header=True,
        index=False,
        na_rep="NA",
    )
    pass


def control_write_product(
    pail_write=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict): collection of information to write to file
        path_parent (str): full path to directory in which to write files

    raises:

    returns:

    """

    # Initialize directories.
    putility.create_directories(
        path=path_directory,
    )
    # Write each table to file.
    for name in pail_write.keys():
        write_product_table(
            name=name,
            table=pail_write[name],
            path_directory=path_directory,
        )
    pass


###############################################################################
# Procedure


# TODO: TCW; 11 January 2024
# TODO: I need to extract the primary-primary genetic correlations...


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
    putility.print_terminal_partition(level=1)
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
        list_sort = putility.read_file_text_list(
            delimiter="\n",
            path_file=path_file_list_sort,
        )
        print("Length of list of indices: " + str(len(list_sort)))
        print(list_sort)
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
                putility.sort_table_rows_by_list_indices(
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
                putility.sort_table_rows_by_list_indices(
                    table=table_heritability_no_liability,
                    list_sort=list_sort,
                    name_column="sort_rows_temporary",
                    report=True,
            ))
        pail_write_heritability["table_heritability_no_liability"] = table_heritability_no_liability
        # Write information to file.
        control_write_product(
            pail_write=pail_write_heritability,
            path_directory=paths["heritability_extraction"],
        )


    # TODO: TCW; 19 January 2024
    # TODO: for the genetic correlation extractions, implement optional feature
    # to sort the rows in each table...

    ##########
    # Manage extraction of information about genetic correlations between
    # pairs of primary and secondary traits.

    if True:
        # Optionally read list of indices by which to sort rows in each
        # extraction table of genetic correlations.
        path_file_list_sort = os.path.join(
            paths["parameters"], "list_sort_correlation_primary_secondary.txt",
        )
        #list_sort = putility.read_file_text_list(
        #    delimiter="\n",
        #    path_file=path_file_list_sort,
        #)
        # Collect information.
        pail_write_correlation = dict()
        # Extract names of child directories within parent directory.
        names_directories = putility.extract_subdirectory_names(
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
            if False:
                table_correlation = (
                    putility.sort_table_rows_by_list_indices(
                        table=table_correlation,
                        list_sort=list_sort,
                        name_column="sort_rows_temporary",
                        report=True,
                ))
            pail_write_correlation[str("table_" + name_directory)] = table_correlation
        # Write information to file.
        control_write_product(
            pail_write=pail_write_correlation,
            path_directory=paths["correlation_primary_secondary_extraction"],
        )

    ##########
    # Manage extraction of information about genetic correlations between
    # pairs of primary and secondary traits.

    if True:
        # Optionally read list of indices by which to sort rows in each
        # extraction table of genetic correlations.
        path_file_list_sort = os.path.join(
            paths["parameters"], "list_sort_correlation_primary.txt",
        )
        #list_sort = putility.read_file_text_list(
        #    delimiter="\n",
        #    path_file=path_file_list_sort,
        #)
        # Collect information.
        pail_write_correlation = dict()
        # Extract names of child directories within parent directory.
        names_directories = putility.extract_subdirectory_names(
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
            if False:
                table_correlation = (
                    putility.sort_table_rows_by_list_indices(
                        table=table_correlation,
                        list_sort=list_sort,
                        name_column="sort_rows_temporary",
                        report=True,
                ))
            pail_write_correlation[str("table_" + name_directory)] = table_correlation
        # Write information to file.
        control_write_product(
            pail_write=pail_write_correlation,
            path_directory=paths["correlation_primary_extraction"],
        )



    pass





#
