
"""
...

This module collects and organizes information about heritability estimates for
metabolites.

Source of GWAS summary statistics is Shin et al, Nature Genetics, 2014
(PubMed:24816252). Metabolite identifiers are from Metabolon Inc.

Method for estimation of heritability from GWAS summary statistics was linkage
disequilibrium score (LDSC) regression in tool LD SCore
(https://github.com/bulik/ldsc).

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

# Relevant

import numpy
import pandas
import scipy.stats


# Custom
import promiscuity.utility as utility


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
    paths["heritability"] = os.path.join(path_dock, "heritability")
    paths["heritability_panyard_2021"] = os.path.join(
        path_dock, "heritability", "33437055_panyard_2021"
    )
    paths["collection"] = os.path.join(
        path_dock, "heritability", "collection"
    )
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["collection"])
    # Initialize directories.
    utility.create_directories(
        path=paths["collection"]
    )
    # Return information.
    return paths


def read_source(
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_table_metabolite_names = os.path.join(
        path_dock, "access", "metabolites", "metaboliteMap.txt"
    )
    path_metabolite_files = os.path.join(
        path_dock, "access", "metabolites", "metabolite_files.txt"
    )
    # Read information from file.
    table_metabolite_names = pandas.read_csv(
        path_table_metabolite_names,
        sep="\t",
        header=0,
        dtype="string",
    )
    metabolite_files = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_metabolite_files,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_metabolite_names)
        utility.print_terminal_partition(level=2)
        print("count of metabolites: " + str(len(metabolite_files)))
        print("example: " + str(metabolite_files[10]))
        utility.print_terminal_partition(level=2)
    # Compile and return information.
    return {
        "table_metabolite_names": table_metabolite_names,
        "metabolite_files": metabolite_files,
    }


def extract_metabolite_file_identifiers(
    metabolite_files=None,
):
    """
    Extracts metabolite identifiers from names of files.

    arguments:
        metabolite_files (list<str>): names of files with metabolite GWAS
            summary statistics

    raises:

    returns:
        (list<str>): identifiers of metabolites

    """

    # Iterate across metabolite file names.
    # Collect metabolite identifiers.
    identifiers = list()
    for file in metabolite_files:
        components = file.split(".")
        identifier = str(components[0]).strip()
        if ((len(identifier) > 1) and (identifier[0] == "M")):
            identifiers.append(identifier)
            pass
        pass
    # Make sure that all identifiers are unique.
    identifiers_unique = utility.collect_unique_elements(
        elements=identifiers,
    )
    # Return information.
    return identifiers_unique


def merge_metabolite_names_heritabilities(
    table_names=None,
    table_heritabilities=None,
    report=None,
):
    """
    Merges tables with information about metabolite names and heritabilities.

    arguments:
        table_names (object): Pandas data frame of metabolite identifiers and
            names
        table_heritabilities (object): Pandas data frame of metabolite
            heritability estimates
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of metabolite identifiers, names, and
            heritability estimates

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_names)
        utility.print_terminal_partition(level=2)
        print(table_heritabilities)
    # Organize data.
    table_names.astype("string")
    table_names.rename(
        columns={
            "metabolonID": "identifier",
            "metabolonDescription": "name",
        },
        inplace=True,
    )
    table_names.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    table_heritabilities["identifier"].astype("string")
    table_heritabilities.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_merge = table_names.merge(
        table_heritabilities,
        how="outer",
        left_on="identifier",
        right_on="identifier",
        suffixes=("_name", "_heritability"),
    )
    # Remove excess columns.

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_merge)
    # Return information.
    return table_merge





def read_extract_metabolite_heritability(
    file=None,
    path_parent=None,
):
    """
    Reads and extracts information from log of LDSC for heritability estimation
    from GWAS summary statistics.

    arguments:
        file (str): name of a file
        path_parent (str): path to parent directory for files with heritability
            estimations for metabolites

    raises:

    returns:
        (dict): information about estimation of a metabolite's heritability

    """

    # Extract metabolite's identifier.
    identifier = str(
        file.replace("heritability_", "").replace(".log", "")
    )
    # Define path to file.
    path_file = os.path.join(
        path_parent, file
    )
    # Initialize variables.
    variants = float("nan")
    heritability = float("nan")
    heritability_error = float("nan")
    ratio = float("nan")
    ratio_error = float("nan")
    # Read relevant lines from file.
    lines = utility.read_file_text_lines(
        path_file=path_file,
        start=22,
        stop=30,
    )
    # Extract information from lines.
    prefix_variants = "After merging with regression SNP LD, "
    suffix_variants = " SNPs remain."
    prefix_heritability = "Total Observed scale h2: "
    prefix_ratio = "Ratio: "
    for line in lines:
        if prefix_variants in line:
            variants = float(
                line.replace(prefix_variants, "").replace(suffix_variants, "")
            )
        elif prefix_heritability in line:
            content = line.replace(prefix_heritability, "")
            contents = content.split(" (")
            heritability = float(contents[0])
            heritability_error = float(
                contents[1].replace(")", "")
            )
            pass
        elif prefix_ratio in line:
            content = line.replace(prefix_ratio, "")
            contents = content.split(" (")
            ratio = float(contents[0])
            ratio_error = float(
                contents[1].replace(")", "")
            )
            pass
        pass
    # Collect information.
    record = dict()
    record["identifier"] = identifier
    record["variants"] = variants
    record["heritability"] = heritability
    record["heritability_standard_error"] = heritability_error
    record["ratio"] = ratio
    record["ratio_standard_error"] = ratio_error
    # Return information.
    return record


def read_collect_metabolites_heritabilities(
    path_parent=None,
    report=None,
):
    """
    Reads, collects, and organizes metabolite heritability estimates.

    arguments:
        path_parent (str): path to parent directory for files with heritability
            estimations for metabolites
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of metabolites' heritability estimates

    """

    # Collect names of files for metabolites' heritabilities.
    files = utility.extract_directory_file_names(path=path_parent)
    files_heritability = list(filter(
        lambda content: ("heritability" in content), files
    ))
    records = list()
    for file in files_heritability:
        record = read_extract_metabolite_heritability(
            file=file,
            path_parent=path_parent,
        )
        records.append(record)
        pass
    # Organize table.
    table = utility.convert_records_to_dataframe(
        records=records
    )
    table.sort_values(
        by=["heritability"],
        axis="index",
        ascending=False,
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(path_parent)
        print(table)
    # Return information.
    return table


def read_collect_metabolites_heritabilities_studies(
    paths=None,
    report=None,
):
    """
    Reads, collects, and organizes metabolite heritability estimates.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): heritability estimations for metabolites from multiple GWAS

    """

    # Collect metabolites' heritabilities from each GWAS.
    pail = dict()
    #pail["table_shin_2014"] =
    pail["table_panyard_2021"] = read_collect_metabolites_heritabilities(
        path_parent=paths["heritability_panyard_2021"],
        report=report,
    )
    # Return information.
    return pail




def write_product_collection(
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_table_metabolite_heritabilities = os.path.join(
        path_parent, "table_metabolite_heritabilities.tsv"
    )
    # Write information to file.
    information["table_metabolite_heritabilities"].to_csv(
        path_or_buf=path_table_metabolite_heritabilities,
        sep="\t",
        header=True,
        index=True,
    )

    pass


def write_product(
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Alcohol consumption.
    write_product_collection(
        information=information["collection"],
        path_parent=paths["collection"],
    )
    pass


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

    utility.print_terminal_partition(level=1)
    print(path_dock)
    print("version check: 1")

    # Initialize directories.
    paths = initialize_directories(
        restore=True,
        path_dock=path_dock,
    )
    if False:
        # Read source information from file.
        source = read_source(
            path_dock=path_dock,
            report=True,
        )

    # Read and collect heritability estimations for metabolites from multiple
    # GWAS.
    pail_collection = read_collect_metabolites_heritabilities_studies(
        paths=paths,
        report=True,
    )

    if False:
        # Extract identifiers of metabolites with GWAS summary statistics.
        metabolite_identifiers = extract_metabolite_file_identifiers(
            metabolite_files=source["metabolite_files"],
        )
        # Collect linkage disequilibrium score regression heritability estimates
        # for each metabolite.
        table_heritabilities = read_collect_metabolite_heritabilities(
            metabolite_identifiers=metabolite_identifiers,
            path_dock=path_dock,
            report=True,
        )
        # Merge metabolite heritabilities with metabolite names.
        table_names_heritabilities = merge_metabolite_names_heritabilities(
            table_names=source["table_metabolite_names"],
            table_heritabilities=table_heritabilities,
            report=True,
        )

        # Collect information.
        information = dict()
        information["collection"] = dict()
        information["collection"]["table_metabolite_heritabilities"] = (
            table_names_heritabilities
        )
        # Write product information to file.
        write_product(
            paths=paths,
            information=information
        )

    pass



if (__name__ == "__main__"):
    execute_procedure()
