
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
    #path_table_reference_shin_2014 = os.path.join(
    #    path_dock, "access", "metabolites", "metaboliteMap.txt"
    #)
    path_table_reference_panyard_2021 = os.path.join(
        path_dock, "parameters", "bipolar_metabolism", "metabolite_reference",
        "33437055_panyard_2021", "table_metabolite_reference.tsv"
    )
    # Read information from file.
    table_reference_panyard_2021 = pandas.read_csv(
        path_table_reference_panyard_2021,
        sep="\t",
        header=0,
        dtype="string",
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_reference_panyard_2021)
        utility.print_terminal_partition(level=2)
    # Compile and return information.
    return {
        "table_reference_panyard_2021": table_reference_panyard_2021,
    }


def organize_metabolite_reference_table(
    table=None,
    identifier=None,
    name=None,
):
    """
    Organizes information about general attributes.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        identifier (str): name of column for metabolite identifier
        name (str): name of column for metabolite biochemical name

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Copy data.
    table = table.copy(deep=True)
    # Translate column names.
    translations = dict()
    translations[identifier] = "identifier"
    translations[name] = "name"
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Select relevant columns.
    table = table.loc[
        :, table.columns.isin(["identifier", "name"])
    ]
    # Organize table.
    table.reset_index(
        level=None,
        inplace=True
    )
    table["identifier"].astype("string")
    table.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Return information.
    return table


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
    table_reference=None,
    path_parent=None,
    report=None,
):
    """
    Reads, collects, and organizes metabolite heritability estimates.

    arguments:
        table_reference (object): Pandas data frame of metabolites' identifiers
            and names from study
        path_parent (str): path to parent directory for files with heritability
            estimations for metabolites
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of metabolites' heritability estimates

    """

    # Organize metabolite reference table.
    table_reference = organize_metabolite_reference_table(
        table=table_reference,
        identifier="identifier_study",
        name="name",
    )

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
    # Organize heritability table.
    table = utility.convert_records_to_dataframe(
        records=records
    )
    table.sort_values(
        by=["heritability"],
        axis="index",
        ascending=False,
        inplace=True,
    )
    table.reset_index(
        level=None,
        inplace=True
    )
    table["identifier"].astype("string")
    table.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_merge = table.merge(
        table_reference,
        how="outer",
        left_on="identifier",
        right_on="identifier",
        suffixes=("_heritability", "_reference"),
    )
    columns_sequence = [
        #"identifier",
        "name",
        "heritability", "heritability_standard_error",
        "ratio", "ratio_standard_error",
        "variants",
    ]
    table_merge = table_merge[[*columns_sequence]]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(path_parent)
        print(table_merge)
    # Return information.
    return table_merge


def read_collect_organize_metabolites_heritabilities_studies(
    table_reference_panyard_2021=None,
    paths=None,
    report=None,
):
    """
    Reads, collects, and organizes metabolite heritability estimates.

    arguments:
        table_reference_panyard_2021 (object): Pandas data frame of metabolites'
            identifiers and names from study
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
        table_reference=table_reference_panyard_2021,
        path_parent=paths["heritability_panyard_2021"],
        report=report,
    )
    # Return information.
    return pail


def write_product_study_table(
    name=None,
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        name (str): base name for file
        information (object): information to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_table = os.path.join(
        path_parent, str(name + ".tsv")
    )
    # Write information to file.
    information.to_csv(
        path_or_buf=path_table,
        sep="\t",
        header=True,
        index=False,
    )
    pass


def write_product_studies(
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

    for name in information.keys():
        write_product_study_table(
            name=name,
            information=information[name],
            path_parent=path_parent,
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

    # Cohort tables in PLINK format.
    write_product_studies(
        information=information["studies"],
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
    # Read source information from file.
    source = read_source(
        path_dock=path_dock,
        report=True,
    )
    # Read and collect heritability estimations for metabolites from multiple
    # GWAS.
    pail_studies = read_collect_organize_metabolites_heritabilities_studies(
        table_reference_panyard_2021=source["table_reference_panyard_2021"],
        paths=paths,
        report=True,
    )

    # Collect information.
    information = dict()
    information["studies"] = pail_studies
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )

    pass



if (__name__ == "__main__"):
    execute_procedure()
