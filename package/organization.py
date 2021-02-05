
"""
...

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
import promiscuity.utility as utility
import uk_biobank.assembly
import uk_biobank.organization

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
    paths["organization"] = os.path.join(path_dock, "organization")
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["organization"])
    # Initialize directories.
    utility.create_directories(
        path=paths["organization"]
    )
    # Return information.
    return paths


##########
# Read


def read_source(
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_table_phenotypes = os.path.join(
        path_dock, "assembly", "table_phenotypes.pickle"
    )
    path_table_metabolites_names = os.path.join(
        path_dock, "access", "24816252_shin_2014", "metaboliteMap.txt"
    )

    # Read information from file.
    table_phenotypes = pandas.read_pickle(
        path_table_phenotypes
    )
    table_metabolites_names = pandas.read_csv(
        path_table_metabolites_names,
        sep="\t",
        header=0,
        dtype="string",
    )
    # Compile and return information.
    return {
        "table_phenotypes": table_phenotypes,
        "table_metabolites_names": table_metabolites_names,
    }


##########
# Metabolites


def determine_metabolite_valid_identity(
    name=None,
):
    """
    Determine whether a single metabolite has a valid identity from Metabolon.

    arguments:
        name (str): name of metabolite from Metabolon reference

    raises:

    returns:
        (float): ordinal representation of person's frequency of alcohol
            consumption

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (len(str(name)) > 2):
        # The variable has a valid value.
        if (str(name).strip().lower().startswith("x-")):
            # Metabolite has an indefinite identity.
            identity = 0
        else:
            # Metabolite has a definite identity.
            identity = 1
    else:
        # Name is empty.
        #identity = float("nan")
        identity = 0
    # Return information.
    return identity


def select_metabolites_with_valid_identities(
    table=None,
    report=None,
):
    """
    Selects identifiers of metabolites from Metabolon with valid identities.

    arguments:
        table (object): Pandas data frame of metabolite identifiers and names
            from Metabolon
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about metabolites, their identifiers,
            and their names

    """

    # Copy information.
    table = table.copy(deep=True)
    # Translate column names.
    translations = dict()
    translations["metabolonID"] = "identifier"
    translations["metabolonDescription"] = "name"
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Determine whether metabolite has a valid identity.
    table["identity"] = table.apply(
        lambda row:
            determine_metabolite_valid_identity(
                frequency=row["name"],
            ),
        axis="columns", # apply across rows
    )
    # Select metabolites with valid identities.
    table_valid = table.loc[
        (table["identity"] > 0.5), :
    ]
    metabolites_valid = table_valid["identifier"].to_list()
    names_valid = table_valid["name"].to_list()
    # Compile information.
    pail = dict()
    pail["table_metabolites_names"] = table
    pail["metabolites_valid"] = list()
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("Report from select_metabolites_with_valid_identities()")
        utility.print_terminal_partition(level=3)
        print(names_valid)
        utility.print_terminal_partition(level=3)
        print(table)
    # Return information.
    return pail


##########
# Write


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

    returns:

    """

    # Specify directories and files.
    path_table_phenotypes = os.path.join(
        paths["organization"], "table_phenotypes.pickle"
    )
    path_table_phenotypes_text = os.path.join(
        paths["organization"], "table_phenotypes.tsv"
    )
    # Write information to file.
    information["table_phenotypes"].to_pickle(
        path_table_phenotypes
    )
    information["table_phenotypes"].to_csv(
        path_or_buf=path_table_phenotypes_text,
        sep="\t",
        header=True,
        index=True,
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
    # Pause procedure.
    time.sleep(5.0)

    # Initialize directories.
    paths = initialize_directories(
        restore=True,
        path_dock=path_dock,
    )
    # Read source information from file.
    # Exclusion identifiers are "eid".
    source = read_source(
        path_dock=path_dock,
        report=True,
    )

    # Select metabolites with valid identities.
    metabolites_valid = select_metabolites_with_valid_identities(
        table_metabolites_names=source["table_metabolites_names"],
    )

    # Organize variables for basic characteristics, genotypes, and hormones
    # across the UK Biobank.
    table_basis = uk_biobank.organization.organize_basic_characteristics(
        table=source["table_phenotypes"],
        report=True,
    )

    # Organize variables for alcohol consumption across the UK Biobank.
    table_alcohol = uk_biobank.organization.organize_alcohol_consumption(
        table=table_basis,
        report=True,
    )

    # TODO: Adapt the ICD9 and ICD10 functionality for depression and bipolar...

    # Collect information.
    information = dict()
    information["metabolites_valid"] = metabolites_valid
    information["table_phenotypes"] = table_alcohol
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )
    pass


if (__name__ == "__main__"):
    execute_procedure()
