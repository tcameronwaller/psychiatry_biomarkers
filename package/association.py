
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
    paths["association"] = os.path.join(path_dock, "association")
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["association"])
    # Initialize directories.
    utility.create_directories(
        path=paths["association"]
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
        path_dock, "organization", "table_phenotypes.pickle"
    )
    path_table_metabolites = os.path.join(
        path_dock, "aggregation", "selection", "table_metabolites_scores.pickle"
    )
    # Read information from file.
    table_phenotypes = pandas.read_pickle(
        path_table_phenotypes
    )
    table_metabolites = pandas.read_pickle(
        path_table_metabolites
    )
    # Compile and return information.
    return {
        "table_phenotypes": table_phenotypes,
        "table_metabolites": table_metabolites,
    }


##########
# Iterate on metabolites...

def organize_regress_metabolites_genetic_scores_against_phenotypes(
    phenotype=None,
    metabolites=None,
    table_metabolites_scores=None,
    table_phenotypes=None,
    report=None,
):
    """
    Organizes information and regresses metabolites' genetic scores against
    phenotypes across the UK Biobank.

    arguments:
        phenotype (str): name of column in phenotype table for variable to set
            as dependent variable in regressions
        metabolites (list<str>): identifiers of metabolites for which to regress
            their genetic scores against phenotypes across UK Biobank
        table_metabolites_scores (object): Pandas data frame of metabolites'
            genetic scores across UK Biobank cohort
        table_phenotypes (object): Pandas data frame of phenotype variables
            across UK Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): information from regressions

    """

    pass


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
    print(source["table_metabolites"])
    print(source["table_phenotypes"])
    # Regress associations between metabolites' genetic scores and phenotypes
    # accross the UK Biobank.
    # M00599: pyruvate
    # M32315: serine
    # M02342: serotonin
    # M00054: tryptophan
    pail_association = (
        organize_regress_metabolites_genetic_scores_against_phenotypes(
            phenotype="body_mass_index",
            metabolites=["M00599", "M32315", "M02342", "M00054"],
            table_metabolites_scores=source["table_metabolites"],
            table_phenotypes=source["table_phenotypes"],
            report=True,
    ))

    # TODO: call function to iterate on metabolites...
    # read metabolite identifiers from table columns...
    # merge metabolite scores with phenotype table
    # drop rows with missing values in relevant variables...
    # z-score standardize metabolite scores and relevant variables...


    if False:
        # Collect information.
        information = dict()
        information["table_phenotypes"] = table_basis
        # Write product information to file.
        write_product(
            paths=paths,
            information=information
        )
    pass


if (__name__ == "__main__"):
    execute_procedure()
