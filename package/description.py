
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
import promiscuity.plot as plot
import uk_biobank.organization as ukb_organization
import uk_biobank.description as ukb_description

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
    paths["description"] = os.path.join(path_dock, "description")
    paths["plots"] = os.path.join(
        path_dock, "description", "plots"
    )

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["description"])
    # Initialize directories.
    utility.create_directories(
        path=paths["description"]
    )
    utility.create_directories(
        path=paths["plots"]
    )
    # Return information.
    return paths


##########
# Read



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

    # Describe variables within cohorts and models.
    pail_summary = (
        ukb_description.execute_describe_cohorts_models_phenotypes(
            table=source["table_phenotypes"],
            genotype_cohorts=False, # genotype cohorts are slow
            set="none", # only relevant for genotype cohorts
            path_dock=path_dock,
            report=True,
    ))

    # Plot figures for cohorts, models, and phenotypes.
    if True:
        pail_plot = ukb_description.execute_plot_cohorts_models_phenotypes(
            table=source["table_phenotypes"],
            report=True,
        )
    else:
        pail_plot = dict()

    # Collect information.
    information = dict()
    information["description"] = dict()
    information["description"]["table_cohorts_measurements_missingness"] = (
        pail_summary["table_cohorts_measurements_missingness"]
    )
    information["description"]["table_cohorts_vitamin_d_deficiency"] = (
        pail_summary["table_cohorts_vitamin_d_deficiency"]
    )
    information["description"]["table_summary_cohorts_models_phenotypes"] = (
        pail_summary["table_summary_cohorts_models_phenotypes"]
    )
    information["description"]["table_summary_cohorts_models_genotypes"] = (
        pail_summary["table_summary_cohorts_models_genotypes"]
    )
    information["plots"] = pail_plot
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )
    pass


if (__name__ == "__main__"):
    execute_procedure()
