
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
import uk_biobank.stratification as ukb_strat



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
    paths["stratification"] = os.path.join(path_dock, "stratification")
    paths["cohorts_models_linear"] = os.path.join(
        path_dock, "stratification", "cohorts_models_linear"
    )
    paths["cohorts_models_logistic"] = os.path.join(
        path_dock, "stratification", "cohorts_models_logistic"
    )

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["stratification"])
    # Initialize directories.
    utility.create_directories(
        path=paths["stratification"]
    )
    utility.create_directories(
        path=paths["cohorts_models_linear"]
    )
    utility.create_directories(
        path=paths["cohorts_models_logistic"]
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
        path_dock, "organization",
        "table_phenotypes.pickle",
    )

    # Read information from file.
    table_phenotypes = pandas.read_pickle(
        path_table_phenotypes
    )
    # Compile and return information.
    return {
        "table_phenotypes": table_phenotypes,
        #"table_ukb_samples": table_ukb_samples,
    }


##########
# Write


def write_product_cohort_model_table(
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


def write_product_cohorts_models(
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
        write_product_cohort_model_table(
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
    write_product_cohorts_models(
        information=information["cohorts_models_linear"],
        path_parent=paths["cohorts_models_linear"],
    )
    write_product_cohorts_models(
        information=information["cohorts_models_logistic"],
        path_parent=paths["cohorts_models_logistic"],
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

    # Select and organize variables across cohorts.
    # Organize phenotypes and covariates in format for analysis in PLINK.
    # else: pail_cohorts_models = dict()
    pail_cohorts_models_linear = (
        ukb_strat.execute_stratify_genotype_cohorts_plink_format_set(
            table=source["table_phenotypes"],
            set="bipolar_body_linear",
            path_dock=path_dock,
            report=True,
    ))
    pail_cohorts_models_logistic = (
        ukb_strat.execute_stratify_genotype_cohorts_plink_format_set(
            table=source["table_phenotypes"],
            set="bipolar_body_logistic",
            path_dock=path_dock,
            report=True,
    ))
    # Collect information.
    information = dict()
    information["cohorts_models_linear"] = pail_cohorts_models_linear
    information["cohorts_models_logistic"] = pail_cohorts_models_logistic
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )

    pass



if (__name__ == "__main__"):
    execute_procedure()
