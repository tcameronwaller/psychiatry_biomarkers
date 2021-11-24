
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
    paths["reference_population"] = os.path.join(
        path_dock, "stratification", "reference_population"
    )
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
        path=paths["reference_population"]
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
    # Reference population.
    if True:
        pail_population = (
            ukb_strat.execute_stratify_genotype_cohorts_plink_format_set(
                table=source["table_phenotypes"],
                set="reference_population",
                path_dock=path_dock,
                report=True,
        ))
    else:
        pail_population = dict()
        pass
    # Cohorts and models for linear genetic analyses.
    if True:
        pail_linear = (
            ukb_strat.execute_stratify_genotype_cohorts_plink_format_set(
                table=source["table_phenotypes"],
                set="bipolar_body_linear",
                path_dock=path_dock,
                report=True,
        ))
    else:
        pail_linear = dict()
        pass
    # Cohorts and models for logistic genetic analyses.
    if True:
        pail_logistic = (
            ukb_strat.execute_stratify_genotype_cohorts_plink_format_set(
                table=source["table_phenotypes"],
                set="bipolar_body_logistic",
                path_dock=path_dock,
                report=True,
        ))
    else:
        pail_logistic = dict()
        pass

    # Describe effect of preferential selection of Bipolar Disorder cases in the
    # Kinship Filter for genetic analyses.
    # Strict Bipolar Disorder.
    table_simple = (
        pail_logistic["table_white_bipolar_control_case_strict"]
    )
    table_priority = (
        pail_logistic["table_white_bipolar_control_case_strict_priority_case"]
    )
    ukb_strat.report_kinship_filter_priority_selection(
        name="... Comparison of case priority for Strict Bipolar Disorder ...",
        priority_values=[1,],
        priority_variable="bipolar_control_case_strict",
        table_full=source["table_phenotypes"],
        table_simple=table_simple,
        table_priority=table_priority,
        report=True,
    )
    # Loose Bipolar Disorder.
    table_simple = (
        pail_logistic["table_white_bipolar_control_case_loose"]
    )
    table_priority = (
        pail_logistic["table_white_bipolar_control_case_loose_priority_case"]
    )
    ukb_strat.report_kinship_filter_priority_selection(
        name="... Comparison of case priority for Loose Bipolar Disorder ...",
        priority_values=[1,],
        priority_variable="bipolar_control_case_loose",
        table_full=source["table_phenotypes"],
        table_simple=table_simple,
        table_priority=table_priority,
        report=True,
    )

    # Collect information.
    information = dict()
    information["reference_population"] = pail_population
    information["cohorts_models_linear"] = pail_linear
    information["cohorts_models_logistic"] = pail_logistic
    # Write product information to file.
    ukb_strat.write_genotype_product(
        paths=paths,
        information=information
    )

    pass



if (__name__ == "__main__"):
    execute_procedure()
