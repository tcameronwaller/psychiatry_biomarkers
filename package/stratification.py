
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
    paths = ukb_strat.initialize_directories(
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

    # Hormones and their regulatory proteins.
    if False:
        pail_hormones_linear = (
            ukb_strat.execute_stratify_genotype_cohorts_plink_format_set(
                table=source["table_phenotypes"],
                set="hormones_linear",
                path_dock=path_dock,
                report=True,
        ))
    else:
        pail_hormones_linear = dict()
        pass
    if False:
        pail_hormones_logistic = (
            ukb_strat.execute_stratify_genotype_cohorts_plink_format_set(
                table=source["table_phenotypes"],
                set="hormones_logistic",
                path_dock=path_dock,
                report=True,
        ))
    else:
        pail_hormones_logistic = dict()
        pass

    # Body mass index (BMI) in Bipolar Disorder.
    if True:
        pail_bipolar_linear = (
            ukb_strat.execute_stratify_genotype_cohorts_plink_format_set(
                table=source["table_phenotypes"],
                set="bipolar_body_linear",
                path_dock=path_dock,
                report=True,
        ))
    else:
        pail_bipolar_linear = dict()
        pass
    if True:
        pail_bipolar_logistic = (
            ukb_strat.execute_stratify_genotype_cohorts_plink_format_set(
                table=source["table_phenotypes"],
                set="bipolar_body_logistic",
                path_dock=path_dock,
                report=True,
        ))
    else:
        pail_bipolar_logistic = dict()
        pass

    # Collect information.
    information = dict()
    information["reference_population"] = pail_population
    information["hormones_linear"] = pail_hormones_linear
    information["hormones_logistic"] = pail_hormones_logistic
    information["body_bipolar_linear"] = pail_bipolar_linear
    information["body_bipolar_logistic"] = pail_bipolar_logistic
    # Write product information to file.
    ukb_strat.write_genotype_product(
        paths=paths,
        information=information
    )

    pass



if (__name__ == "__main__"):
    execute_procedure()
