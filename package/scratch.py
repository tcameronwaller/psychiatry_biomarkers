
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
import partner.utility as utility
import partner.plot as plot
import uk_biobank.organization as ukb_organization
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
    paths["scratch"] = os.path.join(path_dock, "scratch")

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["scratch"])
    # Initialize directories.
    utility.create_directories(
        path=paths["scratch"]
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

    # Scrap...

    #####################
    table_cases = source["table_phenotypes"].loc[
        (
            (source["table_phenotypes"]["bipolar_control_case_strict"] == 1)
        ), :
    ]
    print("CASES!!!")
    print(table_cases.shape[0])
    table_controls = source["table_phenotypes"].loc[
        (
            (source["table_phenotypes"]["bipolar_control_case_strict"] == 0)
        ), :
    ]
    print("CONTROLS!!!")
    print(table_controls.shape[0])
    #######################3


    # Specify directories and files.
    path_table = os.path.join(
        path_dock, "stratification", "cohorts_models_logistic",
        "table_white_bipolar_control_case_strict.tsv",
    )

    # Read information from file.
    table_phenotypes = pandas.read_csv(
        path_table,
        sep="\t",
        header=0,
        dtype={
            "FID": "string",
            "IID": "string",
            "eid": "string",
            "bipolar_control_case_strict": numpy.int32,
        },
    )
    table_phenotypes.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table_phenotypes.set_index(
        "eid",
        append=False,
        drop=True,
        inplace=True
    )
    # Compile and return information.
    return {
        "table_phenotypes": table_phenotypes,
        #"table_ukb_samples": table_ukb_samples,
    }


##########
# Write



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
    print("version check: 21")
    # Pause procedure.
    time.sleep(5.0)

    # Initialize directories.
    paths = initialize_directories(
        restore=True,
        path_dock=path_dock,
    )

    path_table_kinship_pairs = os.path.join(
        path_dock, "access", "ukbiobank_phenotypes",
        "table_kinship_pairs.dat"
    )
    table_kinship_pairs = pandas.read_csv(
        path_table_kinship_pairs,
        sep="\s+",
        header=0,
        dtype={
            "ID1": "string",
            "ID2": "string",
            "HetHet": "float32",
            "IBS0": "float32",
            "Kinship": "float32",
        },
    )
    path_table_kinship_pairs = os.path.join(
        path_dock, "assembly", "table_kinship_pairs.pickle"
    )
    path_table_kinship_pairs_text = os.path.join(
        path_dock, "assembly", "table_kinship_pairs.tsv"
    )
    table_kinship_pairs.to_pickle(
        path_table_kinship_pairs
    )
    table_kinship_pairs.to_csv(
        path_or_buf=path_table_kinship_pairs_text,
        sep="\t",
        header=True,
        index=False,
    )

    # Read source information from file.
    # Read source information from file.
    table_kinship_pairs = ukb_strat.read_source_table_kinship_pairs(
        path_dock=path_dock,
        report=True,
    )
    print(table_kinship_pairs)


    pass



if (__name__ == "__main__"):
    execute_procedure()
