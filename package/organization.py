
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

###############################################################################
# Functionality


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


def select_organize_metabolites_valid_identities_scores(
    table_names=None,
    table_scores=None,
    report=None,
):
    """
    Selects identifiers of metabolites from Metabolon with valid identities.

    arguments:
        table_names (object): Pandas data frame of metabolites' identifiers and
            names from Metabolon
        table_scores (object): Pandas data frame of metabolites' genetic scores
            across UK Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about metabolites, their identifiers,
            and their names

    """

    # Copy information.
    table_names = table_names.copy(deep=True)
    table_scores = table_scores.copy(deep=True)
    # Translate column names.
    translations = dict()
    translations["metabolonID"] = "identifier"
    translations["metabolonDescription"] = "name"
    table_names.rename(
        columns=translations,
        inplace=True,
    )
    # Determine whether metabolite has a valid identity.
    table_names["identity"] = table_names.apply(
        lambda row:
            determine_metabolite_valid_identity(
                name=row["name"],
            ),
        axis="columns", # apply across rows
    )
    # Select metabolites with valid identities.
    table_identity = table_names.loc[
        (table_names["identity"] > 0.5), :
    ]
    metabolites_identity = table_identity["identifier"].to_list()
    names_identity = table_identity["name"].to_list()
    # Organize table.
    table_names["identifier"].astype("string")
    table_names.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Remove table columns for metabolites with null genetic scores.
    table_scores.dropna(
        axis="columns",
        how="all",
        subset=None,
        inplace=True,
    )
    # Select metabolites with valid identities and valid genetic scores.
    metabolites_scores = table_scores.columns.to_list()
    metabolites_valid = utility.filter_common_elements(
        list_minor=metabolites_identity,
        list_major=metabolites_scores,
    )
    # Compile information.
    pail = dict()
    pail["table"] = table_names
    pail["metabolites_valid"] = metabolites_valid
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("Report from select_metabolites_with_valid_identities()")
        utility.print_terminal_partition(level=3)
        print(
            "Count of identifiable metabolites: " +
            str(len(metabolites_identity))
        )
        print(
            "Count of identifiable metabolites with scores: " +
            str(len(metabolites_valid))
        )
        utility.print_terminal_partition(level=3)
        print(table_names)
    # Return information.
    return pail


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

    utility.print_terminal_partition(level=1)
    print(path_dock)
    print("version check: 1")
    # Pause procedure.
    time.sleep(5.0)

    # Execute assembly procedure from uk_biobank package.
    uk_biobank.organization.execute_procedure(path_dock=path_dock)
    utility.print_terminal_partition(level=1)
    print("From package 'uk_biobank', procedure 'organization' is complete.")
    pass


if (__name__ == "__main__"):
    execute_procedure()
