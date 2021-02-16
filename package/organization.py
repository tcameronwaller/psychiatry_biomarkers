
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
    path_table_metabolites_scores = os.path.join(
        path_dock, "aggregation", "selection",
        "table_metabolites_scores_prs_0_1.pickle"
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
    table_metabolites_scores = pandas.read_pickle(
        path_table_metabolites_scores
    )
    # Compile and return information.
    return {
        "table_phenotypes": table_phenotypes,
        "table_metabolites_names": table_metabolites_names,
        "table_metabolites_scores": table_metabolites_scores,
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

    path_table_metabolites_names = os.path.join(
        paths["organization"], "table_metabolites_names.pickle"
    )
    path_table_metabolites_names_text = os.path.join(
        paths["organization"], "table_metabolites_names.tsv"
    )
    path_metabolites_valid = os.path.join(
        paths["organization"], "metabolites_valid.pickle"
    )
    path_table_phenotypes = os.path.join(
        paths["organization"], "table_phenotypes.pickle"
    )
    path_table_phenotypes_text = os.path.join(
        paths["organization"], "table_phenotypes.tsv"
    )

    # Write information to file.
    information["table_metabolites_names"].to_pickle(
        path_table_metabolites_names
    )
    information["table_metabolites_names"].to_csv(
        path_or_buf=path_table_metabolites_names_text,
        sep="\t",
        header=True,
        index=True,
    )
    with open(path_metabolites_valid, "wb") as file_product:
        pickle.dump(information["metabolites_valid"], file_product)
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
    pail_metabolites = select_organize_metabolites_valid_identities_scores(
        table_names=source["table_metabolites_names"],
        table_scores=source["table_metabolites_scores"],
        report=True,
    )

    # Organize variables for persons' genotypes, sex, age, and body mass index
    # across the UK Biobank.
    table_basis = uk_biobank.organization.execute_genotype_sex_age_body(
        table=source["table_phenotypes"],
        report=False,
    )
    # Organize variables for persons' alcohol consumption across the UK Biobank.
    if False:
        table_alcohol = uk_biobank.organization.execute_alcohol(
            table=table_basis,
            report=True,
        )
    # Organize variables for persons' mental health across the UK Biobank.
    table_mental = uk_biobank.organization.execute_mental_health(
        table=table_basis,
        report=True,
    )
    # TODO: Adapt the ICD9 and ICD10 functionality for depression and bipolar...

    # Collect information.
    information = dict()
    information["table_metabolites_names"] = pail_metabolites["table"]
    information["metabolites_valid"] = pail_metabolites["metabolites_valid"]
    information["table_phenotypes"] = table_mental
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )
    pass


if (__name__ == "__main__"):
    execute_procedure()
