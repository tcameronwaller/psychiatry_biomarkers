
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

# Relevant

import numpy
import pandas
import scipy.stats


# Custom
import promiscuity.utility as utility
import promiscuity.plot as plot

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
    paths["aggregation"] = os.path.join(path_dock, "aggregation")
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["aggregation"])
    # Initialize directories.
    utility.create_directories(
        path=paths["aggregation"]
    )
    # Return information.
    return paths


##########
# Read


def filter_files_names_metabolites_genetic_scores(
    pattern=None,
    files=None,
):
    """
    Filters to file names for metabolites' genetic scores.

    arguments:
        pattern (str): string pattern to recognize files
        files (list<str>): names of files

    raises:

    returns:
        (list<str>): unique file names for metabolites' genetic scores

    """

    # Filter to names of files with genetic scores.
    files_scores = list(filter(
        lambda file: (
            (pattern in str(file))
        ), files
    ))
    files_scores_unique = utility.collect_unique_elements(
        elements=files_scores,
    )
    return files_scores_unique


def extract_metabolite_identifier_from_file_name(
    file=None,
):
    """
    Extracts metabolite identifier from a file name.

    arguments:
        file (str): name of a file with genetic scores for a metabolite

    raises:

    returns:
        (str): metabolite's identifier

    """

    # Extract metabolite identifier.
    components_period = file.split(".")
    prefix = str(components_period[0]).strip()
    components_underscore = prefix.split("_")
    identifier = str(components_underscore[1]).strip()
    # Return information.
    return identifier


def associate_metabolites_identifiers_files_paths(
    files=None,
    path_parent=None,
):
    """
    Collects and associates metabolites' identifiers to their paths and file
    names.

    arguments:
        files (list<str>): names of files
        path_parent (str): path to parent directory that contains files

    raises:

    returns:
        (dict<list<str>>): collection of files and paths for metabolites

    """

    # Copy.
    files = copy.deepcopy(files)
    # Filter to files for metabolites' genetic scores.
    files_scores = filter_files_names_metabolites_genetic_scores(
        pattern="all.score.gz",
        files=files,
    )
    # Collect and associate metabolites' identifiers to files' paths and names.
    metabolites_files_paths = dict()
    for file in files_scores:
        file_path = os.path.join(path_parent, file)
        identifier = extract_metabolite_identifier_from_file_name(
            file=file,
        )
        if ((len(identifier) > 1) and (identifier[0] == "M")):
            metabolites_files_paths["identifier"] = dict()
            metabolites_files_paths["identifier"]["metabolite"] = identifier
            metabolites_files_paths["identifier"]["file"] = file
            metabolites_files_paths["identifier"]["path"] = file_path
        pass
    # Return information.
    return metabolites_files_paths


def read_source_initial(
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
    path_metabolite_gene_scores = os.path.join(
        path_dock, "access", "ukbiobank_metabolites"
    )
    # Read information from file.
    files = utility.extract_directory_file_names(
        path=path_metabolite_gene_scores,
    )
    # Associate metabolites' identifiers, file names, and file paths.
    metabolites_files_paths = associate_metabolites_identifiers_files_paths(
        files=files,
        path_parent=path_metabolite_gene_scores,
    )
    # Compile and return information.
    return {
        "metabolites_files_paths": metabolites_files_paths,
    }


def read_source_metabolite_genetic_scores(
    path_file=None,
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
    path_metabolite_gene_scores = os.path.join(
        path_dock, "access", "ukbiobank_metabolites"
    )
    # Read information from file.
    files = utility.extract_directory_file_names(
        path=path_metabolite_gene_scores,
    )
    files_paths = list()
    for file in files:
        file_path = os.path.join(path_metabolite_gene_scores, file)
        files_paths.append(file_path)
    # Compile and return information.
    return {
        "files": files,
        "files_paths": files_paths,
    }


def read_collect_aggregate_metabolites_genetic_scores(
    metabolites_files_paths=None,
):
    """
    Reads metabolites' genetic scores across the UK Biobank from file,
    aggregates scores by Singular Value Decomposition (SVD), and collects these
    within a table.

    arguments:
        metabolites_files_paths (dict<list<str>>): collection of files and paths
            for metabolites

    raises:

    returns:
        (dict): collection of information

    """

    # Compile information.
    pail = dict()
    # Return information.

    pass



##########
# Write


def write_product_quality(
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

    # Specify directories and files.
    path_table_auditc = os.path.join(
        path_parent, "table_auditc.tsv"
    )
    path_table_audit = os.path.join(
        path_parent, "table_audit.tsv"
    )
    path_table_diagnosis = os.path.join(
        path_parent, "table_diagnosis.tsv"
    )
    path_table_alcoholism = os.path.join(
        path_parent, "table_alcoholism.tsv"
    )
    # Write information to file.
    information["table_auditc"].to_csv(
        path_or_buf=path_table_auditc,
        sep="\t",
        header=True,
        index=False,
    )
    information["table_audit"].to_csv(
        path_or_buf=path_table_audit,
        sep="\t",
        header=True,
        index=False,
    )
    information["table_diagnosis"].to_csv(
        path_or_buf=path_table_diagnosis,
        sep="\t",
        header=True,
        index=False,
    )
    information["table_alcoholism"].to_csv(
        path_or_buf=path_table_alcoholism,
        sep="\t",
        header=True,
        index=False,
    )
    pass


def write_product_cohorts_sex_alcoholism_hormone(
    sex_text=None,
    alcoholism=None,
    hormone=None,
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        sex_text (str): textual representation of sex selection
        alcoholism (str): name of column defining alcoholism cases and controls
        hormone (str): name of column for relevant sex hormone
        information (object): information to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_table = os.path.join(
        path_parent[sex_text][alcoholism][hormone],
        "table_phenotypes_covariates.tsv"
    )
    # Write information to file.
    information[sex_text][alcoholism][hormone].to_csv(
        path_or_buf=path_table,
        sep="\t",
        header=True,
        index=False,
    )
    pass


def write_product_cohorts(
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

    sexes = ["female", "male",]
    alcoholisms = [
        "alcohol_auditc", "alcohol_audit",
        "alcoholism_1", "alcoholism_2", "alcoholism_3", "alcoholism_4",
    ]
    hormones = ["oestradiol", "testosterone",]
    for sex in sexes:
        for alcoholism in alcoholisms:
            for hormone in hormones:
                write_product_cohorts_sex_alcoholism_hormone(
                    sex_text=sex,
                    alcoholism=alcoholism,
                    hormone=hormone,
                    information=information,
                    path_parent=path_parent,
                )
    pass


def write_product_trial(
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

    # Specify directories and files.
    path_table_phenotypes = os.path.join(
        path_parent, "table_phenotypes_covariates.pickle"
    )
    path_table_phenotypes_text = os.path.join(
        path_parent, "table_phenotypes_covariates.tsv"
    )
    # Write information to file.
    information["table_phenotypes_covariates"].to_pickle(
        path_table_phenotypes
    )
    information["table_phenotypes_covariates"].to_csv(
        path_or_buf=path_table_phenotypes_text,
        sep="\t",
        header=True,
        index=False,
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

    # Quality control reports.
    write_product_quality(
        information=information["quality"],
        path_parent=paths["quality"],
    )
    # Cohort tables in PLINK format.
    write_product_cohorts(
        information=information["cohorts"],
        path_parent=paths["cohorts"],
    )

    # Trial organization.
    if False:
        write_product_trial(
            information=information["trial"],
            path_parent=paths["trial"],
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
    # Exclusion identifiers are "eid".
    source = read_source_initial(
        path_dock=path_dock,
        report=True,
    )
    print(source["metabolites_files_paths"])
    # Collect metabolites' genetic scores, and aggregate these by singular value
    # decomposition (SVD).
    pail_metabolites_scores = read_collect_aggregate_metabolites_genetic_scores(
        metabolites_files_paths=source["metabolites_files_paths"],
    )
    # Alright... now read in the scores one metabolite at a time...


    # ^^^ read in list of unique metabolite identifiers "M#####"

    # read all file names from directory
    # filter to keep only the file names for scores
    # keep this collection of file names for scores

    # split file names by "."
    # take element 0 from the split
    # split again by "_"
    # take element 1


    # then iterate on that list and read in the actual data 1 metabolite at a time

    # for each metabolite...
    # ... aggregate PRS scores by PRS-PCA
    # ... collect metabolite's PC1 PRS-PCA score with other metabolite scores


    # Collect information.
    information = dict()
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )

    pass



if (__name__ == "__main__"):
    execute_procedure()
