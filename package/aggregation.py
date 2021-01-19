
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
            (str(pattern) in str(file))
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
        if (
            (len(identifier) > 1) and
            (identifier[0] == "M") and
            (identifier not in metabolites_files_paths)
        ):
            metabolites_files_paths[identifier] = dict()
            metabolites_files_paths[identifier]["metabolite"] = identifier
            metabolites_files_paths[identifier]["file"] = file
            metabolites_files_paths[identifier]["path"] = file_path
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

    # Read information from file.
    variables_types = {
        "FID": "string",
        "IID": "string",
        "X5e.08": "float32",
        "X1e.07": "float32",
        "X1e.06": "float32",
        "X1e.05": "float32",
        "X0.0001": "float32",
        "X0.001": "float32",
        "X0.01": "float32",
        "X0.05": "float32",
        "X0.1": "float32",
        "X0.2": "float32",
        "X1": "float32",
    }
    table = pandas.read_csv(
        path_file,
        sep="\s+", # ",", "\t", "\s+"
        header=0,
        dtype=variables_types,
        na_values=["NA", "<NA>"],
        keep_default_na=True,
        compression=None, # "gzip"
    )
    # Report.
    if report:
        # Report only for a few metabolites.
        metabolites = ["M33627", "M22649", "M27710"]
        match = any(list(map(
            lambda metabolite: (metabolite in path_file), metabolites
        )))
        if match:
            utility.print_terminal_partition(level=2)
            print("raw table for example metabolites:")
            print(table)
            utility.print_terminal_partition(level=3)
            print(table.columns.to_list())
            utility.print_terminal_partition(level=3)
            print("variable types:")
            print(table.dtypes)
            utility.print_terminal_partition(level=3)
    # Compile and return information.
    return table


##########
# Aggregation

# raw matrix must have variables across columns and cases (samples) across rows

# PCA on covariance matrix: must center the mean (mean = 0)
# PCA on correlation matrix: must standardize scale (mean = 0, standard deviation = 1)

# use the thin method to exclude unnecessary dimensionality in component matrices


def organize_singular_value_decomposition(
    table=None,
):
    """
    Organizes a Singular Value Decomposition (SVD).

    arguments:
        table (object): Pandas data frame of variable observations across
            columns and cases or samples across rows

    raises:

    returns:
        (dict): collection of information about the SVD

    """

    # Compile information.
    pail = dict()
    # Return.
    return pail



# TODO: in progress...
def organize_aggregate_metabolite_genetic_scores(
    identifier=None,
    column_index=None,
    columns_scores=None,
    table=None,
):
    """
    Aggregates a metabolite's genetic scores across the UK Biobank by Singular
    Value Decomposition (SVD).

    This function drops any columns that neither have specification for index
    nor scores.

    arguments:
        identifier (str):
        column_index (str): name of column with UK Biobank record identifiers
        columns_scores (list<str>): names of columns with scores to include in
            aggregation
        table (object): Pandas data frame of a metabolite's raw genetic scores
            across UK Biobank

    raises:

    returns:
        (object): Pandas data frame of a metabolite's aggregate genetic scores
            across UK Biobank

    """

    # Copy information.
    table = table.copy(deep=True)
    columns_keep = copy.deepcopy(columns_scores)
    columns_keep.append(column_index)
    # Select relevant columns.
    table = table.loc[:, table.columns.isin(columns_keep)]
    # Organize table.
    table[column_index].astype("string")
    table.set_index(
        column_index,
        drop=True,
        inplace=True,
    )
    # Aggregate metabolite's genetic scores.
    #pail_decomposition = organize_singular_value_decomposition()

    # I need a table of principal components on the original scores... then I can "select"
    # PC1 here.
    # The new function should be something like... calculate majority positive PCs

    # TODO: I need to consider all values simultaneously in order to determine the new column...
    # TODO: cannot do this row-by-row

    # TODO: maybe specify an index ("identifier_ukb") before passing to new function

    pass


def read_aggregate_metabolite_genetic_scores(
    metabolite=None,
    metabolites_files_paths=None,
    report=None,
):
    """
    Reads a metabolite's genetic scores across the UK Biobank from file,
    and aggregates these scores by Singular Value Decomposition (SVD).

    This function returns a table for a single metabolite with UK Biobank
    identifiers and a single column of aggregate scores for the metabolite
    across these UK Biobank records.

    arguments:
        metabolite (str): identifier of a metabolite
        metabolites_files_paths (dict<list<str>>): collection of files and paths
            for metabolites
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of a metabolite's aggregate genetic scores
            across UK Biobank

    """

    # Read raw table of metabolite's genetic scores.
    metabolite_file_path = metabolites_files_paths[metabolite]["path"]
    table_raw = read_source_metabolite_genetic_scores(
        path_file=metabolite_file_path,
        report=True,
    )
    # Organize the raw table.
    table_raw.drop(
        labels=["IID",],
        axis="columns",
        inplace=True
    )
    # Translate column names.
    translations = dict()
    translations["FID"] = "identifier_ukb"
    table_raw.rename(
        columns=translations,
        inplace=True,
    )
    # Aggregate scores.
    table_aggregation = organize_aggregate_metabolite_genetic_scores(
        identifier=metabolite,
        column_index="identifier_ukb",
        columns_scores=[
            "X5e.08", "X1e.07", "X1e.06", "X1e.05", "X0.0001", "X0.001",
            "X0.01", "X0.05", "X0.1", "X0.2", "X1",
        ],
        table=table_raw,
    )
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("reporting from: read_aggregate_metabolite_genetic_scores()")
        print("blah blah...")
        print(table_aggregation)
        utility.print_terminal_partition(level=3)
    # Return.
    return table_aggregation


def read_aggregate_test_metabolite_genetic_score(
    metabolite=None,
    metabolites_files_paths=None,
    report=None,
):
    """
    Reads a metabolite's genetic scores across the UK Biobank from file,
    aggregates scores by Singular Value Decomposition (SVD), and tests this
    method.

    arguments:
        metabolite (str): identifier of a metabolite
        metabolites_files_paths (dict<list<str>>): collection of files and paths
            for metabolites
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    table = read_aggregate_metabolite_genetic_scores(
        metabolite=metabolite,
        metabolites_files_paths=metabolites_files_paths,
        report=True,
    )
    table.dropna(
        axis="index",
        how="any",
        subset=["identifier_ukb"],
        inplace=True,
    )
    table.set_index(
        "identifier_ukb",
        drop=True,
        inplace=True,
    )
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("reporting from: read_aggregate_test_metabolite_genetic_score()")
        print("blah blah...")
        print(table)
        utility.print_terminal_partition(level=3)
    # Compile information.
    pail = dict()
    # Return.
    return pail




def read_aggregate_collect_metabolites_genetic_scores(
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

    # TODO: this merge collection strategy does work!

    # Initialize a table for collection.
    table_collection = pandas.DataFrame(columns=["identifier_ukb"])
    # UK Biobank identifier is in column "FID" within the metabolite tables
    # rename to "identifier_ukb"
    table_collection.set_index(
        "identifier_ukb",
        drop=True,
        inplace=True,
    )

    for metabolite in metabolites_files_paths.keys():
        # TODO: function should return a table with UKB ID and metabolite scores
        if False:
            table_metabolite = pandas.DataFrame({
                "identifier_ukb": ["a", "b", "c", "d", "e"],
                metabolite: [1, 2, 3, 4, 5],
            })
        table_metabolite = read_aggregate_metabolite_genetic_scores(
            metabolite=metabolite,
            metabolites_files_paths=metabolites_files_paths,
            report=False,
        )
        table_metabolite.dropna(
            axis="index",
            how="any",
            subset=["identifier_ukb"],
            inplace=True,
        )
        table_metabolite.set_index(
            "identifier_ukb",
            drop=True,
            inplace=True,
        )
        table_collection = table_collection.merge(
            table_metabolite,
            how="outer",
            left_on="identifier_ukb",
            right_on="identifier_ukb",
            suffixes=("_original", "_novel"),
        )
        # TODO: here, I should merge that table into the collection table...

        pass

    # Compile information.
    #pail = dict()
    # Return information.
    return table_collection






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

    # Test the aggregation method for a single metabolite.
    # M00599: pyruvate
    # M32315: serine
    # M02342: serotonin
    # M00054: tryptophan
    pail_test = read_aggregate_test_metabolite_genetic_scores(
        metabolite="M00054",
        metabolites_files_paths=source["metabolites_files_paths"],
        report=True,
    )

    if False:
        # Collect metabolites' genetic scores, and aggregate these by singular value
        # decomposition (SVD).
        # pail_metabolites_scores
        table_scores = read_aggregate_collect_metabolites_genetic_scores(
            metabolites_files_paths=source["metabolites_files_paths"],
        )
        print(table_scores)

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

    if False:
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
