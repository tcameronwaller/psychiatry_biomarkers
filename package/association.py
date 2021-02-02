
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


def select_columns_merge_metabolite_phenotype_tables(
    phenotype=None,
    metabolite=None,
    covariates=None,
    table_metabolites_scores=None,
    table_phenotypes=None,
    report=None,
):
    """
    Merges and organizes information from metabolite and phenotype tables.

    arguments:
        phenotype (str): name of column in phenotype table for variable to set
            as dependent variable in regressions
        metabolite (str): identifiers of a single metabolite for which to
            regress genetic scores against phenotypes across UK Biobank
        covariates (list<str>): names of columns in phenotype table for
            variables to set as covariate independent variables in regressions
        table_metabolites_scores (object): Pandas data frame of metabolites'
            genetic scores across UK Biobank cohort
        table_phenotypes (object): Pandas data frame of phenotype variables
            across UK Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of dependent and independent variables for
            regression

    """

    # Copy information.
    table_metabolites_scores = table_metabolites_scores.copy(deep=True)
    table_phenotypes = table_phenotypes.copy(deep=True)
    # Select relevant columns from metabolites table.
    columns_metabolites = list()
    columns_metabolites.append("identifier_ukb")
    columns_metabolites.append(metabolite)
    table_metabolites_scores.reset_index(
        level=None,
        inplace=True
    )
    table_metabolites_scores = table_metabolites_scores.loc[
        :, table_metabolites_scores.columns.isin(columns_metabolites)
    ]
    # Select relevant columns from phenotypes table.
    columns_phenotypes = list()
    columns_phenotypes.extend(["eid", "IID"])
    columns_phenotypes.append(phenotype)
    columns_phenotypes.extend(covariates)
    table_phenotypes.reset_index(
        level=None,
        inplace=True
    )
    table_phenotypes = table_phenotypes.loc[
        :, table_phenotypes.columns.isin(columns_phenotypes)
    ]
    # Drop any rows with null keys.
    table_metabolites_scores.dropna(
        axis="index",
        how="any",
        subset=["identifier_ukb"],
        inplace=True,
    )
    table_phenotypes.dropna(
        axis="index",
        how="any",
        subset=["IID"],
        inplace=True,
    )
    # Set keys as indices.
    table_metabolites_scores["identifier_ukb"].astype("string")
    table_metabolites_scores.set_index(
        "identifier_ukb",
        drop=True,
        inplace=True,
    )
    table_phenotypes["IID"].astype("string")
    table_phenotypes.set_index(
        "IID",
        drop=True,
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report source: select_columns_merge_metabolite_phenotype_tables()"
        )
        utility.print_terminal_partition(level=3)
        print("Tables before merge:")
        print(table_phenotypes)
        print(table_metabolites_scores)
        pass

    # Merge tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_merge = table_phenotypes.merge(
        table_metabolites_scores,
        how="outer",
        left_on="IID",
        left_index=True,
        right_on="identifier_ukb",
        right_index=True,
        suffixes=("_phenotypes", "_metabolites"),
    )
    # Organize new index.
    table_merge.reset_index(
        level=None,
        inplace=True
    )
    table_merge.drop(
        labels=["IID", "identifier_ukb",],
        axis="columns",
        inplace=True
    )
    table_merge.dropna(
        axis="index",
        how="any",
        subset=["eid"],
        inplace=True,
    )
    table_merge["eid"].astype("string")
    table_merge.set_index(
        "eid",
        drop=True,
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report source: select_columns_merge_metabolite_phenotype_tables()"
        )
        utility.print_terminal_partition(level=3)
        print("Table after merge:")
        print(table_merge)
        pass
    # Return information.
    return table_merge


def remove_null_records_standardize_variables_scales(
    table=None,
    report=None,
):
    """
    Removes records with null values and standardizes variables' scales.

    arguments:
        table (object): Pandas data frame of dependent and independent variables
            for regression
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of dependent and independent variables for
            regression

    """

    # Copy information.
    table = table.copy(deep=True)
    # Drop any rows with null keys.
    table.dropna(
        axis="index",
        how="any",
        subset=None,
        inplace=True,
    )
    # Standardize variables' scales.
    table_scale = utility.standardize_table_values_by_column(
        table=table,
        report=report,
    )
    # Return information.
    return table_scale


def organize_dependent_independent_variables_table(
    phenotype=None,
    metabolite=None,
    covariates=None,
    table_metabolites_scores=None,
    table_phenotypes=None,
    report=None,
):
    """
    Organizes information and regresses metabolites' genetic scores against
    phenotypes across the UK Biobank.

    Metabolites' genetic scores match to persons in the UK Biobank by the "IID"
    for their genotypic information.

    arguments:
        phenotype (str): name of column in phenotype table for variable to set
            as dependent variable in regressions
        metabolite (str): identifiers of a single metabolite for which to
            regress genetic scores against phenotypes across UK Biobank
        covariates (list<str>): names of columns in phenotype table for
            variables to set as covariate independent variables in regressions
        table_metabolites_scores (object): Pandas data frame of metabolites'
            genetic scores across UK Biobank cohort
        table_phenotypes (object): Pandas data frame of phenotype variables
            across UK Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of dependent and independent variables for
            regression

    """

    # Select relevant columns and merge tables.
    table_merge = select_columns_merge_metabolite_phenotype_tables(
        phenotype=phenotype,
        metabolite=metabolite,
        covariates=covariates,
        table_metabolites_scores=table_metabolites_scores,
        table_phenotypes=table_phenotypes,
        report=report,
    )
    # Drop records with null values.
    # Standardize scale of variables.
    table_scale = remove_null_records_standardize_variables_scales(
        table=table_merge,
        report=report,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Report source: organize_dependent_independent_variables_table()")
        print(table_scale)
        pass
    # Return.
    return table_scale


def organize_regress_metabolite_genetic_scores_against_phenotypes(
    phenotype=None,
    metabolite=None,
    covariates=None,
    table_metabolites_scores=None,
    table_phenotypes=None,
    report=None,
):
    """
    Organizes information and regresses metabolites' genetic scores against
    phenotypes across the UK Biobank.

    Metabolites' genetic scores match to persons in the UK Biobank by the "IID"
    for their genotypic information.

    arguments:
        phenotype (str): name of column in phenotype table for variable to set
            as dependent variable in regressions
        metabolite (str): identifiers of a single metabolite for which to
            regress genetic scores against phenotypes across UK Biobank
        covariates (list<str>): names of columns in phenotype table for
            variables to set as covariate independent variables in regressions
        table_metabolites_scores (object): Pandas data frame of metabolites'
            genetic scores across UK Biobank cohort
        table_phenotypes (object): Pandas data frame of phenotype variables
            across UK Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information from regression

    """

    # Organize information for regression.
    table_organization = organize_dependent_independent_variables_table(
        phenotype=phenotype,
        metabolite=metabolite,
        covariates=covariates,
        table_metabolites_scores=table_metabolites_scores,
        table_phenotypes=table_phenotypes,
        report=report,
    )

    # TODO: now call a function to perform the actual regression...
    # TODO: this function should return NAs if the data are not adequate... set a threshold...

    # TODO: report relevant information... but also prepare dict for return

    # Compile information.
    pail = dict()
    # Return information.
    return pail


def organize_regress_metabolites_genetic_scores_against_phenotypes(
    phenotype=None,
    metabolites=None,
    covariates=None,
    table_metabolites_scores=None,
    table_phenotypes=None,
    report=None,
):
    """
    Organizes information and regresses metabolites' genetic scores against
    phenotypes across the UK Biobank.

    Metabolites' genetic scores match to persons in the UK Biobank by the "IID"
    for their genotypic information.

    arguments:
        phenotype (str): name of column in phenotype table for variable to set
            as dependent variable in regressions
        metabolites (list<str>): identifiers of metabolites for which to regress
            their genetic scores against phenotypes across UK Biobank
        covariates (list<str>): names of columns in phenotype table for
            variables to set as covariate independent variables in regressions
        table_metabolites_scores (object): Pandas data frame of metabolites'
            genetic scores across UK Biobank cohort
        table_phenotypes (object): Pandas data frame of phenotype variables
            across UK Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): information from regressions

    """

    # Collect information for metabolites.
    records = list()
    for metabolite in metabolites:
        record = (
            organize_regress_metabolite_genetic_scores_against_phenotypes(
                phenotype="body_mass_index",
                metabolite=metabolite,
                covariates=covariates,
                table_metabolites_scores=table_metabolites_scores,
                table_phenotypes=table_phenotypes,
                report=report,
        ))
        records.append(record)
        pass
    # Organize information in a table.
    table_regression = utility.convert_records_to_dataframe(
        records=records
    )
    print(table_regression)

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
    print("version check: 3")
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
            phenotype="body_mass_index", # "testosterone"
            metabolites=["M00599", "M32315", "M02342", "M00054"],
            covariates=[
                "sex", "age",
                "genotype_pc_1", "genotype_pc_2", "genotype_pc_3",
                "genotype_pc_4", "genotype_pc_5", "genotype_pc_6",
                "genotype_pc_7", "genotype_pc_8", "genotype_pc_9",
                "genotype_pc_10",
            ],
            table_metabolites_scores=source["table_metabolites"],
            table_phenotypes=source["table_phenotypes"],
            report=True,
    ))


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
