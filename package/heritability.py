
"""
...

This module collects and organizes information about heritability estimates for
metabolites.

Source of GWAS summary statistics is Shin et al, Nature Genetics, 2014
(PubMed:24816252). Metabolite identifiers are from Metabolon Inc.

Method for estimation of heritability from GWAS summary statistics was linkage
disequilibrium score (LDSC) regression in tool LD SCore
(https://github.com/bulik/ldsc).

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Import modules from specific path without having to install a general package
# I would have to figure out how to pass a path variable...
# https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path


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
    paths["heritability"] = os.path.join(path_dock, "heritability")
    paths["collection"] = os.path.join(
        path_dock, "heritability", "collection"
    )
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["collection"])
    # Initialize directories.
    utility.create_directories(
        path=paths["collection"]
    )
    # Return information.
    return paths


def read_source(
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_table_metabolite_identifiers = os.path.join(
        path_dock, "access", "metabolites", "metaboliteMap.txt"
    )
    path_metabolite_files = os.path.join(
        path_dock, "access", "metabolites", "metabolite_files.txt"
    )
    # Read information from file.
    table_metabolite_identifiers = pandas.read_csv(
        path_table_metabolite_identifiers,
        sep="\t",
        header=0,
        dtype="string",
    )
    metabolite_files = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_metabolite_files,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_metabolite_identifiers)
        utility.print_terminal_partition(level=2)
        print("count of metabolites: " + str(len(metabolite_files)))
        print("example: " + str(metabolite_files[10]))
        utility.print_terminal_partition(level=2)
    # Compile and return information.
    return {
        "table_metabolite_identifiers": table_metabolite_identifiers,
        "metabolite_files": metabolite_files,
    }


def extract_metabolite_file_identifiers(
    metabolite_files=None,
):
    """
    Extracts metabolite identifiers from names of files.

    arguments:
        metabolite_files (list<str>): names of files with metabolite GWAS
            summary statistics

    raises:

    returns:
        (list<str>): identifiers of metabolites

    """

    # Iterate across metabolite file names.
    # Collect metabolite identifiers.
    identifiers = list()
    for file in metabolite_files:
        components = file.split(".")
        identifier = str(components[0]).strip()
        if ((len(identifier) > 1) and (identifier[0] == "M")):
            identifiers.append(identifier)
            pass
        pass
    # Make sure that all identifiers are unique.
    identifiers_unique = utility.collect_unique_elements(
        elements=identifiers,
    )
    # Return information.
    return identifiers_unique


def read_extract_metabolite_heritability(
    metabolite_identifier=None,
    path_dock=None,
):
    """
    Reads and extracts information from log of LDSC.

    arguments:
        metabolite_identifier (str): identifier of a metabolite
        path_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:
        (dict): information about estimation of a metabolite's heritability

    """

    # Define path to file.
    name_file = str(metabolite_identifier + "_heritability.log")
    path_heritability = os.path.join(
        path_dock, "heritability", "metabolites", name_file
    )
    # Read relevant lines from file.
    lines = utility.read_file_text_lines(
        path_file=path_file,
        start=20,
        stop=25,
    )
    # Extract information from lines.
    prefix_variants = "After merging with regression SNP LD, "
    suffix_variants = " SNPs remain."
    prefix_heritability = "Total Observed scale h2: "
    suffix_heritability = " ("
    suffix_error = ")"
    variants = float("nan")
    heritability = float("nan")
    standard_error = float("nan")
    for line in lines:
        if prefix_variants in line:
            variants = float(
                line.replace(prefix_variants, "").replace(suffix_variants, "")
            )
        elif prefix_heritability in line:
            content = line.replace(prefix_heritability, "")
            contents = content.split(" (")
            heritability = float(contents[0])
            standard_error = float(
                contents[1].replace(")", "")
            )
            pass
        pass
    # Collect information.
    record = dict()
    record["identifier"] = metabolite_identifier
    record["variants"] = variants
    record["heritability"] = heritability
    record["standard_error"] = standard_error
    # Return information.
    return record


def read_collect_metabolite_heritabilities(
    metabolite_identifiers=None,
    path_dock=None,
    report=None,
):
    """
    Reads, collects, and organizes metabolite heritability estimates.

    arguments:
        metabolite_identifiers (list<str>): identifiers of metabolites
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of metabolite heritability estimates

    """

    records = list()
    for identifier in metabolite_identifiers:
        record = read_extract_metabolite_heritability(
            metabolite_identifier=metabolite_identifier,
            path_dock=path_dock,
        )
        records.append(record)
        pass
    # Organize table.
    table = utility.convert_records_to_dataframe(
        records=records
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table)
    # Return information.
    return table








def extract_organize_variable_fields_instances_names(
    data_ukbiobank_variables=None,
    extra_names=None,
):
    """
    Organizes column names for variable fields and instances.

    arguments:
        data_ukbiobank_variables (object): Pandas data frame of information
            about UK Biobank phenotype variables
        extra_names (list<str>): extra names to include

    raises:

    returns:
        (list<str>): column names for variable fields and instances

    """

    # Copy data.
    data_variables = data_ukbiobank_variables.copy(deep=True)
    # Organize information.
    data_variables = data_variables.loc[
        :, data_variables.columns.isin(["field", "instances_keep"])
    ]
    records = utility.convert_dataframe_to_records(data=data_variables)
    # Iterate across records for rows.
    # Collect variables' names and types.
    names = list()
    names.extend(extra_names)
    for record in records:
        field = str(record["field"])
        instances_total_raw = str(record["instances_keep"])
        instances_raw = instances_total_raw.split(",")
        for instance_raw in instances_raw:
            instance = str(instance_raw).strip()
            if len(instance) > 0:
                name = str(field + "-" + instance)
                # Create entry for variable's name and type.
                names.append(name)
                pass
            pass
        pass
    # Return information.
    return names


def remove_data_irrelevant_variable_instances_entries(
    data_ukbiobank_variables=None,
    data_ukb_41826=None,
    data_ukb_43878=None,
    report=None,
):
    """
    Removes irrelevant columns and rows from data.

    arguments:
        data_ukbiobank_variables (object): Pandas data frame of information
            about UK Biobank phenotype variables
        data_ukb_41826 (object): Pandas data frame of variables from UK Biobank
            phenotype accession 41826
        data_ukb_43878 (object): Pandas data frame of variables from UK Biobank
            phenotype accession 43878
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data frames after removal of
            irrelevant columns and rows

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("...before pruning...")
        print("data_ukb_41826 shape: " + str(data_ukb_41826.shape))
        utility.print_terminal_partition(level=4)
        print("data_ukb_43878 shape: " + str(data_ukb_43878.shape))

    # Extract names of columns for relevant variable fields and instances.
    column_names = extract_organize_variable_fields_instances_names(
        data_ukbiobank_variables=data_ukbiobank_variables,
        extra_names=["IID", "eid"],
    )
    print(column_names)
    # Remove all irrelevant columns.
    data_ukb_41826 = data_ukb_41826.loc[
        :, data_ukb_41826.columns.isin(column_names)
    ]
    data_ukb_43878 = data_ukb_43878.loc[
        :, data_ukb_43878.columns.isin(column_names)
    ]
    # Remove rows with all missing values.
    data_ukb_41826.dropna(
        axis="index",
        how="all",
        inplace=True,
    )
    data_ukb_41826.dropna(
        axis="index",
        how="any",
        subset=["eid"],
        inplace=True,
    )
    data_ukb_43878.dropna(
        axis="index",
        how="all",
        inplace=True,
    )
    data_ukb_43878.dropna(
        axis="index",
        how="any",
        subset=["eid"],
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("...after pruning...")
        print("data_ukb_41826 shape: " + str(data_ukb_41826.shape))
        utility.print_terminal_partition(level=4)
        print("data_ukb_43878 shape: " + str(data_ukb_43878.shape))

    # Compile and return information.
    bucket = dict()
    bucket["data_ukb_41826"] = data_ukb_41826
    bucket["data_ukb_43878"] = data_ukb_43878
    return bucket


def merge_data_variables_identifiers(
    data_identifier_pairs=None,
    data_ukb_41826=None,
    data_ukb_43878=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        data_identifier_pairs (object): Pandas data frame of associations
            between "IID" and "eid"
        data_ukb_41826 (object): Pandas data frame of variables from UK Biobank
            phenotype accession 41826
        data_ukb_43878 (object): Pandas data frame of variables from UK Biobank
            phenotype accession 43878
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(data_identifier_pairs)
        utility.print_terminal_partition(level=2)
        print(data_ukb_41826)
        utility.print_terminal_partition(level=2)
        print(data_ukb_43878)
    # Organize data.
    data_identifier_pairs.astype("string")
    data_identifier_pairs.set_index(
        "eid",
        drop=True,
        inplace=True,
    )
    data_ukb_41826["eid"].astype("string")
    data_ukb_41826.set_index(
        "eid",
        drop=True,
        inplace=True,
    )
    data_ukb_43878["eid"].astype("string")
    data_ukb_43878.set_index(
        "eid",
        drop=True,
        inplace=True,
    )

    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    data_merge = data_identifier_pairs.merge(
        data_ukb_41826,
        how="outer",
        left_on="eid",
        right_on="eid",
        suffixes=("_pairs", "_41826"),
    )
    data_merge = data_merge.merge(
        data_ukb_43878,
        how="outer",
        left_on="eid",
        right_on="eid",
        suffixes=("_41826", "_43878"),
    )
    # Remove excess columns.

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(data_merge)
    # Return information.
    return data_merge


def exclude_persons_ukbiobank_consent(
    exclusion_identifiers=None,
    data=None,
    report=None,
):
    """
    Removes entries with specific identifiers from data.

    arguments:
        exclusion_identifiers (list<str>): identifiers of persons who withdrew
            consent from UK Biobank
        data (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Initial data dimensions: " + str(data.shape))
        utility.print_terminal_partition(level=4)
        print("Exclusion of persons: " + str(len(exclusion_identifiers)))
    # Copy data.
    data = data.copy(deep=True)
    # Filter data entries.
    data_exclusion = data.loc[
        ~data.index.isin(exclusion_identifiers), :
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        print("Final data dimensions: " + str(data_exclusion.shape))
    # Return information.
    return data_exclusion


def convert_data_variable_types(
    data=None,
    report=None,
):
    """
    Converts data variable types.

    The UK Biobank encodes several nominal variables with integers. Missing
    values for these variables necessitates some attention to type conversion
    from string to float for analysis.

    arguments:
        data (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Before type conversion")
        utility.print_terminal_partition(level=3)
        print(data.dtypes)
    # Copy data.
    data = data.copy(deep=True)
    # Convert data variable types.
    data["20117-0.0"] = pandas.to_numeric(
        data["20117-0.0"],
        errors="coerce", # force any invalid values to missing or null
        downcast="float",
    )
    if False:
        data["20117-0.0"].fillna(
            value="-3",
            axis="index",
            inplace=True,
        )
        data["20117-0.0"].astype(
            "float32",
            copy=True,
        )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("After type conversion")
        utility.print_terminal_partition(level=3)
        print(data.dtypes)
        utility.print_terminal_partition(level=4)
        print(data["20117-0.0"].value_counts())
    # Return information.
    return data


def calculate_sum_drinks(
    beer_cider=None,
    wine_red=None,
    wine_white=None,
    port=None,
    liquor=None,
    other=None,
):
    """
    Calculates the sum of drinks.

    UK Biobank data coding 100291.
    "do not know": -1
    "prefer not to answer": -3

    arguments:
        beer_cider (int): count of drinks by type
        wine_red (int): count of drinks by type
        wine_white (int): count of drinks by type
        port (int): count of drinks by type
        liquor (int): count of drinks by type
        other (int): count of drinks by type
        report (bool): whether to print reports

    raises:

    returns:
        (int): sum of counts of drinks

    """

    # Consider all types of alcoholic beverage.
    types = [beer_cider, wine_red, wine_white, port, liquor, other]
    # Determine whether all relevant variables have missing values.
    valid = False
    for type in types:
        if ((type != -1) and (type != -3)):
            valid = True
            pass
        pass
    if valid:
        # Variables have valid values.
        # Calculate sum of drinks.
        drinks = 0
        for type in types:
            if ((not pandas.isna(type)) and (type >= 0)):
                drinks = drinks + type
                pass
            pass
    else:
        drinks = float("nan")
    return drinks


def calculate_total_alcohol_consumption_monthly(
    drinks_weekly=None,
    drinks_monthly=None,
    weeks_per_month=None,
):
    """
    Calculate monthly alcohol consumption in drinks.

    arguments:
        drinks_weekly (float): sum of weekly drinks from weekly variables
        drinks_monthly (float): sum of monthly drinks from monthly variables
        weeks_per_month (float): factor to use for weeks per month

    raises:

    returns:
        (float): person's monthly alcohol consumption in drinks

    """

    # Use as much valid information as is available.
    if (not math.isnan(drinks_weekly) and not math.isnan(drinks_monthly)):
        alcohol_drinks_monthly = (
            drinks_monthly + (weeks_per_month * drinks_weekly)
        )
    elif (not math.isnan(drinks_weekly)):
        alcohol_drinks_monthly = (weeks_per_month * drinks_weekly)
    elif (not math.isnan(drinks_monthly)):
        alcohol_drinks_monthly = drinks_monthly
    else:
        # There is no valid information about alcohol consumption quantity.
        alcohol_drinks_monthly = float("nan")
        pass
    return alcohol_drinks_monthly


def determine_total_alcohol_consumption_monthly(
    status=None,
    drinks_weekly=None,
    drinks_monthly=None,
    weeks_per_month=None,
):
    """
    Calculate monthly alcohol consumption in drinks.

    UK Biobank data coding 90.
    "current": 2
    "previous": 1
    "never": 0
    "prefer not to answer": -3

    Accommodate inexact float values.

    arguments:
        status (str): person's status of alcohol consumption
        drinks_weekly (float): sum of weekly drinks from weekly variables
        drinks_monthly (float): sum of monthly drinks from monthly variables
        weeks_per_month (float): factor to use for weeks per month

    raises:

    returns:
        (float): person's monthly alcohol consumption in drinks

    """

    # Calculate alcohol consumption quantity.
    alcohol_monthly = calculate_total_alcohol_consumption_monthly(
        drinks_weekly=drinks_weekly,
        drinks_monthly=drinks_monthly,
        weeks_per_month=weeks_per_month,
    )
    # Consider alcohol consumption status.
    if (not pandas.isna(status)):
        if (-0.5 < status and status < 0.5):
            # Confirm that alcohol consumption is none.
            if (not math.isnan(alcohol_monthly)):
                alcohol_drinks_monthly = alcohol_monthly
            else:
                alcohol_drinks_monthly = 0.0
            pass
        elif (
            (1.5 < status and status < 2.5) or
            (0.5 < status and status < 1.5) or
            (-3.5 < status and status < -2.5)
        ):
            # Determine alcohol consumption quantity.
            alcohol_drinks_monthly = alcohol_monthly
            pass
        else:
            alcohol_drinks_monthly = alcohol_monthly
    else:
        alcohol_drinks_monthly = alcohol_monthly
        pass
    # Return information.
    return alcohol_drinks_monthly


def organize_alcohol_consumption_monthly_drinks(
    data=None,
    report=None,
):
    """
    Organizes information about alcohol consumption in standard drinks monthly.

    arguments:
        data (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about quantity of alcohol consumption

    """

    # Copy data.
    data = data.copy(deep=True)
    # Calculate sum of drinks weekly.
    data["drinks_weekly"] = data.apply(
        lambda row:
            calculate_sum_drinks(
                beer_cider=row["1588-0.0"],
                wine_red=row["1568-0.0"],
                wine_white=row["1578-0.0"],
                port=row["1608-0.0"],
                liquor=row["1598-0.0"],
                other=row["5364-0.0"],
            ),
        axis="columns", # apply across rows
    )
    # Calculate sum of drinks monthly.
    data["drinks_monthly"] = data.apply(
        lambda row:
            calculate_sum_drinks(
                beer_cider=row["4429-0.0"],
                wine_red=row["4407-0.0"],
                wine_white=row["4418-0.0"],
                port=row["4451-0.0"],
                liquor=row["4440-0.0"],
                other=row["4462-0.0"],
            ),
        axis="columns", # apply across rows
    )
    # Determine sum of total drinks monthly.
    data["alcohol_drinks_monthly"] = data.apply(
        lambda row:
            determine_total_alcohol_consumption_monthly(
                status=row["20117-0.0"],
                drinks_weekly=row["drinks_weekly"],
                drinks_monthly=row["drinks_monthly"],
                weeks_per_month=4.345, # 52.143 weeks per year (12 months)
            ),
        axis="columns", # apply across rows
    )
    # Remove columns for variables that are not necessary anymore.
    data_clean = data.copy(deep=True)
    data_clean.drop(
        labels=[
            "1588-0.0", "1568-0.0", "1578-0.0", "1608-0.0", "1598-0.0",
            "5364-0.0",
            "4429-0.0", "4407-0.0", "4418-0.0", "4451-0.0", "4440-0.0",
            "4462-0.0",
            "20117-0.0",
            "drinks_weekly",
            "drinks_monthly",
        ],
        axis="columns",
        inplace=True
    )
    # Organize data for report.
    data_report = data.copy(deep=True)
    data_report = data_report.loc[
        :, data_report.columns.isin([
            "eid", "IID",
            "1588-0.0", "1568-0.0", "1578-0.0", "1608-0.0", "1598-0.0",
            "5364-0.0",
            "4429-0.0", "4407-0.0", "4418-0.0", "4451-0.0", "4440-0.0",
            "4462-0.0",
            "20117-0.0",
            "drinks_weekly",
            "drinks_monthly",
            "alcohol_drinks_monthly",
        ])
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary of alcohol consumption quantity variables: ")
        print(data_report)
    # Collect information.
    bucket = dict()
    bucket["data"] = data
    bucket["data_clean"] = data_clean
    bucket["data_report"] = data_report
    # Return information.
    return bucket



# TODO: organize AUDIT-C variables... sum of AUDIT-C questions 1-3
# TODO: handle the specific missing value codes




# TODO: derive "hysterectomy" and "oophrectomy" separately
# TODO: then derive "menopause"

def organize_female_menopause(
    data=None,
    report=None,
):
    """
    Organizes information about whether female persons have experienced
    menopause.

    arguments:
        data (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy data.
    data = data.copy(deep=True)
    # Calculate sum of drinks weekly.
    data["menopause"] = data.apply(
        lambda row:
            calculate_sum_drinks(
                beer_cider=row["1588-0.0"],
                wine_red=row["1568-0.0"],
                wine_white=row["1578-0.0"],
                port=row["1608-0.0"],
                liquor=row["1598-0.0"],
                other=row["5364-0.0"],
            ),
        axis="columns", # apply across rows
    )
    data["hysterectomy"] = data.apply(
        lambda row:
            calculate_sum_drinks(
                beer_cider=row["1588-0.0"],
                wine_red=row["1568-0.0"],
                wine_white=row["1578-0.0"],
                port=row["1608-0.0"],
                liquor=row["1598-0.0"],
                other=row["5364-0.0"],
            ),
        axis="columns", # apply across rows
    )
    data["oophorectomy"] = data.apply(
        lambda row:
            calculate_sum_drinks(
                beer_cider=row["1588-0.0"],
                wine_red=row["1568-0.0"],
                wine_white=row["1578-0.0"],
                port=row["1608-0.0"],
                liquor=row["1598-0.0"],
                other=row["5364-0.0"],
            ),
        axis="columns", # apply across rows
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        print("Final data dimensions: " + str(data_exclusion.shape))
    # Return information.
    return data_exclusion



def write_product_alcohol_consumption(
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
    path_data_report_quantity = os.path.join(
        path_parent, "data_report_quantity.tsv"
    )
    # Write information to file.
    information["data_report"].to_csv(
        path_or_buf=path_data_report_quantity,
        sep="\t",
        header=True,
        index=True,
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

    # Alcohol consumption.
    write_product_alcohol_consumption(
        information=information["alcohol_consumption"],
        path_parent=paths["alcohol_consumption"],
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
    source = read_source(
        path_dock=path_dock,
        report=True,
    )
    # Extract identifiers of metabolites with GWAS summary statistics.
    metabolite_identifiers = extract_metabolite_file_identifiers(
        metabolite_files=source["metabolite_files"],
    )
    # Collect linkage disequilibrium score regression heritability estimates
    # for each metabolite.
    table_metabolite_heritability = read_collect_metabolite_heritabilities(
        metabolite_identifiers=metabolite_identifiers,
        path_dock=path_dock,
        report=True,
    )


    # Merge heritability table with metabolite name table.

    # 1. read identifiers of all metabolites from list file
    # 2. check that LD Score regression files exist for all metabolites
    # 3. iterate on metabolites with LD Score regression files
    # 4. for each metabolite LDSC result file
    # 5. --- find relevant rows by leading strings
    # 6. --- extract and organize
    # 7. organize summary table
    # 8. save summary table

    if False:
        # Remove data columns for irrelevant variable instances.
        prune = remove_data_irrelevant_variable_instances_entries(
            data_ukbiobank_variables=source["data_ukbiobank_variables"],
            data_ukb_41826=source["data_ukb_41826"],
            data_ukb_43878=source["data_ukb_43878"],
            report=True,
        )

        # Collect information.
        information = dict()
        information["alcohol_consumption"] = dict()
        information["alcohol_consumption"]["data_report"] = (
            bin_consumption["data_report"]
        )
        # Write product information to file.
        write_product(
            paths=paths,
            information=information
        )


    pass



if (__name__ == "__main__"):
    execute_procedure()
