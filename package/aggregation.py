
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
import pandas
pandas.options.mode.chained_assignment = None # default = "warn"

import scipy.stats
import scipy.linalg
import statsmodels.multivariate.pca

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
    paths["selection"] = os.path.join(path_dock, "aggregation", "selection")
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["aggregation"])
    # Initialize directories.
    utility.create_directories(
        path=paths["aggregation"]
    )
    utility.create_directories(
        path=paths["selection"]
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
        metabolites = ["M00599", "M32315", "M02342", "M00054"]
        match = any(list(map(
            lambda metabolite: (metabolite in path_file), metabolites
        )))
        if match:
            utility.print_terminal_partition(level=2)
            print("Path: " + str(path_file))
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


################################################################################
#####
###
##
#

##########
# SVD sign adjustment

# TODO: enable this function to adjust for a specified count of singular vectors
# TODO: (TCW 16 August 2021) I think this function is incomplete...
# TODO: (TCW 17 August 2021) DO NOT use this unless VERY careful about what it's doing
def adjust_singular_value_decomposition_factor_signs(
    matrix=None,
    singular_values=None,
    left_singular_vectors_columns=None,
    right_singular_vectors_rows=None,
    report=None,
):
    """
    Adjusts the otherwise random signs of factors from Singular Value
    Decomposition (SVD) to reduce the directional ambiguity.

    Reference:

    arguments:
        matrix (object): NumPy array matrix of original values across samples
            (rows, dimension 0) and variables (columns, dimension 1)
        singular_values (object): NumPy array of Singular Values from SVD
        left_singular_vectors_columns (object): NumPy array matrix with SVD left
            singular vectors as columns, U
        right_singular_vectors_rows (object): NumPy array matrix with SVD right
            singular vectors as rows, VT or Vh
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the singular value
            decomposition

    """

    # Copy information.
    matrix = numpy.copy(matrix)
    s = numpy.copy(singular_values)
    u = numpy.copy(left_singular_vectors_columns)
    vt = numpy.copy(right_singular_vectors_rows)
    # Organize information.
    matrix_transpose = numpy.transpose(matrix)
    s_diagonal = numpy.diag(s)
    ut = numpy.copy(numpy.transpose(u))
    v = numpy.copy(numpy.transpose(vt))
    # Calculate basic products by matrix multiplication.
    ut_y = numpy.dot(ut, matrix)
    vt_y = numpy.dot(vt, matrix_transpose)
    # Reduce values to indicators of positive and negative signs.
    ut_y_sign = numpy.sign(ut_y)
    vt_y_sign = numpy.sign(vt_y)
    # Calculate squares of matrices.
    # Calculation of square by matrix multiplifcation is only possible for
    # square matrices.
    # Instead calculate the squares of all individual values in the matrices.
    ut_y_square = numpy.square(ut_y)
    vt_y_square = numpy.square(vt_y)
    # Calculate left and right sign matrices.
    signs_left = numpy.dot(ut_y_sign, ut_y_square)
    signs_right = numpy.dot(vt_y_sign, vt_y_square)


    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "adjust_singular_value_decomposition_factor_signs()"
        )
        utility.print_terminal_partition(level=2)
        print("Shape of original matrix: " + str(matrix.shape))
        print("rows (dimension 0): samples (cases, observations)")
        print("columns (dimension 1): variables (features)")
        utility.print_terminal_partition(level=4)
        print("Shape of matrix Sigma (singular values): " + str(s.shape))
        utility.print_terminal_partition(level=4)
        print("Shape of matrix U (left singular vectors): " + str(u.shape))
        print(
            "Shape of matrix UT (transpose left singular vectors): " +
            str(ut.shape)
        )
        utility.print_terminal_partition(level=4)
        print(
            "Shape of matrix VT (transpose right singular vectors): " +
            str(vt.shape)
        )
        print("Shape of matrix V (right singular vectors): " + str(v.shape))
        utility.print_terminal_partition(level=4)
        print("Shape of matrix UT-Y (product): " + str(ut_y.shape))
        print("Shape of matrix UT-Y square: " + str(ut_y_square.shape))
        print("Shape of matrix VT-Y (product): " + str(vt_y.shape))
        print("Shape of matrix VT-Y square: " + str(vt_y_square.shape))
        utility.print_terminal_partition(level=4)
        print("Shape of left signs matrix: " + str(signs_left.shape))
        print("Shape of right signs matrix: " + str(signs_right.shape))
        pass
    # Compile information.
    pail = dict()
    pail["matrix"] = matrix
    pail["left_singular_vectors_columns"] = u_prime
    pail["singular_values"] = s
    pail["right_singular_vectors_rows"] = vt_prime
    # Return.
    return pail


##########
# Principal Component Factors


# TODO: this is the function to derive the PCA Eigenvalues, Eigenvectors, etc from adjusted SVD factors
def calculate_principal_components_from_singular_value_decomposition(
    singular_values=None,
    left_singular_vectors=None,
    right_singular_vectors=None,
    which_singular_vectors=None,
    table=None,
    report=None,
):
    """
    Calculates Principal Components and relevant information from raw factors of
    a Singular Value Decomposition (SVD).

    Reference:
    "https://stats.stackexchange.com/questions/134282/
    relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca"

    arguments:
        threshold_valid_proportion_per_column (float): proportion of rows that
            must have a valid value for a column in order to keep the column
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index, after final scaling and filtering for SVD
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the singular value
            decomposition

    """
    # https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
    # https://towardsdatascience.com/pca-and-svd-explained-with-numpy-5d13b0d2a4d8
    # https://towardsdatascience.com/singular-value-decomposition-and-its-applications-in-principal-component-analysis-5b7a5f08d0bd
    # http://www.math.ucsd.edu/~gptesler/283/slides/pca_18-handout.pdf
    # https://www.cc.gatech.edu/~lsong/teaching/CX4240spring16/pca_wall.pdf

    # Calculate Eigenvalues.
    eigenvalues = (
        calculate_principal_component_eigenvalues_from_singular_values(
            singular_values=s,
            count_samples=pail_svd["count_samples"],
            report=False,
    ))

    # TODO: calculate eigenvectors from the factor specified in
    # "which_singular_vectors"

    # Calculate Eigenvectors.
    # Eigenvectors are the right singular vectors of the original matrix.
    eigenvectors = numpy.copy(numpy.transpose(vh))

    # TODO: I'm not entirely sure that I'm sorting the correct dimension of
    # Eigenvectors...
    # TODO: sort dimension will depend on which singular vector selected

    # Sort Eigenvectors by order of decreasing Eigenvalues.
    pail_sort = sort_eigenvectors_by_decreasing_eigenvalues(
        eigenvectors=eigenvectors,
        eigenvalues=eigenvalues,
        report=False,
    )
    # Calculate loadings.
    loadings = calculate_principal_component_loadings_from_eigen_values_vectors(
        eigenvectors=pail_sort["eigenvectors"],
        eigenvalues=pail_sort["eigenvalues"],
        report=True,
    )
    loadings_direct = (
        calculate_principal_component_loadings_from_direct_factors(
            s=s,
            vh=vh,
            count_samples=pail_organization["count_samples"],
            report=True,
    ))
    # Calculate Principal Components.
    # --> Calculate from U and S
    # or
    # --> Calculate from V and S

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_principal_components_from_singular_value_decomposition()"
        )
        utility.print_terminal_partition(level=2)

        # Original matrix has shape (M, N)
        print(
            "Shape of original matrix: " +
            str(pail_organization["matrix"].shape)
        )
        # Eigenvalues.
        print("Shape of Eigenvalues: " + str(eigenvalues.shape))
        # Eigenvectors.
        print("Shape of Eigenvectors: " + str(eigenvectors.shape))
        # Loadings.
        print("Shape of Loadings: " + str(loadings.shape))
        print(loadings)
        print(
            "Shape of Loadings from SVD factors: " +
            str(loadings_direct.shape)
        )
        print("Loadings nearly equal by both calculations: ")
        print(numpy.allclose(loadings, loadings_direct))
        pass
    # Compile information.
    pail = dict()
    #pail["table_scale"] = pail_organization["table_scale"]
    #pail["u"] = u
    #pail["singular_values"] = s
    #pail["vh"] = vh
    pail["eigenvalues"] = eigenvalues
    pail["eigenvectors"] = eigenvectors
    pail["loadings"] = loadings
    # Return.
    return pail







##########
# Driver function(s)


# TODO: this function is scrap for the PCA loadings method (simple)
def organize_principal_components_positive_sum_loadings(
    threshold_valid_proportion_per_column=None,
    table=None,
    report=None,
):
    """
    Organizes a Principal Components Analysis while forcing loadings to have a
    positive sum.

    arguments:
        threshold_valid_proportion_per_column (float): proportion of rows that
            must have a valid value for a column in order to keep the column
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the singular value
            decomposition

    """

    # Calculate original factors by Singular Value Decomposition (SVD).
    pail_original = (
        organize_principal_components_by_singular_value_decomposition(
            threshold_valid_proportion_per_column=(
                threshold_valid_proportion_per_column
            ),
            table=table,
            report=report,
    ))
    # Determine whether loadings have a positive sum.
    loadings_original = numpy.copy(pail_original["loadings"])
    sum_original = numpy.sum(loadings_original.flatten(order="C"))
    if (sum_original >= 0):
        loading_sign_flip = False
    else:
        loading_sign_flip = True
        # TODO: invert signs of all of Vh (and U???)... then re-calculate derived values.
        loadings_novel = numpy.negative(numpy.copy(pail_components.loadings))
        sum_novel = numpy.sum(loadings_novel.flatten(order="C"))


    # Organize principal component factors within table.
    table_components = organize_principal_component_factor_table(
        factors=pail_components.factors, # TODO: change to factors after sign adjustment
        prefix="component_",
        index=index,
        index_name="identifier_ukb",
        report=True,
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "organize_principal_components_positive_sum_loadings()"
        )
        utility.print_terminal_partition(level=2)
        pass
    # Compile information.
    pail = dict()
    pail["table_scale"] = pail_organization["table_scale"]
    pail["u"] = u
    pail["singular_values"] = s
    pail["vh"] = vh
    pail["eigenvalues"] = eigenvalues
    pail["eigenvectors"] = eigenvectors
    pail["loadings"] = loadings
    # Return.
    return pail

##########
# Aggregation


def organize_principal_component_factor_table(
    factors=None,
    prefix=None,
    index=None,
    index_name=None,
    report=None,
):
    """
    Organizes a table of factors or scores from Principal Component Analysis.

    arguments:
        factors (object): NumPy array matrix of Principal Component factors
        prefix (str): prefix for names of component columns
        index (object): NumPy array of indices for table
        index_name (str): name for table's index column
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    # Copy information.
    factors = numpy.copy(factors)
    # Organize information.
    count = 1
    columns = list()
    for component in range(0, pail_components.factors.shape[1], 1):
        column = str(prefix + str(count))
        columns.append(column)
        count += 1
    table = pandas.DataFrame(
        data=factors,
        index=index,
        columns=columns,
        dtype="float32",
        copy=True,
    )
    table.rename_axis(
        index=index_name,
        axis="index",
        copy=False,
        inplace=True,
    )
    table.reset_index(
        level=None,
        inplace=True
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Matrix before organization:")
        print(factors)
        utility.print_terminal_partition(level=3)
        print("Table after organization:")
        print(table)
    # Return.
    return table


def organize_principal_component_analysis_force_loadings_positive_sum(
    table=None,
    report=None,
):
    """
    Organizes aggregation by modification of Principal Components Analysis
    (PCA).

    Factorize a covariance or correlation matrix into direction (eigenvectors)
    and scale (eigenvalues).
    The eigenvalues impart scale or weight to each eigenvector.
    Each eigenvector has its own eigenvalue, and their sort orders mush match.

    loadings = eigenvectors [dot] square_root(eigenvalues)
    Loadings include aspects of both direction (eigenvectors) and scale
    (eigenvalues).

    Singular Value Decomposition assigns Eigenvector direction (sign, positive
    or negative) at random. To ensure consistency, flip all signs such that the
    sum of loadings (direction and scale) is positive.

    arguments:
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with explicit
            index
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the principal components
            analysis

    """

    # Copy information.
    table = table.copy(deep=True)
    index = copy.deepcopy(table.index)
    # Organize matrix.
    # Matrix format has samples (cases, observations) across dimension 0 (rows)
    # and variables (features) across dimension 1 (columns).
    matrix = table.to_numpy()
    # Calculate Principal Components.
    # If there is specification of "ncomp", then function only returns
    # Eigenvalues, loadings, Eigenvectors, and principal components to this
    # count.
    # Function sorts Eigenvalues in decreasing order.
    # Sort order of Eigenvectors mush match the sort order of Eigenvalues.
    # Statsmodels erroneously returns "loadings" that have identical values and
    # dimensions as the Eigenvectors.
    pail_components = statsmodels.multivariate.pca.PCA(
        matrix,
        ncomp=3, # None # temporarily reduce count to 3
        standardize=True,
        gls=False,
        weights=None,
        method="eig", # "svd", "eig", "nipals"
        missing=None, # None or "drop-row"
    )
    # Calculate loadings.
    loadings = calculate_principal_component_loadings(
        eigenvectors=pail_components.eigenvecs,
        eigenvalues=pail_components.eigenvals,
        report=True,
    )
    print("Statsmodels loadings...")
    print(pail_components.loadings)
    print("Calculated loadings...")
    print(loadings)

    if False:
        # Match sort order of Eigenvectors and Eigenvalues, by decreasing
        # Eigenvalues.
        pail_sort = sort_eigenvectors_by_decreasing_eigenvalues(
            eigenvectors=pail_components.eigenvecs,
            eigenvalues=pail_components.eigenvals,
            report=False,
        )
        # Check the sum of principal component loadings.
        loadings_original = numpy.copy(pail_components.loadings)
        sum_original = numpy.sum(loadings_original.flatten(order="C"))
        if (sum_original <= 0):
            loadings_novel = numpy.negative(numpy.copy(pail_components.loadings))
            sum_novel = numpy.sum(loadings_novel.flatten(order="C"))
            loading_sign_flip = True
        else:
            loadings_novel = loadings_original
            loading_sign_flip = False

        # Organize principal component factors within table.
        table_components = organize_principal_component_factor_table(
            factors=pail_components.factors, # TODO: change to factors after sign adjustment
            prefix="component_",
            index=index,
            index_name="identifier_ukb",
            report=True,
        )

    # numpy.dot
    # PC factors = dot product of (original matrix and Eigenvectors)


    pass


# TODO: TCW 17 August 2021
# TODO: this is the basic driver function for now...
# TODO: note that the main called function is "organize_principal_components_by_singular_value_decomposition()"
def organize_aggregate_metabolite_genetic_scores(
    identifier=None,
    column_index=None,
    columns_scores=None,
    table=None,
    report=None,
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
        report (bool): whether to print reports

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
    pail_decomposition = (
        organize_principal_components_by_singular_value_decomposition(
            threshold_valid_proportion_per_column=0.75,
            table=table,
            report=report,
    ))
    if False:
        pail_aggregation = organize_principal_component_aggregation(
            threshold_valid_proportion_per_column=0.75,
            table=table,
            report=report,
        )
    table_aggregation = pail_decomposition["table_components"].loc[
        :, pail_decomposition["table_components"].columns.isin([
            "identifier_ukb", "component_1"
        ])
    ]
    # Translate column names.
    translations = dict()
    translations["component_1"] = identifier
    table_aggregation.rename(
        columns=translations,
        inplace=True,
    )
    # Return.
    return table_aggregation


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
        report=report,
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
        report=report,
    )
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("Report from: read_aggregate_metabolite_genetic_scores()")
        utility.print_terminal_partition(level=2)
        print("Metabolite: " + str(metabolite))
        print(table_aggregation)
        utility.print_terminal_partition(level=3)
    # Return.
    return table_aggregation


def read_aggregate_test_metabolite_genetic_scores(
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

    # Aggregate metabolite's genetic scores.
    table_aggregation = read_aggregate_metabolite_genetic_scores(
        metabolite=metabolite,
        metabolites_files_paths=metabolites_files_paths,
        report=report,
    )
    # Copy information.
    table = table_aggregation.copy(deep=True)
    # Organize information.
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

    #for metabolite in metabolites_files_paths.keys():
    for metabolite in ["M00599", "M32315", "M02342", "M00054"]:
        # Aggregate metabolite's genetic scores.
        table_aggregation = read_aggregate_metabolite_genetic_scores(
            metabolite=metabolite,
            metabolites_files_paths=metabolites_files_paths,
            report=False,
        )
        # Copy information.
        table_metabolite = table_aggregation.copy(deep=True)
        # Organize information.
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
        # Collect information for metabolite.
        table_collection = table_collection.merge(
            table_metabolite,
            how="outer",
            left_on="identifier_ukb",
            right_on="identifier_ukb",
            suffixes=("_original", "_novel"),
        )

        pass

    # Compile information.
    #pail = dict()
    # Return information.
    return table_collection

#
##
###
#####
################################################################################

##########
# Collection


def read_select_metabolite_genetic_scores(
    metabolite=None,
    selection=None,
    metabolites_files_paths=None,
    report=None,
):
    """
    Reads a metabolite's genetic scores across the UK Biobank from file,
    and selects the scores to keep.

    This function returns a table for a single metabolite with UK Biobank
    identifiers and a single column of selection scores for the metabolite
    across these UK Biobank records.

    arguments:
        metabolite (str): identifier of a metabolite
        selection (str): name of column for selection from Polygenic Score
            thresholds
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
        report=report,
    )
    # Organize the raw table.
    table_raw.drop(
        labels=["IID",],
        axis="columns",
        inplace=True
    )
    # Select scores.
    table_selection = table_raw.loc[
        :, table_raw.columns.isin(["FID", selection])
    ]
    # Translate column names.
    translations = dict()
    translations["FID"] = "identifier_ukb"
    translations[selection] = metabolite
    table_selection.rename(
        columns=translations,
        inplace=True,
    )
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("Report from: read_select_metabolite_genetic_scores()")
        utility.print_terminal_partition(level=2)
        print("Metabolite: " + str(metabolite))
        print(table_selection)
        utility.print_terminal_partition(level=3)
    # Return.
    return table_selection


def read_select_collect_metabolites_genetic_scores(
    selection=None,
    metabolites_files_paths=None,
    report=None,
):
    """
    Reads metabolites' genetic scores across the UK Biobank from file,
    aggregates scores by Singular Value Decomposition (SVD), and collects these
    within a table.

    arguments:
        selection (str): name of column for selection from Polygenic Score
            thresholds
        metabolites_files_paths (dict<list<str>>): collection of files and paths
            for metabolites
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of metabolites' genetic scores across UK
            Biobank cohort

    """

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
        # Select metabolite's genetic scores.
        table_selection = read_select_metabolite_genetic_scores(
            metabolite=metabolite,
            selection=selection,
            metabolites_files_paths=metabolites_files_paths,
            report=False,
        )
        # Copy information.
        table_metabolite = table_selection.copy(deep=True)
        # Organize information.
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
        # Collect information for metabolite.
        table_collection = table_collection.merge(
            table_metabolite,
            how="outer",
            left_on="identifier_ukb",
            right_on="identifier_ukb",
            suffixes=("_original", "_novel"),
        )

        pass
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("Report from: read_select_collect_metabolites_genetic_scores()")
        utility.print_terminal_partition(level=2)
        print("selection: " + str(selection))
        print(table_collection)
        utility.print_terminal_partition(level=3)
    # Compile information.
    #pail = dict()
    # Return information.
    return table_collection



##########
# Write


def write_product_selection(
    selection=None,
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        selection (str): name of column for selection from Polygenic Score
            thresholds
        information (object): information to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_table = os.path.join(
        path_parent,
        str("table_metabolites_scores_prs_" + selection + ".pickle")
    )
    # Write information to file.
    information[str("table_prs_" + selection)].to_pickle(
        path_table
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

    # Simple collection table with single score for each metabolite.
    write_product_selection(
        selection="0_00001",
        information=information,
        path_parent=paths["selection"],
    )
    write_product_selection(
        selection="0_0001",
        information=information,
        path_parent=paths["selection"],
    )
    write_product_selection(
        selection="0_001",
        information=information,
        path_parent=paths["selection"],
    )
    write_product_selection(
        selection="0_01",
        information=information,
        path_parent=paths["selection"],
    )
    write_product_selection(
        selection="0_1",
        information=information,
        path_parent=paths["selection"],
    )

    # Specify directories and files.
    path_metabolites_files_paths = os.path.join(
        paths["aggregation"], "metabolites_files_paths.pickle"
    )
    # Write information to file.
    with open(path_metabolites_files_paths, "wb") as file_product:
        pickle.dump(information["metabolites_files_paths"], file_product)

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
    print("version check: 2")
    # Pause procedure.
    time.sleep(5.0)


    # Initialize directories.
    paths = initialize_directories(
        restore=True,
        path_dock=path_dock,
    )
    # Read source information from file.
    # Exclusion identifiers are "eid".
    source = read_source_initial(
        path_dock=path_dock,
        report=False,
    )

    if False:
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

        # Collect metabolites' genetic scores, and aggregate these by singular value
        # decomposition (SVD).
        # pail_metabolites_scores
        table_scores = read_aggregate_collect_metabolites_genetic_scores(
            metabolites_files_paths=source["metabolites_files_paths"],
        )
        print("printing after read_aggregate_collect_metabolites_genetic_scores")
        print(table_scores)

    # TODO: temporarily by-pass the aggregation process...
    # Collect metabolites' genetic scores at multiple PRS p-value thresholds.
    table_prs_0_00001 = read_select_collect_metabolites_genetic_scores(
        selection="X1e.05",
        metabolites_files_paths=source["metabolites_files_paths"],
        report=True,
    )
    table_prs_0_0001 = read_select_collect_metabolites_genetic_scores(
        selection="X0.0001",
        metabolites_files_paths=source["metabolites_files_paths"],
        report=True,
    )
    table_prs_0_001 = read_select_collect_metabolites_genetic_scores(
        selection="X0.001",
        metabolites_files_paths=source["metabolites_files_paths"],
        report=True,
    )
    table_prs_0_01 = read_select_collect_metabolites_genetic_scores(
        selection="X0.01",
        metabolites_files_paths=source["metabolites_files_paths"],
        report=True,
    )
    table_prs_0_1 = read_select_collect_metabolites_genetic_scores(
        selection="X0.1", # "X0.001", "X0.01", "X0.1"
        metabolites_files_paths=source["metabolites_files_paths"],
        report=True,
    )

    # Collect information.
    information = dict()
    information["metabolites_files_paths"] = source["metabolites_files_paths"]
    information["table_prs_0_00001"] = table_prs_0_00001
    information["table_prs_0_0001"] = table_prs_0_0001
    information["table_prs_0_001"] = table_prs_0_001
    information["table_prs_0_01"] = table_prs_0_01
    information["table_prs_0_1"] = table_prs_0_1
    # TODO: eventually, include a dictionary collection of a table for each
    # metabolite
    #information["pail_metabolites_scores_tables"] = pail
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )

    pass



if (__name__ == "__main__"):
    execute_procedure()
