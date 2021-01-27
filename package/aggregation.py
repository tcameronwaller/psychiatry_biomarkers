
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


##########
# Raw Singular Value Decomposition (SVD)


def copy_organize_table_matrix_for_singular_value_decomposition(
    threshold_valid_proportion_per_column=None,
    table=None,
    report=None,
):
    """
    Organizes a table and relevant information for Singular Value Decomposition
    (SVD).

    arguments:
        threshold_valid_proportion_per_column (float): proportion of rows that
            must have a valid value for a column in order to keep the column
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information for the singular value decomposition

    """

    # Copy information.
    table = table.copy(deep=True)
    # Drop any columns with inadequate valid values across rows.
    rows = table.shape[0]
    threshold = round(rows*threshold_valid_proportion_per_column)
    table.dropna(
        axis="columns",
        thresh=threshold,
        subset=None,
        inplace=True,
    )
    # Drop any rows with null values in any columns.
    table.dropna(
        axis="index",
        how="any",
        subset=None,
        inplace=True,
    )
    # Principal components analysis assumptions require at least centering the
    # means (mean = 0) of variables (features).
    # Standardizing the scale of variables (features) is equivalent to
    # calculation on correlation matrix instead of covariance matrix.
    # Standardize scale across variables (features) to mean zero (mean = 0) and
    # standard deviation one (standard deviation = 1).
    table_scale = utility.standardize_table_values_by_column(
        table=table,
        report=report,
    )
    # Copy information.
    index = copy.deepcopy(table_scale.index)
    # Organize matrix.
    # Matrix format has samples (cases, observations) across rows (dimension 0)
    # and variables (features) across columns (dimension 1).
    #matrix = numpy.transpose(table_scale.to_numpy())
    matrix = table_scale.to_numpy()

    # Compile information.
    pail = dict()
    pail["table_valid"] = table
    pail["table_valid_scale"] = table_scale
    pail["index"] = index
    pail["matrix"] = matrix
    pail["count_samples"] = copy.deepcopy(matrix.shape[0])
    pail["count_variables"] = copy.deepcopy(matrix.shape[1])
    # Return.
    return pail


def calculate_initial_raw_singular_value_decomposition_factors(
    threshold_valid_proportion_per_column=None,
    table=None,
    report=None,
):
    """
    Calculates the initial, raw factors of a Singular Value Decomposition (SVD).

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

    # Organize information.
    # Matrix format has samples (cases, observations) across rows (dimension 0)
    # and variables (features) across columns (dimension 1).
    pail_organization = (
        copy_organize_table_matrix_for_singular_value_decomposition(
            threshold_valid_proportion_per_column=(
                threshold_valid_proportion_per_column
            ),
            table=table,
            report=False,
    ))
    # Organize original matrix.

    # Calculate Singular Value Decomposition (SVD).
    # u: unitary matrix with left singular vectors as columns
    # s: singular values
    # vt: unitary matrix with right singular vectors as rows
    u, s, vt = scipy.linalg.svd(
        pail_organization["matrix"],
        full_matrices=False, # Full matrices do not convey more information.
        compute_uv=True,
        overwrite_a=False,
        check_finite=True,
        lapack_driver="gesdd",
    )
    # Calculate the original data matrix as a product of the SVD factors.
    s_diagonal = numpy.diag(s)
    matrix_product = numpy.dot(u, numpy.dot(s_diagonal, vt))

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_initial_raw_singular_value_decomposition_factors()"
        )
        utility.print_terminal_partition(level=2)

        # Original matrix has shape (M, N)
        print(
            "Shape of original matrix: " +
            str(pail_organization["matrix"].shape)
        )
        print("Shape of product matrix: " + str(matrix_product.shape))
        print("rows (dimension 0): samples (cases, observations)")
        print("columns (dimension 1): variables (features)")
        # M: count of samples (cases)
        # N: count of variables (features)
        # K: minimum of M or N
        # Matrix "u" has shape (M, K)
        print("Shape of matrix U (left singular vectors): " + str(u.shape))
        # Matrix "s" has shape (K, )
        # Matrix "s" is basically a one-dimensional array.
        print("Shape of matrix S (singular values): " + str(s.shape))
        # Matrix "vt" has shape (K, N)
        print(
            "Shape of matrix VT (transpose right singular vectors): " +
            str(vt.shape)
        )
        # Compare original matrix to matrix calculation from SVD factors.
        print("Compare original matrix to product of SVD factors: ")
        print(numpy.allclose(
            pail_organization["matrix"], matrix_product,
            rtol=1e-2,
            atol=1e-3,
            equal_nan=False,
        ))

        pass
    # Compile information.
    pail = dict()
    pail["table_valid_scale"] = pail_organization["table_valid_scale"]
    pail["index"] = pail_organization["index"]
    pail["matrix"] = pail_organization["matrix"]
    pail["count_samples"] = pail_organization["count_samples"]
    pail["count_variables"] = pail_organization["count_variables"]
    pail["left_singular_vectors_columns"] = u
    pail["singular_values"] = s
    pail["right_singular_vectors_rows"] = vt
    # Return.
    return pail


##########
# SVD sign adjustment

# TODO: enable this function to adjust for a specified count of singular vectors
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
    s_diagonal = numpy.diag(s)
    ut = numpy.copy(numpy.transpose(u))
    v = numpy.copy(numpy.transpose(vt))
    # Calculate basic products by matrix multiplication.
    ut_y = numpy.dot(ut, matrix)
    ut_y_square = numpy.dot(ut_y, ut_y)
    vt_y = numpy.dot(vt, matrix)
    vt_y_square = numpy.dot(vt_y, vt_y)



    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "adjust_singular_value_decomposition_factor_signs()"
        )
        utility.print_terminal_partition(level=2)
        print(
            "Shape of original matrix: " +
            str(pail_organization["matrix"].shape)
        )
        print("rows (dimension 0): samples (cases, observations)")
        print("columns (dimension 1): variables (features)")
        utility.print_terminal_partition(level=4)
        print("Shape of matrix S (singular values): " + str(s.shape))
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
        print(
            "Shape of matrix V (right singular vectors): " +
            str(v.shape)
        )
        utility.print_terminal_partition(level=4)
        print(
            "Shape of matrix UT-Y (product): " +
            str(ut_y.shape)
        )
        print(
            "Shape of matrix VT-Y (product): " +
            str(vt_y.shape)
        )
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


def calculate_principal_component_eigenvalues_from_singular_values(
    singular_values=None,
    count_samples=None,
    report=None,
):
    """
    Calculates Principal Components Analysis (PCA) Eigenvalues from Singular
    Values of Singular Value Decomposition.

    arguments:
        singular_values (object): NumPy array of Singular Values
        count_samples (float): count of samples in the original Singular Value
            Decomposition
        report (bool): whether to print reports

    raises:

    returns:
        (object): NumPy array of Eigenvalues

    """

    def divide_by_sample_count(value, count_samples):
        return (value / (count_samples - 1))
    array_divide_by_sample_count = numpy.vectorize(divide_by_sample_count)

    # Copy information.
    singular_values = numpy.copy(singular_values)
    # Calculate Eigenvalues.
    singular_values_square = numpy.square(singular_values)
    eigenvalues = array_divide_by_sample_count(
        singular_values_square, count_samples
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_principal_component_eigenvalues_from_singular_values()"
        )
        utility.print_terminal_partition(level=2)
        print("Singular values...")
        print(singular_values)
        print("Eigenvalues...")
        print(eigenvalues)
    # Return.
    return eigenvalues


def sort_eigenvectors_by_decreasing_eigenvalues(
    eigenvectors=None,
    eigenvalues=None,
    report=None,
):
    """
    Sorts Principal Components Analysis (PCA) Eigenvectors in order of
    decreasing Eigenvalues.

    arguments:
        eigenvectors (object): NumPy array matrix of Eigenvectors
        eigenvalues (object): NumPy array matrix of Eigenvalues
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    # Copy information.
    eigenvectors = numpy.copy(eigenvectors)
    eigenvalues = numpy.copy(eigenvalues)
    # Calculate sort indices for Eigenvalues in decreasing order.
    # Reverse an increasing sort order with either "numpy.flip" or "[::-1]".
    indices_sort_increasing = numpy.argsort(
        eigenvalues,
        axis=-1,
        kind="stable",
    )
    indices_sort = numpy.flip(
        indices_sort_increasing,
        axis=0,
    )
    # Apply sort order indices to Eigenvalues and Eigenvectors.
    # TODO: I'm not entirely sure that I'm sorting the correct dimension of
    # Eigenvectors...
    eigenvectors = eigenvectors[:,indices_sort]
    eigenvalues = eigenvalues[indices_sort]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Matrices before sort:")
        print("Eigenvalues...")
        print(eigenvalues)
        print("Eigenvectors...")
        print(eigenvectors)
        utility.print_terminal_partition(level=3)
        print("Matrices after sort:")
        print("Eigenvalues...")
        print(eigenvalues)
        print("Eigenvectors...")
        print(eigenvectors)
    # Compile information.
    pail = dict()
    pail["eigenvalues"] = eigenvalues
    pail["eigenvectors"] = eigenvectors
    # Return.
    return pail


def calculate_principal_component_loadings_from_eigen_values_vectors(
    eigenvectors=None,
    eigenvalues=None,
    report=None,
):
    """
    Calculates Principal Components Analysis (PCA) loadings from Eigenvectors
    and Eigenvalues.

    Statsmodels erroneously returns "loadings" that have identical values and
    dimensions as the Eigenvectors; however, Eigenvectors are loadings are not
    equivalent.

    loadings = eigenvectors [dot] square_root(eigenvalues)
    Loadings include aspects of both direction (eigenvectors) and scale
    (eigenvalues).

    Reference:
    "https://stats.stackexchange.com/questions/143905/
    loadings-vs-eigenvectors-in-pca-when-to-use-one-or-another"

    arguments:
        eigenvectors (object): NumPy array matrix of Eigenvectors
        eigenvalues (object): NumPy array matrix of Eigenvalues
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    # Copy information.
    eigenvectors = numpy.copy(eigenvectors)
    eigenvalues = numpy.copy(eigenvalues)
    # Calculate square roots of Eigenvalues.
    # Organize a diagonal matrix of square roots of Eigenvalues.
    eigenvalues_square_root = numpy.sqrt(eigenvalues)
    eigenvalues_root_diagonal = numpy.diag(eigenvalues_square_root)
    # Calculate loadings.
    loadings = numpy.dot(eigenvectors, eigenvalues_root_diagonal)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_principal_component_loadings_from_eigen_values_vectors()"
        )
        utility.print_terminal_partition(level=2)
        print("Loadings = Eigenvectors [dot] square_root(diagonal Eigenvalues)")
        print(loadings)
    # Return.
    return loadings


def calculate_principal_component_loadings_from_direct_factors(
    s=None,
    vh=None,
    count_samples=None,
    report=None,
):
    """
    Calculates Principal Components Analysis (PCA) loadings from direct factors
    of Singular Value Decomposition.

    Reference:
    "https://stats.stackexchange.com/questions/134282/
    relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca"

    arguments:
        s (object): NumPy array of Singular Values
        vh (object): NumPy unitary matrix with right singular vectors as rows
        count_samples (float): count of samples in the original Singular Value
            Decomposition
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    def divide_by_sample_count(value, count_samples):
        return (value / math.sqrt(count_samples - 1))
    array_divide_by_sample_count = numpy.vectorize(divide_by_sample_count)


    # Copy information.
    s = numpy.copy(s)
    vh = numpy.copy(vh)
    # Calculate loadings.
    v = numpy.transpose(vh)
    s_diagonal = numpy.diag(s)
    product = numpy.dot(v, s_diagonal)
    loadings = array_divide_by_sample_count(
        product, count_samples
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_principal_component_loadings_from_direct_factors()"
        )
        utility.print_terminal_partition(level=2)
        print("Loadings = (V [dot] S) / square_root(samples - 1)")
        print(loadings)
    # Return.
    return loadings


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



def organize_principal_components_by_singular_value_decomposition(
    threshold_valid_proportion_per_column=None,
    table=None,
    report=None,
):
    """
    Organizes a Singular Value Decomposition (SVD).

    Reference:
    "https://stats.stackexchange.com/questions/134282/
    relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca"

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

    # Calculate initial, raw Singular Value Decomposition factors.
    pail_svd = calculate_initial_raw_singular_value_decomposition_factors(
        threshold_valid_proportion_per_column=(
            threshold_valid_proportion_per_column
        ),
        table=table,
        report=True,
    )
    #print(pail_svd["table_valid_scale"])
    #print(pail_svd["index"])
    #print(pail_svd["count_samples"])
    #print(pail_svd["count_variables"])
    # Adjust the signs of singular vectors.
    pail_svd_sign = adjust_singular_value_decomposition_factor_signs(
        matrix=pail_svd["matrix"],
        singular_values=pail_svd["singular_values"],
        left_singular_vectors_columns=pail_svd["left_singular_vectors_columns"],
        right_singular_vectors_rows=pail_svd["right_singular_vectors_rows"],
        report=True,
    )

    # TODO: call the function to derive PCA info from the adjusted SVD factors

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "organize_principal_components_by_singular_value_decomposition()"
        )
        utility.print_terminal_partition(level=2)

        pass

    if False:
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
    pass



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



# TODO: PCA loadings adjustment...
# TODO: organize this better...
# TODO: give useful metrics on the PCA
# TODO: return loadings
# TODO: return variances (Eigenvalues?)
def organize_principal_component_aggregation(
    threshold_valid_proportion_per_column=None,
    table=None,
    report=None,
):
    """
    Organizes aggregation by modification of Principal Components Analysis
    (PCA).

    Factorize covariance or correlation into direction (eigenvectors) and scale
    (eigenvalues).
    The eigenvalues impart scale or weight to each eigenvector.
    Each eigenvector has its own eigenvalue, and their sort orders mush match.

    loadings = eigenvectors [dot] square_root(eigenvalues)
    Loadings include aspects of both direction (eigenvectors) and scale
    (eigenvalues).

    Singular Value Decomposition assigns Eigenvector direction (sign, positive
    or negative) at random. To ensure consistency, flip all signs such that the
    sum of loadings (direction and scale) is positive.

    arguments:
        threshold_valid_proportion_per_column (float): proportion of rows that
            must have a valid value for a column in order to keep the column
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with explicit
            index
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the principal components
            aggregation

    """

    # Copy information.
    table = table.copy(deep=True)
    # Drop any columns with inadequate valid values across rows.
    rows = table.shape[0]
    threshold = round(rows*threshold_valid_proportion_per_column)
    table.dropna(
        axis="columns",
        thresh=threshold,
        subset=None,
        inplace=True,
    )
    # Drop any rows with null values in any columns.
    table.dropna(
        axis="index",
        how="any",
        subset=None,
        inplace=True,
    )
    # Calculate principal components, forcing loadings to have positive sum.
    pail_components = (
        organize_principal_component_analysis_force_loadings_positive_sum(
            table=table,
            report=True,
    ))
    # Eigenvalues.

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Report from: organize_principal_component_aggregation()")
        utility.print_terminal_partition(level=2)
        print(
            "Threshold on rows with valid values per column: " + str(threshold)
        )
        print("Shape of original matrix: " + str(matrix.shape))
        print(
            "Shape of Principal Components: " +
            str(pail_components.factors.shape)
        )
        print("Shape of loadings: " + str(pail_components.loadings.shape))
        print("sum of loadings original...")
        print(sum_original)
        print("sum of loadings novel...")
        print(sum_novel)
        print("Shape of Eigenvalues: " + str(pail_components.eigenvals.shape))
        print("Are the Eigenvalues in decreasing order?")
        print(pail_components.eigenvals)
        print("Shape of Eigenvectors: " + str(pail_components.eigenvecs.shape))
        utility.print_terminal_partition(level=3)
        print("table_components")
        utility.print_terminal_partition(level=4)
        print(table_components)
        pass

    # Compile information.
    pail = dict()
    pail["loadings"] = []
    pail["variances"] = []
    pail["table_components"] = table_components
    # Return.
    return pail


# TODO: this is the basic driver function for now...
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
    source = read_source_initial(
        path_dock=path_dock,
        report=False,
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
        print("printing after read_aggregate_collect_metabolites_genetic_scores")
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
