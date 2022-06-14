"""
...
"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard.
import argparse
import textwrap

# Relevant.

# Custom.

#import genetic_correlation
#import aggregation


import assembly
import importation
import organization
import stratification
import description
import regression
import collection
import scratch

# Call subpackage module execution functions directly.
import bipolar_biobank.assembly
import bipolar_biobank.organization
import bipolar_biobank.regression

#import plot
#import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def define_interface_parsers():
    """
    Defines and parses arguments from terminal's interface.

    arguments:

    raises:

    returns:
        (object): arguments from terminal

    """

    # Define description.
    description = define_description_general()
    # Define epilog.
    epilog = define_epilog_general()
    # Define arguments.
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(title="procedures")
    parser_main = define_subparser_main(subparsers=subparsers)
    parser_bipolar_biobank = define_subparser_bipolar_biobank(
        subparsers=subparsers
    )
    # TODO: add other subparsers here...
    # Parse arguments.
    return parser.parse_args()


def define_description_general():
    """
    Defines description for terminal interface.

    arguments:

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------

        Access data from UK Biobank and do other stuff.

        --------------------------------------------------
    """)
    return description


def define_epilog_general():
    """
    Defines epilog for terminal interface.

    arguments:

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\

        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """)
    return epilog


def define_subparser_main(subparsers=None):
    """
    Defines subparser for procedures that adapt a model of human metabolism.

    arguments:
        subparsers (object): reference to subparsers' container

    raises:

    returns:
        (object): reference to parser

    """

    # Define description.
    description = define_description_main()
    # Define epilog.
    epilog = define_epilog_main()
    # Define parser.
    parser = subparsers.add_parser(
        name="main",
        description=description,
        epilog=epilog,
        help="Help for main routine.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Define arguments.
    parser.add_argument(
        "-path_dock", "--path_dock", dest="path_dock", type=str, required=True,
        help=(
            "Path to dock directory for source and product " +
            "directories and files."
        )
    )
    parser.add_argument(
        "-genetic_correlation", "--genetic_correlation",
        dest="genetic_correlation",
        action="store_true",
        help=(
            "Genetic correlations for metabolites from multiple GWAS."
        )
    )
    parser.add_argument(
        "-aggregation", "--aggregation", dest="aggregation",
        action="store_true",
        help=(
            "Aggregation of genetic scores for metabolites across UK Biobank."
        )
    )
    parser.add_argument(
        "-assembly", "--assembly", dest="assembly",
        action="store_true",
        help=(
            "Assemble phenotype information from UK Biobank."
        )
    )
    parser.add_argument(
        "-importation", "--importation", dest="importation",
        action="store_true",
        help=(
            "Assemble phenotype information from UK Biobank."
        )
    )
    parser.add_argument(
        "-organization", "--organization", dest="organization",
        action="store_true",
        help=(
            "Organize phenotype information from UK Biobank."
        )
    )
    parser.add_argument(
        "-stratification", "--stratification",
        dest="stratification",
        action="store_true",
        help=(
            "Stratification of cohorts and formatting tables of phenotypes " +
            "and covariates for genetic analyses (especially GWAS in PLINK2)."
        )
    )
    parser.add_argument(
        "-description", "--description",
        dest="description",
        action="store_true",
        help=(
            "Description of cohorts and phenotypes with summary statistics " +
            "and plots."
        )
    )
    parser.add_argument(
        "-regression", "--regression",
        dest="regression",
        action="store_true",
        help=(
            "Regression analyses of phenotypes within cohorts."
        )
    )
    parser.add_argument(
        "-collection", "--collection",
        dest="collection",
        action="store_true",
        help=(
            "Collection and summary of reports from genotypic analyses."
        )
    )
    parser.add_argument(
        "-scratch", "--scratch",
        dest="scratch",
        action="store_true",
        help=(
            "Scratch analyses."
        )
    )
    # Define behavior.
    parser.set_defaults(func=evaluate_parameters_main)
    # Return parser.
    return parser


def define_subparser_bipolar_biobank(subparsers=None):
    """
    Defines subparser for procedures that adapt a model of human metabolism.

    arguments:
        subparsers (object): reference to subparsers' container

    raises:

    returns:
        (object): reference to parser

    """

    # Define description.
    description = define_description_bipolar_biobank()
    # Define epilog.
    epilog = define_epilog_bipolar_biobank()
    # Define parser.
    parser = subparsers.add_parser(
        name="bipolar_biobank",
        description=description,
        epilog=epilog,
        help="Help for bipolar_biobank routine.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Define arguments.
    parser.add_argument(
        "-path_dock", "--path_dock", dest="path_dock", type=str, required=True,
        help=(
            "Path to dock directory for source and product " +
            "directories and files."
        )
    )
    parser.add_argument(
        "-assembly", "--assembly", dest="assembly",
        action="store_true",
        help=(
            "Assemble phenotype information from Bipolar Biobank."
        )
    )
    parser.add_argument(
        "-organization", "--organization", dest="organization",
        action="store_true",
        help=(
            "Organize phenotype information from Bipolar Biobank."
        )
    )
    parser.add_argument(
        "-regression", "--regression", dest="regression",
        action="store_true",
        help=(
            "Regression analyses on phenotype information from Bipolar Biobank."
        )
    )

    # Define behavior.
    parser.set_defaults(func=evaluate_parameters_bipolar_biobank)
    # Return parser.
    return parser


def define_description_main():
    """
    Defines description for terminal interface.

    arguments:

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------

        Package's main procedure

        Do stuff.

        --------------------------------------------------
    """)
    return description


def define_epilog_main():
    """
    Defines epilog for terminal interface.

    arguments:

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\

        --------------------------------------------------
        main routine

        --------------------------------------------------
        additional notes...


        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """)
    return epilog


def define_description_bipolar_biobank():
    """
    Defines description for terminal interface.

    arguments:

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------

        Package's main procedure

        Do stuff.

        --------------------------------------------------
    """)
    return description


def define_epilog_bipolar_biobank():
    """
    Defines epilog for terminal interface.

    arguments:

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\

        --------------------------------------------------
        main routine

        --------------------------------------------------
        additional notes...


        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """)
    return epilog


def evaluate_parameters_main(arguments):
    """
    Evaluates parameters for model procedure.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    print("--------------------------------------------------")
    print("... call to main routine ...")
    # Execute procedure.
    if arguments.genetic_correlation:
        # Report status.
        print("... executing genetic_correlation procedure ...")
        # Execute procedure.
        genetic_correlation.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.aggregation:
        # Report status.
        print("... executing aggregation procedure ...")
        # Execute procedure.
        aggregation.execute_procedure(
            path_dock=arguments.path_dock
        )

    if arguments.assembly:
        # Report status.
        print("... executing assembly procedure ...")
        # Execute procedure.
        assembly.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.importation:
        # Report status.
        print("... executing importation procedure ...")
        # Execute procedure.
        importation.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.organization:
        # Report status.
        print("... executing organization procedure ...")
        # Execute procedure.
        organization.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.stratification:
        # Report status.
        print("... executing stratification procedure ...")
        # Execute procedure.
        stratification.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.description:
        # Report status.
        print("... executing description procedure ...")
        # Execute procedure.
        description.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.regression:
        # Report status.
        print("... executing regression procedure ...")
        # Execute procedure.
        regression.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.collection:
        # Report status.
        print("... executing collection procedure ...")
        # Execute procedure.
        collection.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.scratch:
        # Report status.
        print("... executing scratch procedure ...")
        # Execute procedure.
        scratch.execute_procedure(
            path_dock=arguments.path_dock
        )
    pass


def evaluate_parameters_bipolar_biobank(arguments):
    """
    Evaluates parameters for model procedure.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    print("--------------------------------------------------")
    print("... call to bipolar_biobank routine ...")
    # Execute procedure.
    if arguments.assembly:
        # Report status.
        print("... executing assembly procedure ...")
        # Execute procedure.
        bipolar_biobank.assembly.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.organization:
        # Report status.
        print("... executing organization procedure ...")
        # Execute procedure.
        bipolar_biobank.organization.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.regression:
        # Report status.
        print("... executing regression procedure ...")
        # Execute procedure.
        bipolar_biobank.regression.execute_procedure(
            path_dock=arguments.path_dock
        )
    pass



###############################################################################
# Procedure


def execute_procedure():
    """
    Function to execute module's main behavior.

    arguments:

    returns:

    raises:

    """

    # TODO: I want 2 separate procedures: 1. definition, 2. analysis

    # Parse arguments from terminal.
    arguments = define_interface_parsers()
    # Call the appropriate function.
    arguments.func(arguments)


if (__name__ == "__main__"):
    execute_procedure()
