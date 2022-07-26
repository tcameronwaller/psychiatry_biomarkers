"""
...
"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Import modules from specific path without having to install a general package
# I would have to figure out how to pass a path variable...
# https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path

# Standard.
import argparse
import textwrap

# Relevant.

# Custom.

import scratch
import uk_biobank.interface
import stragglers.interface

#dir()
#importlib.reload()

###############################################################################
# Functionality


# Parser and management of subparsers.


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
    subparsers = parser.add_subparsers(title="routines")
    parser_main = define_subparser_main(subparsers=subparsers)
    parser_uk_biobank = uk_biobank.interface.define_subparser_main(
        subparsers=subparsers
    )
    parser_stragglers = stragglers.interface.define_subparser_main(
        subparsers=subparsers
    )
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

        Manage routines and procedures for Psychiatric Metabolism project.

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


# Package's main subparser.


def define_subparser_main(subparsers=None):
    """
    Defines subparser and parameters.

    arguments:
        subparsers (object): reference to subparsers' container

    raises:

    returns:
        (object): reference to subparser

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
    if arguments.scratch:
        # Report status.
        print("... executing 'scratch' procedure ...")
        # Execute procedure.
        scratch.execute_procedure(
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

    # Parse arguments from terminal.
    arguments = define_interface_parsers()
    # Call the appropriate functions.
    arguments.func(arguments)


if (__name__ == "__main__"):
    execute_procedure()
