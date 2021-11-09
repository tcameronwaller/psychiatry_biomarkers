
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
import promiscuity.regression as regression
import uk_biobank.organization as ukb_organ
import uk_biobank.stratification as ukb_strat
import uk_biobank.description as ukb_descr
import uk_biobank.regression as ukb_regre

###############################################################################
# Functionality


################################################################################
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

    # Execute assembly procedure from uk_biobank package.
    ukb_regre.execute_procedure(path_dock=path_dock)
    utility.print_terminal_partition(level=1)
    print("From package 'uk_biobank', procedure 'regression' is complete.")

    pass


if (__name__ == "__main__"):
    execute_procedure()
