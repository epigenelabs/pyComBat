# -----------------------------------------------------------------------------
# Copyright (C) 2019-2020 A. Behdenna, A. Nordor, J. Haziza and A. Gema

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# For more information, please contact Abdelkader Behdenna <abdelkader@epigenelabs.com>/<kaderbehdenna@gmail.com>

# file 	pycombat.py
# author A. Behdenna, J. Haziza, A. Gema, A. Nordor
# date Sept 2020
# -----------------------------------------------------------------------------


# this file is only used for unit testing
# We import the function that will be tested one by one, incrementally
# Generates a report about the functions tested

import os
import sys

sys.path.append(os.getcwd())

import numpy as np
import pandas as pd

##########
# import function used for unit testing
from combat import pycombat

##########
print("\n#### Unit Testing for pyComBat ####\n")

##########
# Define constants for unit testing

batch = np.asarray([1, 1, 1, 2, 2, 3, 3, 3, 3])
matrix = np.transpose(
    [
        np.random.normal(size=1000, loc=3, scale=1),
        np.random.normal(size=1000, loc=3, scale=1),
        np.random.normal(size=1000, loc=3, scale=1),
        np.random.normal(size=1000, loc=2, scale=0.6),
        np.random.normal(size=1000, loc=2, scale=0.6),
        np.random.normal(size=1000, loc=4, scale=1),
        np.random.normal(size=1000, loc=4, scale=1),
        np.random.normal(size=1000, loc=4, scale=1),
        np.random.normal(size=1000, loc=4, scale=1),
    ]
)

matrix = pd.DataFrame(
    data=matrix,
    columns=["sample_" + str(i + 1) for i in range(9)],
    index=["gene_" + str(i + 1) for i in range(1000)],
)

print("Matrix and batch generated.")

matrix_adjusted = pycombat(matrix, batch)

print("Adjusted matrix generated.")

##########
# local tests before unit testing


def test_means():
    print(f"mean matrix: {np.mean(matrix)}")
    print(f"mean matrix (adjusted): {np.mean(matrix_adjusted)}")
    print("**********")
    print(f"var matrix: {np.var(matrix)}")
    print(f"var matrix (adjusted): {np.var(matrix_adjusted)}")


# general tests on pyComBat
def test_pycombat():
    assert np.shape(matrix) == np.shape(matrix_adjusted)
    assert (
        abs(np.mean(matrix.values) - np.mean(matrix_adjusted.values))
        <= np.mean(matrix.values) * 0.05
    )
    assert np.var(matrix_adjusted.values) <= np.var(matrix.values)
