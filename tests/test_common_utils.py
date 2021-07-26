import pytest

import numpy as np
import pandas as pd

from combat.utils.common_utils import check_all_ones, check_mean_only, check_NAs


# tests for check_all_ones function
def test_check_all_ones():
    assert check_all_ones(np.array([1, 1, 1, 1, 1])) == True

    assert check_all_ones(np.array([1, 1, 1, 1, 0])) == False
    assert check_all_ones(np.array([0, 0, 0, 0, 0])) == False

    assert check_all_ones(np.array([1.5, 0.5, 1, 1, 1])) == False # This test to show the limit of the method we use


# test for check_mean_only
def test_check_mean_only():
    check_mean_only(True)
    check_mean_only(False)
    print("Only one line of text should have been printed above.")


# test for check_NAs
def test_check_NAs():
    assert check_NAs([0,1,2]) == False
    assert check_NAs([0,np.nan,2]) == True