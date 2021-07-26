import numpy as np


def check_all_ones(list_of_elements:list) -> bool:
    """
    Check if all elements in a list are 1s

    Args:
        list_of_elements (list): list of elements

    Returns:
        bool: True if all elements of the list are 1s
    """
    return all(list_of_elements == 1)


def check_mean_only(mean_only: bool) -> None:
    """
    Check mean_only option
    
    Args:
        mean_only (bool): user's choice regarding the mean_only option
    """
    if mean_only == True:
        print("Using mean only version")


def check_NAs(data):
    """check if NaNs - in theory, we construct the data without NAs

    Arguments:
        data {matrix} -- the data matrix

    Returns:
        NAs {bool} -- boolean characterising the presence of NaNs in the data matrix
    """
    # NAs = True in (np.isnan(dat))
    NAs = np.isnan(np.sum(data))  # Check if NaN exists
    if NAs:
        print("Found missing data values. Please remove all missing values before proceeding with pyComBat.")
    return NAs