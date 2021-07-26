from typing_extensions import (
    Literal,
    Tuple
)
from typing import List, Optional

import numpy as np
from functools import partial
import mpmath as mp
import pandas as pd

from .common_utils import check_all_ones


def model_matrix(info: list, intercept: bool = True, drop_first: bool = True) -> np.array:
    """
    Creates the model_matrix from batch list

    Args:
        info (list): list info with batch or covariates data
        intercept (bool, optional): boolean for intercept in model matrix. Defaults to True.
        drop_first (bool, optional): boolean for drop the first row. Defaults to True.

    Returns:
        np.array: Model matrix generate from batch list
    """
    if not isinstance(info[0],list) :
        info = [info]
    else:
        info = info
    info_dict = {}
    for i in range(len(info)):
        info_dict[f"col{str(i)}"] = list(map(str,info[i]))
    df = pd.get_dummies(pd.DataFrame(info_dict), drop_first=drop_first, dtype=float)
    if intercept:
        df["intercept"] = 1.0
    return df.to_numpy()


def compute_prior(prior: Literal["a", "b"], gamma_hat: np.array, mean_only: bool) -> float:
    """
    aprior and bprior are useful to compute "hyper-prior values"
    -> prior parameters used to estimate the prior gamma distribution for multiplicative batch effect
    aprior - calculates empirical hyper-prior values

    Arguments:
        prior {char} -- 'a' or 'b' depending of the prior to be calculated
        gamma_hat {matrix} -- matrix of additive batch effect
        mean_only {bool} -- True iff mean_only selected

    Returns:
        float -- [the prior calculated (aprior or bprior)

    Args:
        prior (str): 'a' or 'b' depending of the prior to be calculated
        gamma_hat (np.array): [description]
        mean_only (bool): [description]

    Returns:
        float: [description]
    """
    if mean_only:
        return 1

    m = np.mean(gamma_hat)
    s2 = np.var(gamma_hat)

    if prior == 'a':
        calculated_prior = (2*s2+m*m)/s2
    elif prior == 'b':
        calculated_prior = (m*s2+m*m*m)/s2

    return calculated_prior


def postmean(g_bar: np.array, d_star: np.array, t2_n: np.array, t2_n_g_hat: np.array) -> np.array:
    """
    Estimates additive batch effect

    Args:
        g_bar (np.array): matrix of additive batch effect
        d_star (np.array): matrix of multiplicative batch effect
        t2_n (np.array): [description]
        t2_n_g_hat (np.array): [description]

    Returns:
        np.array: matrix of estimated additive batch effect
    """

    return np.divide(t2_n_g_hat + d_star * g_bar, np.asarray(t2_n + d_star))


def postvar(sum2: np.array, n: float, a: float, b: float) -> np.array:
    """
    Estimates multiplicative batch effect

    Args:
        sum2 (np.array): Vector
        n (float): [description]
        a (float): aprior
        b (float): bprior

    Returns:
        np.array: matrix of estimated multiplicative batch effect
    """
    return(np.divide((np.multiply(0.5, sum2)+b), (np.multiply(0.5, n)+a-1)))


def it_sol(sdat: np.array, g_hat: np.array, d_hat: np.array, g_bar: np.array, t2: np.array, a: float, b: float, conv: float = 0.0001, exit_iteration: float = 10e5) -> np.array:
    """
    Iterative solution for Empirical Bayesian method

    Args:
        sdat (np.array): data matrix
        g_hat (np.array): matrix of the average additive batch effect
        d_hat (np.array): matrix of the average multiplicative batch effect
        g_bar (np.array): matrix of the additive batch effect
        t2 (np.array): [description]
        a (float): aprior
        b (float): bprior
        conv (float, optional): convergence criterion. Defaults to 0.0001.
        exit_iteration (float, optional): maximum number of iterations before exit. Defaults to 10e5.

    Returns:
        np.array: An array of array with the estimated additive and multiplicative batch effect
    """

    n = [len(i) for i in np.asarray(sdat)]
    t2_n = np.multiply(t2, n)
    t2_n_g_hat = np.multiply(t2_n, g_hat)
    g_old = np.ndarray.copy(g_hat)
    d_old = np.ndarray.copy(d_hat)
    change = 1
    count = 0  # number of steps needed (for diagnostic only)
    
    # convergence criteria, if new-old < conv, then stop
    while (change > conv) and (count < exit_iteration):
        g_new = postmean(g_bar, d_old, t2_n, t2_n_g_hat)  # updated additive batch effect
        sum2 = np.sum(np.asarray(np.square(
            sdat-np.outer(g_new[0][0], np.ones(np.ma.size(sdat, axis=1))))), axis=1)
        d_new = postvar(sum2, n, a, b)  # updated multiplicative batch effect
        change = max(np.amax(np.absolute(g_new-np.asarray(g_old))/np.asarray(g_old)), np.amax(
            np.absolute(d_new-d_old)/d_old))  # maximum difference between new and old estimate
        g_old = np.ndarray.copy(g_new)  # save value for g
        d_old = np.ndarray.copy(d_new)  # save value for d
        count += 1
    adjust = np.asarray([g_new, d_new])
    
    return adjust


def int_eprior(sdat: np.array, g_hat: np.array, d_hat: np.array, precision: float) -> np.array:
    """
    int_eprior - Monte Carlo integration function to find nonparametric adjustments
    Johnson et al (Biostatistics 2007, supp.mat.) show that we can estimate the multiplicative and additive batch effects with an integral
    This integral is numerically computed through Monte Carlo inegration (iterative method)

    Args:
        sdat (np.array): data matrix
        g_hat (np.array): matrix of the average additive batch effect
        d_hat (np.array): matrix of the average multiplicative batch effect
        precision (float): level of precision for precision computing

    Returns:
        np.array: An array of array with the estimated additive and multiplicative batch effect
    """
    g_star = []
    d_star = []
    # use this variable to only print error message once if approximation used
    test_approximation = 0
    for i in range(len(sdat)):
        # additive batch effect
        g = np.asarray(np.delete(np.transpose(g_hat), i))

        # multiplicative batch effect
        d = np.asarray(np.delete(np.transpose(d_hat), i))
        x = np.asarray(np.transpose(sdat[i]))
        n = len(x)
        j = [1]*n
        dat = np.repeat(x, len(np.transpose(g)), axis=1)
        resid2 = np.square(dat-g)
        sum2 = np.dot(np.transpose(resid2), j)
        
        # /begin{handling high precision computing}
        temp_2d = 2*d
        if (precision == None):
            LH = np.power(1 / (np.pi * temp_2d), n / 2) * np.exp(np.negative(sum2) / (temp_2d))
        else:
            # only if precision parameter informed
            # increase the precision of the computing (if negative exponential too close to 0)
            mp.dps = precision
            buf_exp = list(map(mp.exp, np.negative(sum2)/(temp_2d)))
            buf_pow = list(map(partial(mp.power, y=n/2), 1/(np.pi*temp_2d)))
            LH = buf_pow * buf_exp  # likelihood
        # /end{handling high precision computing}
        LH = np.nan_to_num(LH)  # corrects NaNs in likelihood
        if np.sum(LH) == 0 and test_approximation == 0:  # correction for LH full of 0.0
            test_approximation = 1  # this message won't appear again
            print("###\nValues too small, approximation applied to avoid division by 0.\nPrecision mode can correct this problem, but increases computation time.\n###")
            LH[LH == 0] = np.exp(-745)
            g_star.append(np.sum(g * LH) / np.sum(LH))
            d_star.append(np.sum(d * LH) / np.sum(LH))
        else:
            g_star.append(np.sum(g * LH) / np.sum(LH))
            d_star.append(np.sum(d * LH) / np.sum(LH))
    adjust = np.asarray([np.asarray(g_star), np.asarray(d_star)])
    
    return adjust


def param_fun(i: int, standardised_data: np.array, batches: List[list], mean_only: bool, gamma_hat: np.array, gamma_bar: np.array, delta_hat: np.array, t2: np.array, a_prior: float, b_prior: float) -> Tuple[np.array, np.array]:
    """
    Parametric estimation of batch effects

    Args:
        i (int): column index
        standardised_data (np.array): data matrix
        batches (List[list]): list of list of batches' elements
        mean_only (bool): True if mean_only selected
        gamma_hat (np.array): average additive batch effect
        gamma_bar (np.array): estimated additive batch effect
        delta_hat (np.array): average multiplicative batch effect
        t2 (np.array): [description]
        a_prior (float): aprior
        b_prior (float): bprior

    Returns:
        Tuple[np.array, np.array]: estimated adjusted additive (gamma_star) and multiplicative (delta_star) batch effect
    """
    if mean_only:  # if mean_only, no need for complex method: batch effect is immediately calculated
        t2_n = np.multiply(t2[i], 1)
        t2_n_g_hat = np.multiply(t2_n, gamma_hat[i])
        gamma_star = postmean(gamma_bar[i], 1, t2_n, t2_n_g_hat)  # additive batch effect
        delta_star = [1]*len(standardised_data)  # multiplicative batch effect
    else:  # if not(mean_only) then use it_solve
        temp = it_sol(np.transpose(np.transpose(standardised_data)[
                      batches[i]]), gamma_hat[i], delta_hat[i], gamma_bar[i], t2[i], a_prior[i], b_prior[i])
        gamma_star = temp[0]  # additive batch effect
        delta_star = temp[1]  # multiplicative batch effect
    return (gamma_star, delta_star)


# FIXME: Params order is not the same
def nonparam_fun(i: int, mean_only: bool, delta_hat: np.array, standardised_data: np.array, batches: List[list], gamma_hat: np.array, precision: float) -> Tuple[np.array, np.array]:
    """
    Non-parametric estimation

    Args:
        i (int): column index
        mean_only (bool): True if mean_only selected
        delta_hat (np.array): estimated multiplicative batch effect
        standardised_data (np.array): data matrix
        batches (List[list]): list of list of batches' elements
        gamma_hat (np.array): estimated additive batch effect
        precision (float): level of precision for precision computing

    Returns:
        Tuple[np.array, np.array]: estimated adjusted additive and multiplicative batch effect
    """
    if mean_only:  # if mean only, change delta_hat to vector of 1s
        delta_hat[i] = [1]*len(delta_hat[i])
    # use int_eprior for non-parametric estimation
    temp = int_eprior(np.transpose(np.transpose(standardised_data)[
                      batches[i]]), gamma_hat[i], delta_hat[i], precision)
    return (temp[0], temp[1])


def check_ref_batch(ref_batch: int, batch: List[int], batchmod: np.array) -> Tuple[Optional[List[int]], np.array]:
    """
    check ref_batch option and treat it if needed

    Args:
        ref_batch (int): the reference batch
        batch (List[int]): list of batch id
        batchmod (np.array): model matrix related to batches

    Raises:
        ValueError: Reference level ref_batch must be one of the levels of batch

    Returns:
        Tuple[Optional[List[int]], np.array]: 
            - the corresponding positions of the reference batch in the batch list
            - updated model matrix related to batches, with reference
    """
    if ref_batch is not None:
        if ref_batch not in batch:
            raise ValueError("Reference level ref_batch must be one of the levels of batch.")

        print(f"Using batch {str(ref_batch)} as a reference batch.")
        # ref keeps in memory the columns concerned by the reference batch
        ref = np.where(np.unique(batch) == ref_batch)[0][0]
        # updates batchmod with reference
        batchmod[:, ref] = 1
    else:
        ref = None  # default settings

    return (ref, batchmod)


def treat_covariates(batchmod: np.array, mod: np.array, ref: int, n_batch: int) -> np.array:
    """
    treat covariates

    Args:
        batchmod (np.array): model matrix for batch
        mod (np.array): model matrix for other covariates
        ref (int): index of reference batch
        n_batch (int): number of batches

    Raises:
        ValueError: The covariate is confunded with batch. Remove the covariate and rerun pyComBat.
        ValueError: The covariates are confounded! Please remove one or more of the covariates so the design is not confounded.
        ValueError: At least one covariate is confounded with batch. Please remove confounded covariates and rerun pyComBat.

    Returns:
        np.array: model matrix for all covariates, including batch
    """
    # design matrix for sample conditions
    if mod == []:
        design = batchmod
    else:
        mod_matrix = model_matrix(mod, intercept=True)
        design = np.concatenate((batchmod, mod_matrix), axis=1)
    check = list(map(check_all_ones, np.transpose(design)))
    if ref is not None:
        check[ref] = False  # the reference in not considered as a covariate
    design = design[:, ~np.array(check)]
    design = np.transpose(design)

    print(f"Adjusting for {str(len(design)-len(np.transpose(batchmod)))} covariate(s) or covariate level(s).")

    # if matrix cannot be invertible, different cases
    if np.linalg.matrix_rank(design) < len(design):
        if len(design) == n_batch + 1:  # case 1: covariate confunded with a batch
            raise ValueError("The covariate is confunded with batch. Remove the covariate and rerun pyComBat.")
        if len(design) > n_batch + 1:  # case 2: multiple covariates confunded with a batch
            if np.linalg.matrix_rank(np.transpose(design)[:n_batch]) < len(design):
                raise ValueError("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded.")
            else:  # case 3: at least a covariate confunded with a batch
                raise ValueError("At least one covariate is confounded with batch. Please remove confounded covariates and rerun pyComBat")
    return design


def treat_batches(batch: list) -> Tuple[int, List[int], List[int], int]:
    """
    treat batches

    Args:
        batch (list): batch list

    Returns:
        Tuple[int, List[int], List[int], int]:
            - number of batches
            - list of unique batches
            - list of batches lengths 
            - total size of dataset
    """
    n_batch = len(np.unique(batch))  # number of batches
    print(f"Found {n_batch} batches.")

    batches = []  # list of lists, contains the list of position for each batch
    for i in range(n_batch):
        batches.append(np.where(batch == np.unique(batch)[i])[0])
    batches = np.asarray(batches)
    n_batches = list(map(len, batches))
    if 1 in n_batches:
        #mean_only = True  # no variance if only one sample in a batch - mean_only has to be used
        print("\nOne batch has only one sample, try setting mean_only=True.\n")
    n_array = sum(n_batches)
    return (n_batch, batches, n_batches, n_array)


def calculate_mean_var(design: np.array, batches: List[int], ref: int, dat: np.array, NAs: bool, ref_batch: int, n_batches: List[int], n_batch: int, n_array: int) -> Tuple[np.array, np.array, np.array]:
    """
    calculates the Normalisation factors

    Args:
        design (np.array): model matrix for all covariates
        batches (List[int]): list of unique batches
        ref (int): index of reference batch
        dat (np.array): data matrix
        NAs (bool): presence of NaNs in the data matrix
        ref_batch (int): reference batch
        n_batches (List[int]): list of batches lengths
        n_batch (int): number of batches
        n_array (int): total size of dataset

    Returns:
        Tuple[np.array, np.array, np.array]:
            - regression coefficients corresponding to the design matrix
            - Mean for each gene and each batch
            - Variance for each gene and each batch
    """
    print("Standardizing Data across genes.")
    if not(NAs):  # NAs not supported
        # B_hat is the vector of regression coefficients corresponding to the design matrix
        B_hat = np.linalg.solve(np.dot(design, np.transpose(
            design)), np.dot(design, np.transpose(dat)))

    # Calculates the general mean
    if ref_batch is not None:
        grand_mean = np.transpose(B_hat[ref])
    else:
        grand_mean = np.dot(np.transpose(
            [i / n_array for i in n_batches]), B_hat[0:n_batch])
    
    # Calculates the general variance
    if not NAs:  # NAs not supported
        if ref_batch is not None:  # depending on ref batch
            ref_dat = np.transpose(np.transpose(dat)[batches[ref]])
            var_pooled = np.dot(np.square(ref_dat - np.transpose(np.dot(np.transpose(
                design)[batches[ref]], B_hat))), [1/n_batches[ref]]*n_batches[ref])
        else:
            var_pooled = np.dot(np.square(
                dat - np.transpose(np.dot(np.transpose(design), B_hat))), [1/n_array]*n_array)

    return (B_hat, grand_mean, var_pooled)


def calculate_stand_mean(grand_mean: np.array, n_array: int, design: np.array, n_batch: int, B_hat: np.array) -> np.array:
    """
    transform the format of the mean for substraction

    Args:
        grand_mean (np.array): Mean for each gene and each batch
        n_array (int): total size of dataset
        design (np.array): design matrix for all covariates including batch
        n_batch (int): number of batches
        B_hat (np.array): regression coefficients corresponding to the design matrix

    Returns:
        np.array: standardised mean
    """
    stand_mean = np.dot(np.transpose(np.mat(grand_mean)), np.mat([1]*n_array))
    # corrects the mean with design matrix information
    if design is not None:
        tmp = np.ndarray.copy(design)
        tmp[0:n_batch] = 0
        stand_mean = stand_mean + \
            np.transpose(np.dot(np.transpose(tmp), B_hat))
    return stand_mean


def standardise_data(dat: np.array, stand_mean: np.array, var_pooled: np.array, n_array: int) -> np.array:
    """
    standardise the data: substract mean and divide by variance

    Args:
        dat (np.array): data matrix
        stand_mean (np.array): standardised mean
        var_pooled (np.array): Variance for each gene and each batch
        n_array (int): total size of dataset

    Returns:
        np.array: standardised data matrix
    """
    standardised_data = (dat - stand_mean) / \
        np.dot(np.transpose(np.mat(np.sqrt(var_pooled))), np.mat([1]*n_array))
    return standardised_data


def fit_model(design: np.array, n_batch: int, standardised_data: np.array, batches: List[list], mean_only: bool, par_prior: bool, precision: float, ref_batch: int, ref: List[int], NAs: bool) -> Tuple[np.array, np.array, np.array]:
    """
    Fitting L/S model and finding priors.

    Args:
        design (np.array): model matrix for all covariates, including batch
        n_batch (int): number of batches
        standardised_data (np.array): standardised data matrix
        batches (List[list]): list of list of batches' elements
        mean_only (bool): True if just adjusting the means and not individual batch effects
        par_prior (bool): False for non-parametric estimation of batch effects
        precision (float): level of precision for precision computing
        ref_batch (int): reference batch selected
        ref (List[int]): the corresponding positions of the reference batch in the batch list
        NAs (bool): boolean characterising the presence of NaNs in the data matrix

    Returns:
        Tuple[np.array, np.array, np.array]:
            - estimated additive batch effect
            - estimated multiplicative batch effect
            - information about batches in design matrix
    """
    print("Fitting L/S model and finding priors.")

    # fraction of design matrix related to batches
    batch_design = design[0:n_batch]

    if not NAs:  # CF SUPRA FOR NAs
        # gamma_hat is the vector of additive batch effect
        gamma_hat = np.linalg.solve(np.dot(batch_design, np.transpose(batch_design)),
                                    np.dot(batch_design, np.transpose(standardised_data)))

    delta_hat = []  # delta_hat is the vector of estimated multiplicative batch effect

    if (mean_only):
        # no variance if mean_only == True
        delta_hat = [np.asarray([1]*len(standardised_data))] * len(batches)
    else:
        for i in batches:  # feed incrementally delta_hat
            list_map = np.transpose(np.transpose(standardised_data)[i]).var(
                axis=1)  # variance for each row
            delta_hat.append(np.squeeze(np.asarray(list_map)))

    gamma_bar = list(map(np.mean, gamma_hat))  # vector of means for gamma_hat
    t2 = list(map(np.var, gamma_hat))  # vector of variances for gamma_hat

    # calculates hyper priors for gamma (additive batch effect)
    a_prior = list(
        map(partial(compute_prior, 'a', mean_only=mean_only), delta_hat))
    b_prior = list(
        map(partial(compute_prior, 'b', mean_only=mean_only), delta_hat))

    # initialise gamma and delta for parameters estimation
    gamma_star = np.empty((n_batch, len(standardised_data)))
    delta_star = np.empty((n_batch, len(standardised_data)))

    if par_prior:
        # use param_fun function for parametric adjustments (cf. function definition)
        print("Finding parametric adjustments.")
        results = list(map(partial(param_fun,
                                   standardised_data=standardised_data,
                                   batches=batches,
                                   mean_only=mean_only,
                                   gamma_hat=gamma_hat,
                                   gamma_bar=gamma_bar,
                                   delta_hat=delta_hat,
                                   t2=t2,
                                   a_prior=a_prior,
                                   b_prior=b_prior), range(n_batch)))
    else:
        # use nonparam_fun for non-parametric adjustments (cf. function definition)
        print("Finding nonparametric adjustments")
        results = list(map(partial(nonparam_fun, mean_only=mean_only, delta_hat=delta_hat,
                                   standardised_data=standardised_data, batches=batches, gamma_hat=gamma_hat, precision=precision), range(n_batch)))

    for i in range(n_batch):  # store the results in gamma/delta_star
        results_i = results[i]
        gamma_star[i], delta_star[i] = results_i[0], results_i[1]

    # update if reference batch (the reference batch is not supposed to be modified)
    if ref_batch:
        len_gamma_star_ref = len(gamma_star[ref])
        gamma_star[ref] = [0] * len_gamma_star_ref
        delta_star[ref] = [1] * len_gamma_star_ref

    return (gamma_star, delta_star, batch_design)