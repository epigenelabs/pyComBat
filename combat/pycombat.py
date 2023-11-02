#-----------------------------------------------------------------------------
# Copyright (C) 2019-2022 A. Behdenna, A. Nordor, J. Haziza and A. Gema

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
# date 	Sept 2020
#-----------------------------------------------------------------------------



import numpy as np
from functools import partial
import mpmath as mp
import pandas as pd

from combat.exceptions import ConfoundingVariablesError

#import unittest


def model_matrix(info, intercept=True, drop_first=True):
    """Creates the model_matrix from batch list

    Arguments:
        info {list} -- list info with batch or covariates data
        intercept {bool} -- boolean for intercept in model matrix

    Returns:
        matrix -- model matrix generate from batch list
    """
    if not isinstance(info[0], list):
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


def all_1(list_of_elements):
    """checks if all elements in a list are 1s

    Arguments:
        list_of_elements {list} -- list of elements

    Returns:
        bool -- True iff all elements of the list are 1s
    """
    return((list_of_elements == 1).all())


# aprior and bprior are useful to compute "hyper-prior values"
# -> prior parameters used to estimate the prior gamma distribution for multiplicative batch effect
# aprior - calculates empirical hyper-prior values

def compute_prior(prior, gamma_hat, mean_only):
    """[summary]

    Arguments:
        prior {char} -- 'a' or 'b' depending of the prior to be calculated
        gamma_hat {matrix} -- matrix of additive batch effect
        mean_only {bool} -- True iff mean_only selected

    Returns:
        float -- [the prior calculated (aprior or bprior)
    """
    if mean_only:
        return 1
    m = np.mean(gamma_hat)
    s2 = np.var(gamma_hat)
    if prior == 'a':
        return (2*s2+m*m)/s2
    elif prior == 'b':
        return (m*s2+m*m*m)/s2


def postmean(g_bar, d_star, t2_n, t2_n_g_hat):
    """estimates additive batch effect

    Arguments:
        g_bar {matrix} -- additive batch effect
        d_star {matrix} -- multiplicative batch effect
        t2_n {matrix} --
        t2_n_g_hat {matrix} --

    Returns:
        matrix -- estimated additive batch effect
    """
    return np.divide(t2_n_g_hat+d_star*g_bar, np.asarray(t2_n+d_star))


def postvar(sum2, n, a, b):
    """estimates multiplicative batch effect

    Arguments:
        sum2 {vector} --
        n {[type]} --
        a {float} -- aprior
        b {float} -- bprior

    Returns:
        matrix -- estimated multiplicative batch effect
    """
    return(np.divide((np.multiply(0.5, sum2)+b), (np.multiply(0.5, n)+a-1)))


def it_sol(sdat, g_hat, d_hat, g_bar, t2, a, b, conv=0.0001, exit_iteration=10e5):
    """iterative solution for Empirical Bayesian method

    Arguments:
        sdat {matrix} --
        g_hat {matrix} -- average additive batch effect
        d_hat {matrix} -- average multiplicative batch effect
        g_bar {matrix} -- additive batch effect
        t2 {matrix} --
        a {float} -- aprior
        b {float} -- bprior

    Keyword Arguments:
        conv {float} -- convergence criterion (default: {0.0001})
        exit_iteration {float} -- maximum number of iterations before exit (default: {10e5})

    Returns:
        array list -- estimated additive and multiplicative batch effect
    """
    n = [len(i) for i in sdat]
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
    adjust = np.asarray([g_new, d_new], dtype=object)
    return(adjust)  # remove parenthesis in returns

# int_eprior - Monte Carlo integration function to find nonparametric adjustments
# Johnson et al (Biostatistics 2007, supp.mat.) show that we can estimate the multiplicative and additive batch effects with an integral
# This integral is numerically computed through Monte Carlo inegration (iterative method)


def int_eprior(sdat, g_hat, d_hat, precision):
    """ int_eprior - Monte Carlo integration function to find nonparametric adjustments
        Johnson et al (Biostatistics 2007, supp.mat.) show that we can estimate the multiplicative and additive batch effects with an integral
        This integral is numerically computed through Monte Carlo inegration (iterative method)

    Arguments:
        sdat {matrix} -- data matrix
        g_hat {matrix} -- average additive batch effect
        d_hat {matrix} -- average multiplicative batch effect
        precision {float} -- level of precision for precision computing

    Returns:
        array list -- estimated additive and multiplicative batch effect
    """
    g_star = []
    d_star = []
    # use this variable to only print error message once if approximation used
    test_approximation = 0
    for i in range(len(sdat)):
        # additive batch effect
        g = np.delete(np.transpose(g_hat), i)
        # multiplicative batch effect
        d = np.asarray(np.delete(np.transpose(d_hat), i))
        x = np.expand_dims(np.transpose(sdat[i]), axis=1)
        n = len(x)
        j = [1]*n
        dat = np.repeat(x, len(np.transpose(g)), axis=1)
        resid2 = np.square(dat-g)
        sum2 = np.dot(np.transpose(resid2), j)
        # /begin{handling high precision computing}
        temp_2d = 2*d
        if (precision == None):
            LH = np.power(1/(np.pi*temp_2d), n/2)*np.exp(np.negative(sum2)/(temp_2d))

        else:  # only if precision parameter informed
            # increase the precision of the computing (if negative exponential too close to 0)
            mp.dps = precision
            buf_exp = np.array(list(map(mp.exp, np.negative(sum2)/(temp_2d))))
            buf_pow = np.array(list(map(partial(mp.power, y=n/2), 1/(np.pi*temp_2d))))
            #print(buf_exp.dtype, buf_pow.dtype)
            LH = buf_pow*buf_exp  # likelihood
        # /end{handling high precision computing}
        LH = np.nan_to_num(LH)  # corrects NaNs in likelihood
        if np.sum(LH) == 0 and test_approximation == 0:
            test_approximation = 1  # this message won't appear again
            print("###\nValues too small, approximation applied to avoid division by 0.\nPrecision mode can correct this problem, but increases computation time.\n###")

        if np.sum(LH) == 0: # correction for LH full of 0.0
            LH[LH == 0] = np.exp(-745)
            g_star.append(np.sum(g*LH)/np.sum(LH))
            d_star.append(np.sum(d*LH)/np.sum(LH))
        else:
            g_star.append(np.sum(g*LH)/np.sum(LH))
            d_star.append(np.sum(d*LH)/np.sum(LH))
    adjust = np.asarray([np.asarray(g_star), np.asarray(d_star)])
    return(adjust)


def param_fun(i, s_data, batches, mean_only, gamma_hat, gamma_bar, delta_hat, t2, a_prior, b_prior):
    """parametric estimation of batch effects

    Arguments:
        i {int} -- column index
        s_data {matrix} --
        batches {list list} -- list of list of batches' elements
        mean_only {bool} -- True iff mean_only selected
        gamma_hat {matrix} -- average additive batch effect
        gamma_bar {matrix} -- estimated additive batch effect
        delta_hat {matrix} -- average multiplicative batch effect
        t2 {matrix} --
        a_prior {float} -- aprior
        b_prior {float} -- bprior

    Returns:
        array list -- estimated adjusted additive and multiplicative batch effect
    """
    g_hat = np.expand_dims(gamma_hat[i], axis=0)
    if mean_only:  # if mean_only, no need for complex method: batch effect is immediately calculated
        t2_n = np.multiply(t2[i], 1)
        t2_n_g_hat = np.multiply(t2_n, g_hat)
        gamma_star = postmean(gamma_bar[i], 1, t2_n, t2_n_g_hat)  # additive batch effect
        delta_star = [1]*len(s_data)  # multiplicative batch effect
    else:  # if not(mean_only) then use it_solve
        temp = it_sol(np.transpose(np.transpose(s_data)[
                      batches[i]]), g_hat, delta_hat[i], gamma_bar[i], t2[i], a_prior[i], b_prior[i])
        gamma_star = temp[0]  # additive batch effect
        delta_star = temp[1]  # multiplicative batch effect
    return [gamma_star, delta_star]

def nonparam_fun(i, mean_only, delta_hat, s_data, batches, gamma_hat, precision):
    """non-parametric estimation

    Arguments:
        i {int} -- column index
        mean_only {bool} -- True iff mean_only selected
        delta_hat {matrix} -- estimated multiplicative batch effect
        s_data {matrix} --
        batches {list list} -- list of list of batches' elements
        gamma_hat {matrix} -- estimated additive batch effect
        precision {float} -- level of precision for precision computing

    Returns:
        array list -- estimated adjusted additive and multiplicative batch effect
    """
    if mean_only:  # if mean only, change delta_hat to vector of 1s
        delta_hat[i] = [1]*len(delta_hat[i])
    # use int_eprior for non-parametric estimation
    temp = int_eprior(np.transpose(np.transpose(s_data)[
                      batches[i]]), np.expand_dims(gamma_hat[i], axis=0), delta_hat[i], precision)
    return [temp[0], temp[1]]

############
# pyComBat #
############


def check_mean_only(mean_only):
    """checks mean_only option

    Arguments:
        mean_only {boolean} -- user's choice about mean_only

    Returns:
        ()
    """
    if mean_only == True:
        print("Using mean only version")


def define_batchmod(batch):
    """generates model matrix

    Arguments:
        batch {list} -- list of batch id

    Returns:
        batchmod {matrix} -- model matrix for batches
    """
    batchmod = model_matrix(list(batch), intercept=False, drop_first=False)
    return(batchmod)


def check_ref_batch(ref_batch, batch, batchmod):
    """check ref_batch option and treat it if needed

    Arguments:
        ref_batch -- the reference batch
        batch {list} -- list of batch id
        batchmod {matrix} -- model matrix related to batches

    Returns:
        ref {int} -- the index of the reference batch in the batch list
        batchmod {matrix} -- updated model matrix related to batches, with reference
    """
    if ref_batch is not None:
        if ref_batch not in batch:
            raise ValueError("Reference level ref.batch must be one of the levels of batch.")
        print("Using batch "+str(ref_batch) +
              " as a reference batch.")
        # ref keeps in memory the columns concerned by the reference batch
        ref = np.where(np.unique(batch) == ref_batch)[0][0]
        # updates batchmod with reference
        batchmod[:,ref] = 1
    else:
        ref = None  # default settings
    return(ref, batchmod)


def treat_batches(batch):
    """treat batches

    Arguments:
        batch {list} -- batch list

    Returns:
        n_batch {int} -- number of batches
        batches {int list} -- list of unique batches
        n_batches {int list} -- list of batches lengths
        n_array {int} -- total size of dataset
    """
    batch = pd.Series(batch)
    n_batch = len(np.unique(batch))  # number of batches
    print("Found "+str(n_batch)+" batches.")
    batches = []  # list of lists, contains the list of position for each batch
    for i in range(n_batch):
        batches.append(np.where(batch == np.unique(batch)[i])[0].astype(np.int32))
    n_batches = list(map(len, batches))
    if 1 in n_batches:
        #mean_only = True  # no variance if only one sample in a batch - mean_only has to be used
        print("\nOne batch has only one sample, try setting mean_only=True.\n")
    n_array = sum(n_batches)
    return(n_batch, batches, n_batches, n_array)


def treat_covariates(batchmod, mod, ref, n_batch):
    """treat covariates

    Arguments:
        batchmod {matrix} -- model matrix for batch
        mod {matrix} -- model matrix for other covariates
        ref {int} -- index of reference batch
        n_batch {int} -- number of batches

    Returns:
        check {bool list} -- a list characterising all covariates
        design {matrix} -- model matrix for all covariates, including batch
    """
    # design matrix for sample conditions
    if mod == []:
        design = batchmod
    else:
        mod_matrix = model_matrix(mod, intercept=True)
        design = np.concatenate((batchmod, mod_matrix), axis=1)
    check = list(map(all_1, np.transpose(design)))
    if ref is not None:  # if ref
        check[ref] = False  # the reference batch is not considered a covariate
    design = design[:, ~np.array(check)]
    design = np.transpose(design)

    print("Adjusting for "+str(len(design)-len(np.transpose(batchmod))) +
          " covariate(s) or covariate level(s).")

    # if matrix cannot be invertible, different cases
    if np.linalg.matrix_rank(design) < len(design):
        if len(design) == n_batch + 1:  # case 1: covariate confunded with a batch
            raise ConfoundingVariablesError("Covariate is confounded with batch. Try removing the covariates.")
        if len(design) > n_batch + 1:  # case 2: multiple covariates confunded with a batch
            if np.linalg.matrix_rank(np.transpose(design)[:n_batch]) < len(design):
                raise ConfoundingVariablesError("Confounded design. Try removing one or more covariates.")
            else:  # case 3: at least a covariate confunded with a batch
                raise ConfoundingVariablesError("At least one covariate is confounded with batch. Try removing confounded covariates.")
    return(design)


def check_NAs(dat):
    """check if NaNs - in theory, we construct the data without NAs

    Arguments:
        dat {matrix} -- the data matrix

    Returns:
        NAs {bool} -- boolean characterising the presence of NaNs in the data matrix
    """
    # NAs = True in (np.isnan(dat))
    NAs = np.isnan(np.sum(dat))  # Check if NaN exists
    if NAs:
        print("Found missing data values. Please remove all missing values before proceeding with pyComBat.")
    return(NAs)


def calculate_mean_var(design, batches, ref, dat, NAs, n_batches, n_batch, n_array):
    """ calculates the Normalisation factors

    Arguments:
        design {matrix} -- model matrix for all covariates
        batches {int list} -- list of unique batches
        ref {int} -- reference batch index
        dat {matrix} -- data matrix
        NAs {bool} -- presence of NaNs in the data matrix
        n_batches {int list} -- list of batches lengths
        n_array {int} -- total size of dataset

    Returns:
        B_hat {matrix} -- regression coefficients corresponding to the design matrix
        grand_mean {matrix} -- Mean for each gene and each batch
        var_pooled {matrix} -- Variance for each gene and each batch
    """
    print("Standardizing Data across genes.")
    if not(NAs):  # NAs not supported
        # B_hat is the vector of regression coefficients corresponding to the design matrix
        B_hat = np.linalg.solve(np.dot(design, np.transpose(
            design)), np.dot(design, np.transpose(dat)))

    # Calculates the general mean
    if ref is not None:
        grand_mean = np.transpose(B_hat[ref])
    else:
        grand_mean = np.dot(np.transpose(
            [i / n_array for i in n_batches]), B_hat[0:n_batch])
    # Calculates the general variance
    if not NAs:  # NAs not supported
        if ref is not None:  # depending on ref batch
            ref_dat = np.transpose(np.transpose(dat)[batches[ref]])
            var_pooled = np.dot(np.square(ref_dat - np.transpose(np.dot(np.transpose(
                design)[batches[ref]], B_hat))), [1/n_batches[ref]]*n_batches[ref])
        else:
            var_pooled = np.dot(np.square(
                dat - np.transpose(np.dot(np.transpose(design), B_hat))), [1/n_array]*n_array)

    return(B_hat, grand_mean, var_pooled)


def calculate_stand_mean(grand_mean, n_array, design, n_batch, B_hat):
    """ transform the format of the mean for substraction

    Arguments:
        grand_mean {matrix} -- Mean for each gene and each batch
        n_array {int} -- total size of dataset
        design {[type]} -- design matrix for all covariates including batch
        n_batch {int} -- number of batches
        B_hat {matrix} -- regression coefficients corresponding to the design matrix

    Returns:
        stand_mean {matrix} -- standardised mean
    """
    stand_mean = np.dot(np.transpose(np.expand_dims(grand_mean, axis=0)), np.array([[1]*n_array]))
    # corrects the mean with design matrix information
    if design is not None:
        tmp = np.ndarray.copy(design)
        tmp[0:n_batch] = 0
        stand_mean = stand_mean + \
            np.transpose(np.dot(np.transpose(tmp), B_hat))
    return(stand_mean)


def standardise_data(dat, stand_mean, var_pooled, n_array):
    """standardise the data: substract mean and divide by variance

    Arguments:
        dat {matrix} -- data matrix
        stand_mean {matrix} -- standardised mean
        var_pooled {matrix} -- Variance for each gene and each batch
        n_array {int} -- total size of dataset

    Returns:
        s_data {matrix} -- standardised data matrix
    """
    s_data = np.divide((dat - stand_mean), \
        np.dot(np.transpose(np.expand_dims(np.sqrt(var_pooled), axis=0)), np.array([[1]*n_array])))
    return(s_data)


def fit_model(design, n_batch, s_data, batches, mean_only, par_prior, precision, ref, NAs):
    print("Fitting L/S model and finding priors.")

    # fraction of design matrix related to batches
    batch_design = design[0:n_batch]

    if not NAs:  # CF SUPRA FOR NAs
        # gamma_hat is the vector of additive batch effect
        gamma_hat = np.linalg.solve(np.dot(batch_design, np.transpose(batch_design)),
                                    np.dot(batch_design, np.transpose(s_data)))

    delta_hat = []  # delta_hat is the vector of estimated multiplicative batch effect

    if (mean_only):
        # no variance if mean_only == True
        delta_hat = [np.asarray([1]*len(s_data))] * len(batches)
    else:
        for i in batches:  # feed incrementally delta_hat
            list_map = np.transpose(np.transpose(s_data)[i]).var(
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
    gamma_star = np.empty((n_batch, len(s_data)))
    delta_star = np.empty((n_batch, len(s_data)))

    if par_prior:
        # use param_fun function for parametric adjustments (cf. function definition)
        print("Finding parametric adjustments.")
        results = list(map(partial(param_fun,
                                   s_data=s_data,
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
                                   s_data=s_data, batches=batches, gamma_hat=gamma_hat, precision=precision), range(n_batch)))

    for i in range(n_batch):  # store the results in gamma/delta_star
        results_i = results[i]
        gamma_star[i], delta_star[i] = results_i[0], results_i[1]

    # update if reference batch (the reference batch is not supposed to be modified)
    if ref is not None:
        len_gamma_star_ref = len(gamma_star[ref])
        gamma_star[ref] = [0] * len_gamma_star_ref
        delta_star[ref] = [1] * len_gamma_star_ref

    return(gamma_star, delta_star, batch_design)


def adjust_data(s_data, gamma_star, delta_star, batch_design, n_batches, var_pooled, stand_mean, n_array, ref, batches, dat):
    """Adjust the data -- corrects for estimated batch effects

    Arguments:
        s_data {matrix} -- standardised data matrix
        gamma_star {matrix} -- estimated additive batch effect
        delta_star {matrix} -- estimated multiplicative batch effect
        batch_design {matrix} -- information about batches in design matrix
        n_batches {int list} -- list of batches lengths
        stand_mean {matrix} -- standardised mean
        var_pooled {matrix} -- Variance for each gene and each batch
        n_array {int} -- total size of dataset
        ref {int} -- the index of the reference batch in the batch list
        batches {int list} -- list of unique batches
        dat

    Returns:
        bayes_data [matrix] -- data adjusted for correction of batch effects
    """
    # Now we adjust the data:
    # 1. substract additive batch effect (gamma_star)
    # 2. divide by multiplicative batch effect (delta_star)
    print("Adjusting the Data")
    bayes_data = np.transpose(s_data)
    j = 0
    for i in batches:  # for each batch, specific correction
        bayes_data[i] = (bayes_data[i] - np.dot(np.transpose(batch_design)[i], gamma_star)) / \
            np.transpose(
                np.outer(np.sqrt(delta_star[j]), np.asarray([1]*n_batches[j])))
        j += 1

    # renormalise the data after correction:
    # 1. multiply by variance
    # 2. add mean
    bayes_data = np.multiply(np.transpose(bayes_data), np.outer(
        np.sqrt(var_pooled), np.asarray([1]*n_array))) + stand_mean

    # correction for reference batch
    if ref is not None:
        bayes_data[batches[ref]] = dat[batches[ref]]

    # returns the data corrected for batch effects
    return bayes_data


def pycombat(data, batch, mod=[], par_prior=True, prior_plots=False, mean_only=False, ref_batch=None, precision=None, **kwargs):
    """Corrects batch effect in microarray expression data. Takes an gene expression file and a list of known batches corresponding to each sample.

    Arguments:
        data {matrix} -- The expression matrix (dataframe). It contains the information about the gene expression (rows) for each sample (columns).

        batch {list} -- List of batch indexes. The batch list describes the batch for each sample. The batches list has as many elements as the number of columns in the expression matrix.

    Keyword Arguments:
        mod {list} -- List (or list of lists) of covariate(s) indexes. The mod list describes the covariate(s) for each sample. Each mod list has as many elements as the number of columns in the expression matrix (default: {[]}).

        par_prior {bool} -- False for non-parametric estimation of batch effects (default: {True}).

        prior_plots {bool} -- True if requires to plot the priors (default: {False} -- Not implemented yet!).

        mean_only {bool} -- True iff just adjusting the means and not individual batch effects (default: {False}).

        ref_batch -- reference batch selected (default: {None}).

        precision {float} -- level of precision for precision computing (default: {None}).

    Returns:
        bayes_data_df -- The expression dataframe adjusted for batch effects.
    """

    list_samples = data.columns
    list_genes = data.index
    dat = data.values

    check_mean_only(mean_only)

    batchmod = define_batchmod(batch)
    ref, batchmod = check_ref_batch(ref_batch, batch, batchmod)
    n_batch, batches, n_batches, n_array = treat_batches(batch)
    design = treat_covariates(batchmod, mod, ref, n_batch)
    NAs = check_NAs(dat)
    if not(NAs):
        B_hat, grand_mean, var_pooled = calculate_mean_var(
            design, batches, ref, dat, NAs, n_batches, n_batch, n_array)
        stand_mean = calculate_stand_mean(
            grand_mean, n_array, design, n_batch, B_hat)
        s_data = standardise_data(dat, stand_mean, var_pooled, n_array)
        gamma_star, delta_star, batch_design = fit_model(
            design, n_batch, s_data, batches, mean_only, par_prior, precision, ref, NAs)
        bayes_data = adjust_data(s_data, gamma_star, delta_star, batch_design,
                                n_batches, var_pooled, stand_mean, n_array, ref, batches, dat)

        bayes_data_df = pd.DataFrame(bayes_data,
                    columns = list_samples,
                    index = list_genes)

        return(bayes_data_df)
    else:
        raise ValueError("NaN value is not accepted")
