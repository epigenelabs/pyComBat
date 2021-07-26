import pytest

import numpy as np
import pandas as pd

from combat.utils.combat_utils import model_matrix, compute_prior, postmean, postvar, it_sol, int_eprior, check_ref_batch, treat_batches, treat_covariates, calculate_mean_var, calculate_stand_mean, standardise_data, fit_model
from combat.utils.common_utils import check_NAs

"""
Load the shared variables
"""
batch = np.asarray([1, 1, 1, 2, 2, 3, 3, 3, 3])
matrix = np.transpose([np.random.normal(size=1000, loc=3, scale=1), np.random.normal(size=1000, loc=3, scale=1), np.random.normal(size=1000, loc=3, scale=1),
                      np.random.normal(size=1000, loc=2, scale=0.6), np.random.normal(size=1000, loc=2, scale=0.6),
                      np.random.normal(size=1000, loc=4, scale=1), np.random.normal(size=1000, loc=4, scale=1), np.random.normal(size=1000, loc=4, scale=1), np.random.normal(size=1000, loc=4, scale=1)])

matrix = pd.DataFrame(data=matrix,columns=["sample_"+str(i+1) for i in range(9)],index=["gene_"+str(i+1) for i in range(1000)])
ref_batch = None
mean_only = False
par_prior = False
precision = None
mod = []
dat = matrix.values
batchmod = model_matrix(list(batch), intercept=False, drop_first=False)
ref,batchmod = check_ref_batch(ref_batch,batch,batchmod)
n_batch, batches, n_batches, n_array = treat_batches(batch)
design = treat_covariates(batchmod, mod, ref, n_batch)
NAs = check_NAs(dat)
B_hat, grand_mean, var_pooled = calculate_mean_var(design, batches, ref, dat, NAs, ref_batch, n_batches, n_batch, n_array)
stand_mean = calculate_stand_mean(grand_mean, n_array, design, n_batch,B_hat)
s_data = standardise_data(dat, stand_mean, var_pooled, n_array)
gamma_star, delta_star, batch_design = fit_model(design,n_batch,s_data, batches, mean_only, par_prior, precision, ref_batch, ref, NAs)


# test for compute_prior
def test_compute_prior():
    print("aprior", compute_prior("a", gamma_star, False))
    assert compute_prior("a", gamma_star, True) == 1
    print("bprior", compute_prior("b", gamma_star, False))
    assert compute_prior("b", gamma_star, True) == 1


# test for model_matrix
def test_model_matrix():
    model_matrix_test = model_matrix([1, 1, 0, 1, 0])
    assert np.shape(model_matrix_test) == (5, 2)
    assert list(model_matrix_test[0]) == [1.0, 1.0]



# test for compute_prior
def test_compute_prior():
    print("aprior", compute_prior("a", gamma_star, False))
    assert compute_prior("a", gamma_star, True) == 1
    print("bprior", compute_prior("b", gamma_star, False))
    assert compute_prior("b", gamma_star, True) == 1


# test for postmean
def test_postmean():
    assert np.shape(postmean(gamma_star, delta_star, gamma_star, delta_star)) == np.shape(gamma_star)


# test for postvar
def test_postvar():
    assert np.sum(postvar([2, 4, 6], 2, 1, 1) - [2, 3, 4]) == 0


# test for it_sol
def test_it_sol():
    assert 1 == 1


# test for int_eprior
def test_int_eprior():
    assert 1 == 1


# test for model_matrix
def test_model_matrix():
    model_matrix_test = model_matrix([1, 1, 0, 1, 0])
    assert np.shape(model_matrix_test) == (5, 2)
    assert list(model_matrix_test[0]) == [1.0, 1.0]


# test for check_ref_batch
def test_check_ref_batch():
    assert check_ref_batch(1, batch, batchmod) == (0, batchmod)
    assert check_ref_batch(2, batch, batchmod) == (1, batchmod)
    print("Using batch 1 then 2. Above lines should inform on that.")
    assert check_ref_batch(None, batch, batchmod) == (None, batchmod)


# test for treat_batches 
def test_treat_batches():
    assert n_batch == 3
    assert batches[0].tolist() == [0, 1, 2]
    assert batches[1].tolist() == [3, 4]
    assert batches[2].tolist() == [5, 6, 7, 8]
    assert n_batches == [3, 2, 4]
    assert n_array == 9
 

# test for treat_covariates
def test_treat_covariates():
    batchmod = model_matrix(list(batch), intercept=False, drop_first=False)
    assert np.sum(design - np.transpose(batchmod)) == 0


# test for calculate_mean_var
def test_calculate_mean_var():
    assert np.shape(B_hat)[0] == np.shape(design)[0]
    assert np.shape(grand_mean)[0] == np.shape(dat)[0]
    assert np.shape(var_pooled)[0] == np.shape(dat)[0]


# test for calculate_stand_mean
def test_calculate_stand_mean():
    assert np.shape(stand_mean) == np.shape(dat)


# test for standardise_data
def test_standardise_data():
    assert np.shape(s_data) == np.shape(dat)


# test for fit_model
def test_fit_model():
    assert np.shape(gamma_star)[1] == np.shape(dat)[0]
    assert np.shape(delta_star)[1] == np.shape(dat)[0]
    assert np.shape(batch_design) == np.shape(design)
