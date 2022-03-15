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




# this file is only used for unit testing
# We import the function that will be tested one by one, incrementally
# Generates a report about the functions tested

import numpy as np
import pandas as pd
#from patsy import dmatrix

##########
# import function used for unit testing
from .pycombat import model_matrix, all_1
#covariate_model_matrix
from .pycombat import compute_prior, postmean, postvar, it_sol, int_eprior
from .pycombat import check_mean_only, define_batchmod, check_ref_batch, treat_batches, treat_covariates, check_NAs
from .pycombat import calculate_mean_var, calculate_stand_mean
from .pycombat import standardise_data, fit_model, adjust_data
from .pycombat import pycombat

##########
print("\n#### Unit Testing for pyComBat ####\n")

##########
# Define constants for unit testing

batch = np.asarray([1,1,1,2,2,3,3,3,3])
# matrix = np.transpose(np.asmatrix([np.random.normal(size=1000,loc=3,scale=1),np.random.normal(size=1000,loc=3,scale=1),np.random.normal(size=1000,loc=3,scale=1),
#                       np.random.normal(size=1000,loc=2,scale=0.6),np.random.normal(size=1000,loc=2,scale=0.6),
#                       np.random.normal(size=1000,loc=4,scale=1),np.random.normal(size=1000,loc=4,scale=1),np.random.normal(size=1000,loc=4,scale=1),np.random.normal(size=1000,loc=4,scale=1)]))
matrix = np.transpose([np.random.normal(size=1000,loc=3,scale=1),np.random.normal(size=1000,loc=3,scale=1),np.random.normal(size=1000,loc=3,scale=1),
                      np.random.normal(size=1000,loc=2,scale=0.6),np.random.normal(size=1000,loc=2,scale=0.6),
                      np.random.normal(size=1000,loc=4,scale=1),np.random.normal(size=1000,loc=4,scale=1),np.random.normal(size=1000,loc=4,scale=1),np.random.normal(size=1000,loc=4,scale=1)])

matrix = pd.DataFrame(data=matrix,columns=["sample_"+str(i+1) for i in range(9)],index=["gene_"+str(i+1) for i in range(1000)])

print("Matrix and batch generated.")

matrix_adjusted = pycombat(matrix,batch)

list_samples = matrix_adjusted.columns
list_genes = matrix_adjusted.index

print("Adjusted matrix generated.")

##########
# local tests before unit testing

def test_means():
    print("mean matrix:",np.mean(matrix))
    print("mean matrix (adjusted):",np.mean(matrix_adjusted))
    print("**********")
    print("var matrix:",np.var(matrix))
    print("var matrix (adjusted):",np.var(matrix_adjusted))

###########
# useful constants for unit testing
ref_batch = None
mean_only = False
par_prior = False
precision = None
mod = []
dat = matrix.values
#batchmod = define_batchmod(batch)
batchmod = model_matrix(list(batch), intercept=False, drop_first=False)
ref,batchmod = check_ref_batch(ref_batch,batch,batchmod)
n_batch, batches, n_batches, n_array = treat_batches(batch)
design = treat_covariates(batchmod, mod, ref, n_batch)
NAs = check_NAs(dat)
B_hat, grand_mean, var_pooled = calculate_mean_var(design, batches, ref, dat, NAs, ref_batch, n_batches, n_batch, n_array)
stand_mean = calculate_stand_mean(grand_mean, n_array, design, n_batch,B_hat)
s_data = standardise_data(dat, stand_mean, var_pooled, n_array)
gamma_star, delta_star, batch_design = fit_model(design,n_batch,s_data, batches, mean_only, par_prior, precision, ref_batch, ref, NAs)
bayes_data = adjust_data(s_data, gamma_star, delta_star, batch_design, n_batches, var_pooled, stand_mean, n_array, ref_batch, ref, batches, dat)


##########
# Unit tests

# test for compute_prior
def test_compute_prior():
    print("aprior",compute_prior("a",gamma_star,False))
    assert compute_prior("a",gamma_star,True) == 1
    print("bprior",compute_prior("b",gamma_star,False))
    assert compute_prior("b",gamma_star,True) == 1

# test for postmean
def test_postmean():
    assert np.shape(postmean(gamma_star,delta_star,gamma_star,delta_star)) == np.shape(gamma_star)

# test for postvar
def test_postvar():
    assert np.sum(postvar([2,4,6],2,1,1) - [2,3,4]) == 0

# test for it_sol
def test_it_sol():
    ()

# test for int_eprior
def test_int_eprior():
    ()

# test for model_matrix
def test_model_matrix():
    model_matrix_test = model_matrix([1,1,0,1,0])
    assert np.shape(model_matrix_test) == (5,2)
    assert list(model_matrix_test[0]) == [1.0,1.0]

# tests for all_1 function
def test_all_1():
    assert all_1(np.array([1,1,1,1,1])) == True

    assert all_1(np.array([1,1,1,1,0])) == False
    assert all_1(np.array([0,0,0,0,0])) == False

    assert all_1(np.array([1.5,0.5,1,1,1])) == False # This test to show the limit of the method we use


# test for check_mean_only
def test_check_mean_only():
    check_mean_only(True)
    check_mean_only(False)
    print("Only one line of text should have been printed above.")

# test for define_batchmode
def test_define_batchmod():
    assert np.shape(define_batchmod(batch)) == (9,3)

# test for check_ref_batch
def test_check_ref_batch():
    assert check_ref_batch(1,batch,batchmod) == (0, batchmod)
    assert check_ref_batch(2,batch,batchmod) == (1, batchmod)
    print("Using batch 1 then 2. Above lines should inform on that.")
    assert check_ref_batch(None,batch,batchmod) == (None, batchmod)

# test for treat_batches 
def test_treat_batches():
    assert n_batch == 3
    assert batches[0].tolist() == [0,1,2]
    assert batches[1].tolist() == [3,4]
    assert batches[2].tolist() == [5,6,7,8]
    assert n_batches == [3,2,4]
    assert n_array == 9
 
# test for treat_covariates
def test_treat_covariates():
    batchmod = define_batchmod(batch)
    assert np.sum(design - np.transpose(batchmod)) == 0

# test for check_NAs
def test_check_NAs():
    assert check_NAs([0,1,2]) == False

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

#test for adjust_data
def test_adjust_data():
    assert np.shape(bayes_data) == np.shape(dat)

# general tests on pyComBat
def test_pycombat():
    assert np.shape(matrix) == np.shape(matrix_adjusted)
    assert abs(np.mean(matrix.values)-np.mean(matrix_adjusted.values)) <= np.mean(matrix.values)*0.05
    assert np.var(matrix_adjusted.values) <= np.var(matrix.values)

if __name__ == '__main__':
    test_means()
    print("\n*** UNIT TESTING ***\n")
    # unittest.main()
