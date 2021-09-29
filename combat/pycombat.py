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


from __future__ import annotations

import numpy as np
import pandas as pd

from sklearn.base import BaseEstimator, TransformerMixin

from .utils.combat_utils import (
    model_matrix,
    check_ref_batch,
    treat_covariates,
    treat_batches,
    calculate_mean_var,
    calculate_stand_mean,
    standardise_data,
    fit_model,
)
from .utils.common_utils import (
    check_NAs,
    check_mean_only,
)


class PyCombat(BaseEstimator, TransformerMixin):
    def __init__(self):
        """
        Corrects batch effect in microarray expression data.
        Takes an gene expression file and a list of known batches corresponding to each sample.
        """
        self.ref = None

        self.ref_batch = None

        self.gamma_star = None
        self.delta_star = None
        self.batch_design = None
        self.n_batches = None
        self.var_pooled = None
        self.stand_mean = None
        self.n_array = None
        self.ref_batch = None
        self.ref = None
        self.batches = None

    def check_if_fitted(self):
        error_message = "is not initialized, run .fit_transform() or .fit()"

        assert (
            self.gamma_star is not None
        ), f"gamma_star (estimated additive batch effect) {error_message}"
        assert (
            self.delta_star is not None
        ), f"delta_star (estimated multiplicative batch effect) {error_message}"
        assert (
            self.batch_design is not None
        ), f"batch_design (information about batches in design matrix) {error_message}"
        assert (
            self.n_batches is not None
        ), f"n_batches (list of batches lengths ) {error_message}"
        assert (
            self.var_pooled is not None
        ), f"var_pooled (standardised mean) {error_message}"
        assert (
            self.stand_mean is not None
        ), f"stand_mean (Variance for each gene and each batch) {error_message}"
        assert (
            self.n_array is not None
        ), f"n_array (total size of dataset) {error_message}"
        assert (
            self.batches is not None
        ), f"batches (list of unique batches) {error_message}"

    def fit(
        self,
        data: pd.DataFrame,
        batch: list,
        mod: list = [],
        par_prior: bool = True,
        prior_plots: bool = False,
        mean_only: bool = False,
        ref_batch: int = None,
        precision: float = None,
    ) -> PyCombat:
        """
        Initialized the important variables for the batch effect correction

        Args:
            data (pd.DataFrame): The expression matrix (dataframe). It contains the information about the gene expression (rows) for each sample (columns).
            batch (list): List of batch indexes. The batch list describes the batch for each sample. The batches list has as many elements as the number of columns in the expression matrix.
            mod (list, optional): List (or list of lists) of covariate(s) indexes. The mod list describes the covariate(s) for each sample. Each mod list has as many elements as the number of columns in the expression matrix. Defaults to [].
            par_prior (bool, optional): False for non-parametric estimation of batch effects. Defaults to True.
            prior_plots (bool, optional): True if requires to plot the priors. Defaults to False.
            mean_only (bool, optional): True if just adjusting the means and not individual batch effects. Defaults to False.
            ref_batch (int, optional): reference batch selected. Defaults to None.
            precision (float, optional): level of precision for precision computing. Defaults to None.

        Raises:
            ValueError: Error will be trigerred if NaN is found in the data

        Returns:
            PyCombat: The fitted object with all needed components for adjusting the data
        """
        self.list_samples = data.columns
        self.list_genes = data.index
        self.data_value = data.values

        self.ref_batch = ref_batch

        check_mean_only(mean_only)

        # Generate model matrix
        batchmod = model_matrix(list(batch), intercept=False, drop_first=False)

        # self.ref (List[int]): the corresponding positions of the reference batch in the batch list
        self.ref, batchmod = check_ref_batch(self.ref_batch, batch, batchmod)

        # self.batches (List[int]): list of unique batches
        # self.n_batches (List[int]): list of batches lengths
        # self.n_array (int): total size of dataset
        n_batch, self.batches, self.n_batches, self.n_array = treat_batches(batch)

        design = treat_covariates(batchmod, mod, self.ref, n_batch)

        NAs = check_NAs(self.data_value)
        if not NAs:
            B_hat, grand_mean, self.var_pooled = calculate_mean_var(
                design,
                self.batches,
                self.ref,
                self.data_value,
                NAs,
                self.ref_batch,
                self.n_batches,
                n_batch,
                self.n_array,
            )
            self.stand_mean = calculate_stand_mean(
                grand_mean, self.n_array, design, n_batch, B_hat
            )
            self.standardised_data = standardise_data(
                self.data_value, self.stand_mean, self.var_pooled, self.n_array
            )

            # self.gamma_star (np.array): estimated additive batch effect
            # self.delta_star (np.array): estimated multiplicative batch effect
            # self.batch_design (np.array): information about batches in design matrix
            self.gamma_star, self.delta_star, self.batch_design = fit_model(
                design,
                n_batch,
                self.standardised_data,
                self.batches,
                mean_only,
                par_prior,
                precision,
                self.ref_batch,
                self.ref,
                NAs,
            )
        else:
            raise ValueError("NaN value is not accepted")

        return self

    def transform(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Adjust the data, correct the data given the estimated batch effects

        Args:
            data (pd.DataFrame): The expression matrix (dataframe). It contains the information about the gene expression (rows) for each sample (columns).

        Returns:
            pd.DataFrame: The expression dataframe adjusted given the batch effects.
        """
        # Check if the fit function has been initiated before transform
        self.check_if_fitted()

        print("Adjusting the Data")
        bayes_data = np.transpose(self.standardised_data)
        for i, batch_index in enumerate(
            self.batches
        ):  # for each batch, specific correction
            bayes_data[batch_index] = (
                bayes_data[batch_index]
                - np.dot(np.transpose(self.batch_design)[batch_index], self.gamma_star)
            ) / np.transpose(
                np.outer(
                    np.sqrt(self.delta_star[i]), np.asarray([1] * self.n_batches[i])
                )
            )

        # renormalise the data after correction:
        # 1. multiply by variance
        # 2. add mean
        bayes_data = (
            np.multiply(
                np.transpose(bayes_data),
                np.outer(np.sqrt(self.var_pooled), np.asarray([1] * self.n_array)),
            )
            + self.stand_mean
        )

        # correction for reference batch
        if self.ref_batch:
            bayes_data[self.batches[self.ref]] = self.data_value[self.batches[self.ref]]

        # returns the data corrected for batch effects
        return pd.DataFrame(
            bayes_data, columns=self.list_samples, index=self.list_genes
        )


def pycombat(
    data: pd.DataFrame,
    batch: list,
    mod: list = [],
    par_prior: bool = True,
    prior_plots: bool = False,
    mean_only: bool = False,
    ref_batch: int = None,
    precision: float = None,
    **kwargs,
):
    return PyCombat().fit_transform(
        data,
        batch,
        mod=mod,
        par_prior=par_prior,
        prior_plots=prior_plots,
        mean_only=mean_only,
        ref_batch=ref_batch,
        precision=precision,
    )
