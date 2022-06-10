# pyComBat releases changelog

## [0.1.7] - 18-05-2020

* First release

## [0.2.0] - 02-09-2020

pyComBat 0.2.0

* Batches can now be defined with non-integers lists (e.g. strings)

* Covariates can now be taken into account for the batch-effects correction through the use of a model matrix

* Corrected a bug forbidding the first batch to be a reference batch

* Added unit testing

## [0.2.1] - 12-02-2021

pyComBat 0.2.1

* Covariates can now be defined by a simple list instead of a model matrix

* Performances enhancement

## [0.3.0] - 12-03-2021

pyComBat 0.3.0

* Corrected a bug preventing the use of multiple covariates at the same time

* Corrected a bug preventing the use of batch 0 as a reference batch

* pyComBat now exits with an informative message when nans are detected

* Removed soon deprecated functions calls

## [0.3.1] - 16-03-2022

pyComBat 0.3.1

* Updated copyright

* Updated dependency versions

* Corrected a bug preventing the definition of batches with lists of strings

* Corrected typos in the documentation

## [0.3.2] - 10-06-2022

pyComBat 0.3.2

* Enhanced error handling for confounded variables
