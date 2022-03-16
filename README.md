# pyComBat

pyComBat [1] is a Python 3 implementation of ComBat [2], one of the most widely used tool for correcting technical biases, called batch effects, in microarray expression data.

More detailed documentation can be found at [this address](https://epigenelabs.github.io/pyComBat/).

## TO DO

## Minimum dependencies

We list here the versions of the packages that have been used for development/testing of pyComBat, as well as for writing the documentation.

### pyComBat dependencies

* python 3.6

* numpy 1.18.5

* mpmath 1.1.0

* pandas 0.24.2

* patsy 0.5.1

### Documentation

* sphinx 2.1.2

## Usage example

### Installation

You can install pyComBat directly with:

```python
pip install combat
```

You can upgrade pyComBat to its latest version with:

```python
pip install combat --upgrade
```

### Running pyComBat

The simplest way of using pyComBat is to first import it, and then simply use the pycombat function with default parameters:

```python
from combat.pycombat import pycombat
data_corrected = pycombat(data,batch)
```

* data: The expression matrix as a dataframe. It contains the information about the gene expression (rows) for each sample (columns).

* batch: List of batch indexes. The batch list describes the batch for each sample. The list of batches contains as many elements as the number of columns in the expression matrix.

## How to contribute

Please refer to [CONTRIBUTING.md](https://github.com/epigenelabs/pyComBat/blob/master/CONTRIBUTING.md) to learn more about the contribution guidelines.

## References

[1] Behdenna A, Haziza J, Azencot CA and Nordor A. (2020) pyComBat, a Python tool for batch effects correction in high-throughput molecular data using empirical Bayes methods. bioRxiv doi: 10.1101/2020.03.17.995431

[2] Johnson W E, et al. (2007) Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics, 8, 118–127
