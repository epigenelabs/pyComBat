# pyComBat

pyComBat is a Python 3 implementation of ComBat (1), one of the most widely used tool for correcting technical biases, called batch effects, in microarray expression data.

## Minimum dependencies

We list here the versions of the paquages that have been used for development/testing of pyComBat, as well as for writing the documentation.

### pyCombat dependencies

* python 3.6

* numpy 1.16.4

* mpmath 1.1.0

* pandas 0.24.2

### Documentation

* sphinx 2.1.2

## Usage example

The simplest way of using pyComBat is to first import it, and simply use the pycombat function with default parameters:

```python
import pycombat
pycombat(data,batch)
```

* data: The expression matrix. It contains the information about the gene expression (rows) for each sample (columns). The first column (resp. row) is dedicated for the gene (resp. sample) names.

* batch: List of batch indexes. The batch list describes the batch for each sample. The list of batches contains as many elements as the number of columns in the expression matrix.

## How to contribute

Please refer to [CONTRIBUTING.md](https://github.com/epigenelabs/pyComBat/blob/master/CONTRIBUTING.md) to learn more about the contribution guidelines.

## References

(1) Johnson,W.E. et al. (2007) Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics, 8, 118–127