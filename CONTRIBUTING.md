# Contributing to pyComBat

pyComBat has been thought from the beginning as a project that would be made available to the community as soon as possible. Also we already produced some nice results, the project is still in development, there is still room for improvement, and any contribution is important!

## How you can contribute

You have **many** ways to contribute to the pyComBat project, and make it successful by:

* Adding more useful features

* Detecting and fixing issues

* generating new interesting examples (and make the data available)

## Contributing in practice

### Adding a feature

If you want to implement a new interesting feature, please start by posting about your idea as an issue with `enhancement` tag. You are encourage to discuss it with the community and experts! Once all the necessary information is gathered, and that the idea is good (accepted by the community or urgently important) and feasible, you can start to contribute through pull requests!

### Signaling (and/or solving) an issue

Please start by searching amongst the issues list if the bug has already been signaled, and post about it if not.
Vigilante is very much appreciated, help us to spot if any issues is a duplicate. If you have any question, don't hesitate to contact us.

*PLEASE USE PULL REQUEST IF YOU WANT TO CONTRIBUTE.*

### Testing pyComBat on new datasets

We believe that the best way of increasing pyComBat robustness is to use it on as many various datasets as possible.

If you have data that can be of use in terms of testing pyComBat (and maybe discovering new possible bufs), please provide any benchmarking results you can produce, as well as the corresponding data, if possible.

In the case where you find any new bug, please refer to the previous section of this document.

Again, all actions should go through the pull request process!

## More ways of contributing

### Unit testing

A unit testing framework (under development) is ready in the form of the script `test_unit.py`. Simply run it by using `pytest` directly in the terminal.

You can install pytest by running the followinfg command in your command line:

```python
pip install -U pytest
```

Ideas to enhance and complete unit testing for pyComBat are welcome. Please follow the same procedure as for classic contribution over the main code.

### Writing documentation

pyComBat uses [Sphinx](https://www.sphinx-doc.org/en/master/usage/quickstart.html) to automatically document the main scripts. You can use the following command to build the documentation:

```bash
sphinx-build -b html source/ .
```
