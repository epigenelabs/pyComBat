.. pyComBat documentation master file, created by
   sphinx-quickstart on Wed Feb 19 17:44:47 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyComBat's documentation!
====================================

General Overview
----------------

Variability in datasets are not only the product of biological processes: they are also the product of technical biases (Lander et al, 1999).
ComBat is one of the most widely used tool for correcting those technical biases called batch effects.

pyComBat (Behdenna et al, 2020) is a new Python implementation of ComBat (Johnson et al, 2007), a software widely used for the adjustment of batch effects in microarray data.
While the mathematical framework is strictly the same, pyComBat:

    1. has similar results in terms of batch effects correction;

    2. is as fast or faster than the R implementation of ComBat and;

    3. offers new tools for the community to participate in its development.

Implementation
--------------
pyComBat is an open source program written in Python language version 3.7.3. It can only be run as a command line. It is available at https://github.com/epigenelabs/pyComBat.

License
-------

pyComBat is implemented in the Python language and is available under `GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.en.html>`_ license.

You can find more detailed information in the LICENSE file.

Installing pyComBat
-------------------

pyComBat is usable as a Python package, which is installed with:

.. code-block:: Python

    pip install combat

To upgrade the package to the latest version, please use:

.. code-block:: Python

    pip install combat --upgrade

Then, you can import the pycombat function:

.. code-block:: Python

    from combat.pycombat import pycombat

Using pyComBat
--------------

Minimal usage example:

.. code-block:: Python

    from combat.pycombat import pycombat
    import numpy as np
    import pandas as pd

    # prepare simulated data
    matrix = np.transpose([np.random.normal(size=1000,loc=3,scale=1),np.random.normal(size=1000,loc=3,scale=1),np.random.normal(size=1000,loc=3,scale=1),
                      np.random.normal(size=1000,loc=2,scale=0.6),np.random.normal(size=1000,loc=2,scale=0.6),
                      np.random.normal(size=1000,loc=4,scale=1),np.random.normal(size=1000,loc=4,scale=1),np.random.normal(size=1000,loc=4,scale=1),np.random.normal(size=1000,loc=4,scale=1)])

    df = pd.DataFrame(matrix,
                 columns=["sample"+str(i) for i in range(9)],
                 index=["gene"+str(i) for i in range(1000)])

    batch = np.asarray([1,1,1,2,2,3,3,3,3])

    # run pyComBat
    df_corrected = pycombat(df,batch)

    # visualise results
    import matplotlib.pyplot as plt
    plt.boxplot(df.transpose())
    plt.boxplot(df_corrected.transpose())



Documentation for the code
==========================
.. toctree::
   :maxdepth: 2
   :caption: Contents:


pyComBat utils
--------------
.. automodule:: pycombat
    :members: pycombat


Contributing to pyComBat
========================

Contribution guidelines
-----------------------

Contribution guidelines can be found in `CONTRIBUTING.md <https://github.com/epigenelabs/pyComBat/blob/master/CONTRIBUTING.md>`_.

Unit Testing
------------

Most of the subfunctions can be tested separately. 
The "unit_test" script implements all of them, and can be used to check the good functioning of the whole pyComBat software.


Authors
=======
Contact
-------
To  ask  a  question  on  pyComBat,  report  a  suggestion  (e.g.   why  not  including  other options) or if you think you have discovered a bug (if any?), please contact:

Abdelkader Behdenna at abdelkader@epigenelabs.com

Citing pyComBat
---------------
A. Behdenna, J. Haziza, C.-A. Azencott and A. Nordor. 2020. pyComBat, a Python tool for batch effects correction in high-throughput molecular data using empirical Bayes methods. bioRxiv doi: 10.1101/2020.03.17.995431


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`