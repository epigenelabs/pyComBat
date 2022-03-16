import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    install_requires=[
        "numpy>=1.18.5,<=1.19.5",
        "mpmath>=1.1.0",
        "pandas>=0.24.2,<=1.1.5",
        "patsy==0.5.1"
    ],
    name="combat",
    version="0.3.1",
    author="Abdelkader Behdenna",
    author_email="abdelkader@epigenelabs.com",
    description="pyComBat, a Python tool for batch effects correction in high-throughput molecular data using empirical Bayes methods",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/epigenelabs/pyComBat",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
