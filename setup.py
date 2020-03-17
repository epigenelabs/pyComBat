import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="kaderbehdenna",
    version="0.1.0",
    author="Abdelkader Behdenna",
    author_email="abdelkader@epigenelabs.com",
    description="pyComBat, a Python tool for batch effects correction in high-throughput molecular data using empirical Bayes methods",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/epigenelabs/pyComBat",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0 License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)