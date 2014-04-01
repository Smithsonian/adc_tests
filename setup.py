import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "adc5g",
    version = "0.0.1",
    author = "Rurik A Primiani",
    author_email = "rprimian@cfa.harvard.edu",
    description = "Test Utilities for the ASIAA 5 GSps CASPER ADC",
    long_description=read('README.md'),
    packages=['adc5g',],
)
