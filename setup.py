# -*- coding: utf-8 -*-                                                                                
                                                                                                       
###                                                                                                    
# © 2018 The Board of Trustees of the Leland Stanford Junior University                                
# Nathaniel Watson                                                                                     
# nathankw@stanford.edu                                                                                
###

# For some useful documentation, see
# https://docs.python.org/2/distutils/setupscript.html.
# This page is useful for dependencies: 
# http://python-packaging.readthedocs.io/en/latest/dependencies.html.

# PSF tutorial for packaging up projects:
# https://packaging.python.org/tutorials/packaging-projects

import glob
import os
from setuptools import setup, find_packages

SCRIPTS_DIR = "pulsarpy_to_dcc/scripts/"
scripts = glob.glob(os.path.join(SCRIPTS_DIR,"*.py"))
scripts.remove(os.path.join(SCRIPTS_DIR,"__init__.py"))

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
  author = "Nathaniel Watson",
  author_email = "nathankw@stanford.edu",
  classifiers = [
      "Programming Language :: Python :: 3",
      "License :: OSI Approved :: MIT License",
      "Operating System :: OS Independent",
  ],
  description = "Client for submitting experiments in Pulsar LIMS to the ENCODE Portal",
  install_requires = [
    "encode-utils",
    "pulsarpy",
    "inflection",
    "requests",
  ],
  long_description = long_description,
  long_description_content_type = "text/markdown",
  name = "pulsarpy-to-encodedcc",
  packages = find_packages(),
  project_urls = {
      "Read the Docs": "https://pulsarpy-to-encodedcc.readthedocs.io/en/latest",
  },
  scripts = scripts,
  url = "https://github.com/nathankw/pulsarpy-to-encodedcc", # home page
  version = "0.1.0",
)
