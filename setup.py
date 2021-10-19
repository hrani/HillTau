""" setup.py : Script for FindSim """
__author__      = "HarshaRani"
__copyright__   = "Copyright 2018 HillTau, NCBS"
__maintainer__  = "HarshaRani"
__email__       = "hrani@ncbs.res.in"

import os
import sys
sys.path.append(os.path.dirname(__file__))
import setuptools
from distutils.core import setup, Extension

with open("README.md","r") as fh:
        description = fh.read()

setuptools.setup(
        name="HillTau",
        description="A fast, compact abstraction for model reduction in biochemical signaling networks",
        version="1.0",
        long_description = description,
        packages=["HillTau","HillTau.CppCode"],
        package_dir={'HillTau': "."},
        ext_modules=[Extension('ht', ['CppCode/htbind.cpp','CppCode/ht.cpp'], include_dirs=['extern/pybind11/include','extern/exprtk'])],
        install_requires = ['numpy','matplotlib','pydot','simplesbml','graphviz'],
        
        #dependency_links=[
        # Make sure to include the `#egg` portion so the `install_requires` recognizes the package
        #'git+ssh://git@github.com/pybind/pybind11.git#egg=ExampleRepo-0.1'],
	url ="http://github.com/Bhallalab/HillTau",
	package_data = {"HillTau" : ['HT_MODELS/*.json','*.xml','CppCode/hillTau.py']},
	license="GPLv3",
	entry_points = {
	
		'console_scripts' : [ 'HillTau = HillTau.__main__:run',
		                      'htgraph = HillTau.__main__:run_htgraph',
				       'ht2sbml = HillTau.__main__:run_ht2sbml'
				   ]
			},
		)

