""" setup.py : Script for FindSim """
__author__      = "HarshaRani"
__copyright__   = "Copyright 2021 HillTau, NCBS"
__maintainer__  = "HarshaRani"
__email__       = "hrani@ncbs.res.in"

import os
import sys
sys.path.append(os.path.dirname(__file__))
import setuptools
#from distutils.core import setup, Extension
from setuptools import setup, Extension
setuptools.setup(
        name="HillTau",
        description="A fast, compact abstraction for model reduction in biochemical signaling networks",
        version="1.0",
        author= "Upinder S. Bhalla",
        author_email="bhalla@ncbs.res.in",
        maintainer= "HarshaRani",
        maintainer_email= "hrani@ncbs.res.in",
        long_description = open('README.md', encoding='utf-8').read(),
        long_description_content_type='text/markdown',
        packages=["HillTau","HillTau.CppCode"],
        package_dir={'HillTau': "."},
        ext_modules=[Extension('ht', ['CppCode/htbind.cpp','CppCode/ht.cpp'], include_dirs=['extern/pybind11/include','extern/exprtk'])],
        install_requires = ['numpy','matplotlib','graphviz','pydot-ng','pydot','simplesbml'],
        
        #dependency_links=[
        # Make sure to include the `#egg` portion so the `install_requires` recognizes the package
        #'git+ssh://git@github.com/pybind/pybind11.git#egg=ExampleRepo-0.1'],
	url ="http://github.com/Bhallalab/HillTau",
	package_data = {"HillTau" : ['HT_MODELS/*.json','*.xml','CppCode/hillTau.py']},
	license="GPLv3",
	entry_points = {
	
		'console_scripts' : [ 'HillTau = HillTau.__main__:run',
				       'Hilltau = HillTau.__main__:run',
				       'hillTau = HillTau.__main__:run',
				       'hilltau = HillTau.__main__:run',
		                      'htgraph = HillTau.__main__:run_htgraph',
				       'ht2sbml = HillTau.__main__:run_ht2sbml'
				   ]
			},
		)

