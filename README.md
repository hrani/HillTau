![alt text](./Images/HillTau_Logo4_360px.png?raw=true "HillTau logo")

# HillTau
Sparse, efficient, non-ODE approximation for biochemical signaling

Copyright (C) Upinder S. Bhalla and NCBS-TIFR, 2020


## Contents of this repository:

README.md: this file

LICENSE: License terms. GPL 3.0.

PythonCode: The code to run HillTau.

CppCode: The code to run HillTau using python bindings via Pybind11 to C++.

MASH: Model Abstraction from SBML to HillTau. Utility program for model 
reduction, see 

[Documentation for MASH](Mashdoc.md)

ht2sbml: Utility program to convert HillTau (JSON format) model to a near-
	approximation in SBML, suitable for running on several ODE simulators
	like COPASI.

[Documentation for HillTau](Documentation.md)

[Preprint for HillTau](https://www.biorxiv.org/content/10.1101/2020.09.20.305250v1), which discusses many aspects of its design and use.

## Installation
Copy the two files hillTau.py and hillTauSchema.json from PythonCode to your
target directory. 

COMING SOON:

pip install hillTau<br>
This will help install the _much_ faster C++ version of HillTau, with the same
friendly Python interface.

## Versions
The Python version of HillTau has been tested with Python 2.7.17 and Python 3.6.9<br>
The C++ version of HillTau uses Python 3.6 or higher. It is the recommended
version.

## Examples
Examples: Directory with examples

Examples/HT_MODELS: Examples of HillTau models

Examples/KKIT_MODELS: Examples of KKIT models which map to the HillTau ones.
	KKIT models are an old ODE biochemical model format compatible with
	GENESIS and MOOSE.

Examples/SBML_MODELS: Examples of SBML models which map to the HillTau ones.
	SBML is the Systems Biology Markup Language and is a standard for defining
	many kinds of systems biology models, including the current mass-action
	ones.

Examples/PaperFigures: Using the HillTau form to generate the figures for the
	reference paper for HillTau. Most of these require MOOSE to be 
	installed, but fig1.py just requires HillTau.

Other projects and papers that relate to HillTau: [Resources.md](Resources.md)
