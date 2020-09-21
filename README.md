![alt text](./Images/HillTau_Logo4_360px.png?raw=true "HillTau logo")

# HillTau
Sparse, efficient, non-ODE approximation for biochemical signaling

Copyright (C) Upinder S. Bhalla and NCBS-TIFR, 2020


## Contents of this repository:

README.md: this file

LICENSE: License terms. GPL 3.0.

PythonCode: The code to run HillTau.

[Documentation for HillTau](Documentation.md)

[Preprint for HillTau](https://www.biorxiv.org/content/10.1101/2020.09.20.305250v1), which discusses many aspects of its design and use.

## Installation
Copy the two files hillTau.py and hillTauSchema.json from PythonCode to your
target directory.

## Versions
HillTau has been tested with Python 2.7.17 and Python 3.6.9

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
