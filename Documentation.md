# Documentation for HillTau

HillTau has two major usage modes: as a standalone simulator, and as a
library from which HillTau calls are used by other programs. Both these 
are supported by the same two HillTau files: 
hillTau.py and hillTauSchema.json

## Standalone Simulator examples

- Getting the syntax: 

	python hillTau.py -h

	> usage: hillTau.py [-h] [-r RUNTIME] [-s STIMULUS [STIMULUS ...]] [-p PLOTS] model

	> hillTau simulator This program simulates abstract kinetic/neural models
	> defined in the HillTau formalism. HillTau is an event-driven JSON form to
	> represent dynamics of mass-action chemistry and neuronal activity in a fast,
	> reduced form. The hillTau program loads and checks HillTau models, and
	> optionally does rudimentary stimulus specification and plotting
	> 
	> positional arguments:
  	> model                 Required: filename of model, in JSON format.
	> 
	> optional arguments:
  	> -h, --help            show this help message and exit
  	> -r RUNTIME, --runtime RUNTIME
	>                        Optional: Run time for model, in seconds. If flag is
	>                        not set the model is not run and there is no display
	>  -s STIMULUS [STIMULUS ...], --stimulus STIMULUS [STIMULUS ...]
	>                        Optional: Deliver stimulus as follows: --stimulus
	>                        molecule conc [start [stop]]. Any number of stimuli
	>                        may be given, each indicated by --stimulus. By
	>                        default: start = 0, stop = runtime
	>  -p PLOTS, --plots PLOTS
	>                        Optional: plot just the specified molecule(s). The
	>                        names are specified by a comma-separated list.


## Use of HillTau as a library

HillTau is a Python script, so there are a lot of functions you can access
if you want. Only four are needed or recommended:

1. loadHillTau( filename )
	Returns: Dictionary of HillTau model as loaded from JSON.
2. getQuantityScale( jsonDict )
	Argument: the dictionary of the model as loaded from JSON.
	Returns: scale factor to use for all quantity conversions. Reference
	units are mM.
3. scaleDict( jsonDict, qs )
	This function goes through the loaded JSON dictionary and rescales all
	concentration quantities by the provided quantityScale qs.
	Argument 1: the dictionary of the model as loaded from JSON.
	Argument 2: the quantity scale.
	There is no return value.
4. parseModel( jsonDict )
	This function parses the JSON dictionary and converts it to a HillTau
	Model object. The Model object is your handle for running the model
	Argument: the dictionary of the model as loaded from JSON.
	Returns: HillTau Model object
	
Once you have your model, you can run HillTau simulations.

1.	model.advance( advanceTime, settle = False )
	This advances the simulation by the specified time. The optional 
	_settle_ flag, when True, tells HillTau that intermediate 
	time-points are not needed and to jump very quickly to the steady-state.
2. 	model.reinit()
	Reinitializes the simulation time to zero, reinitializes all the
	state variables to their starting values.
