![alt text](./Images/HillTau_Logo4_360px.png?raw=true "HillTau logo")
# Documentation for HillTau

HillTau has two major usage modes: as a standalone simulator, and as a
library from which HillTau calls are used by other programs. Both these 
are supported by the same two HillTau files: 

	hillTau.py and hillTauSchema.json


## Standalone HillTau syntax

In the present document we provide some examples of how to use HillTau in
standalone mode.

### Getting the syntax: 

```
python hillTau.py -h

	usage: hillTau.py [-h] [-r RUNTIME] [-s STIMULUS [STIMULUS ...]] [-p PLOTS] model

	This is the hillTau simulator. This program simulates abstract kinetic/neural
	models defined in the HillTau formalism. HillTau is an event-driven JSON form to
	represent dynamics of mass-action chemistry and neuronal activity in a fast,
	reduced form. The hillTau program loads and checks HillTau models, and
	optionally does simple stimulus specification and plotting
	 
	positional arguments:
  	model                 Required: filename of model, in JSON format.
	
	optional arguments:
  	-h, --help            show this help message and exit
  	-r RUNTIME, --runtime RUNTIME
	                       Optional: Run time for model, in seconds. If flag is
	                       not set the model is not run and there is no display
	-s STIMULUS [STIMULUS ...], --stimulus STIMULUS [STIMULUS ...]
	                       Optional: Deliver stimulus as follows: --stimulus
	                       molecule conc [start [stop]]. Any number of stimuli
	                       may be given, each indicated by --stimulus. By
	                       default: start = 0, stop = runtime
	 -p PLOTS, --plots PLOTS
	                       Optional: plot just the specified molecule(s). The
	                       names are specified by a comma-separated list.
```


### Running a HillTau model
	
Go into the Examples directory and type:

	python ../PythonCode/hillTau.py HT_MODELS/osc.json -r 5000

![alt text](./Images/osc.png?raw=true "Oscillatory model")

This runs the oscillatory model with runtime 5000 seconds, and plots it.

### Selecting only a subset of plots

	python ../PythonCode/hillTau.py HT_MODELS/osc.json -r 5000 -p output
![alt text](./Images/osc_output.png?raw=true "Display only one plot")

	python ../PythonCode/hillTau.py HT_MODELS/osc.json -r 5000 -p output,nfb
![alt text](./Images/osc_output_nfb.png?raw=true "Display two selected plots")

### Giving a stimulus

	python ../PythonCode/hillTau.py HT_MODELS/exc.json -r 20 -s input 1e-3 5 10
![alt text](./Images/singleStim_exc.png?raw=true "Apply single stimulus")

Here we assign molecule 'input' to value 1 uM at time 5 seconds, and
the stimulus turns off at 10 seconds. This is shown in the blue plot.
Note that the units of the model are the default millimolar, so the
stimulus must also be in the same units.

### Giving multiple stimuli

	python ../PythonCode/hillTau.py HT_MODELS/bcm_bistable.json -r 100 -s Ca 2 10 11 -s Ca 0.3 50 80 -p Ca,synAMPAR,on_CaMKII
![alt text](./Images/doubleStim_synapse.png?raw=true "Apply double stimulus")

This model has a bistable mediated by CaMKII feedback, which is triggered by
a Ca input. There is also a Ca-activated phosphatase CaN which acts to 
inhibit and hence turn off the activity of the CaMKII. Thus the Ca input 
both excites the CaMKII, and via CaN, inhibits it. These two counteracting
processes act at different speeds and concentrations.

Here we give a brief strong Ca stimulus to turn on the bistable, and a long,
low Ca stimulus to turn it off again. The CaMKII in turn activates the
synaptic AMPA receptor (synAMPAR) as a readout of synaptic weight. Note that 
the model units are in uM (micromolar), and so the stimulus units are 
handled also in uM.


## Use of HillTau as a library

HillTau is a Python script, so there are a lot of functions you can access
if you want. Only a few are needed or recommended:

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
	Model object. The Model object is your handle for running the model.

	Argument: the dictionary of the model as loaded from JSON.

	Returns: HillTau Model object
	
Once you have your model, you can run HillTau simulations.

1. 	model.reinit()

	Reinitializes the simulation time to zero, reinitializes all the
	state variables to their starting values.

2.	model.advance( advanceTime, settle = False )

	This advances the simulation by the specified time. The optional 
	_settle_ flag, when True, tells HillTau that intermediate 
	time-points are not needed and to jump very quickly to the steady-state.

### Frequently used classes

There are a couple of frequently used classes.

**MolInfo**: 

This class contains information about each molecule species. Relevant
fields are name and index.

**Model**:

The Model class exposes a few additional fields and functions.

1.	model.molInfo. This is a dict of MolInfos. It is indexed by the name
	of the molecule. The most common use is of the form
	
		```myIndex = input.molInfo.get( "MyMoleculeName" ).index```

2.	model.conc. This is a numpy array of molecule concentrations, indexed
	as using the molecule index. You can get or set it.

		```
		myConc = model.conc[myIndex]
		model.conc[myIndex] = myConc * 2
		```


3.	model.concInit. This is an array of molecule initial concentrations, 
	indexed by the molecule index as above. At 
	```model.reinit()```
	the model.conc vector is initialized to model.concInit.

4.	model.plotvec. This is a list of time-series values of all the
	molecules in the simulation. Every time-step, the entire model.conc
	array is appended to the plotvec. This is how you would get the 
	vector of
	values for myMolecule:

	```myVec = np.transpose( np.array( model.plotvec ) )[myIndex]```

5.	model.dt: This is the timestep of the simulation. User can set it.

6.	model.currentTime: Current time of simulation. User must not set it.




## HillTau model specification format

HillTau uses a JSON format and this is fully specified by a schema file:

	hillTauSchema.json

It would in principle be possible to write HillTau models using SBML, though
most SBML models are ODE models. In other words, regular SBML simulators would
not be able to compute HillTau models, nor would HillTau know what to do with
most existing SBML models. Should there be sufficient interest in this
I could take it up.

### Units

HillTau uses seconds for time units.

The **quantityUnits** property specifies one of ["M", "mM", "uM", "nM", "pM"]
as an optional unit for concentration. The default is mM.

### Groups

HillTau organizes all reaction sets into groups. There can be as many groups
as the reaction system requires. This is purely an organizational feature and
has no computational implications, though for systems of any complexity it 
is essential to group reactions in order to keep track and to map properly to
known signaling pathways.

In each group there can be further JSON objects, namely

-	Species

	This is a list of name:value pairs, to specify species name and its 
	initial concentration.
-	Reacs

	Each of these is an object with:
	-	Name

		The Reaction name automatically becomes the name of the product
		of the reaction. This product is like any other molecule and can
		be used as a substrate in other reactions or equations.
	-	subs[sub1, sub2...]

		Required list of substrates (array of names)
	-	Reaction required parameters *KA* and *tau*, both floats.
	-	Reaction optional parameters *tau2* and *baseline*, both floats.
	-	Optional flag for *inhibit*: 0 or 1. Default is 0.
-	Eqns
	-	Name

		The Equation name automatically becomes the name of a molecule
		whose value is defined by the equation. This can be used as
		a substrate in other reactions or equations.

	-	Equation string

		This is a string expressing an algebraic function to evaluate.
		The function can use any named molecule, standard python 
		operations, and numbers.

### Species names and namespaces

HillTau creates a molecular species for each of the following:
-	Every species defined with the *Species* keyword
-	Every reaction name
-	Every equation name

HillTau namespace is flat and global, that is, any species defined anywhere
in the system is accessible anywhere in the system.

One can use the *Species* list to initialize the values for molecules 
subsequently defined as a *Reac* or as an *Eqn*. If this initialization
is not done, then initial values the *Reac* and *Eqn* molecules are computed
from the values of the input molecules at initialization time.


## HillTau calculations

HillTau calculations are best explained in the code and in the paper,
mentioned in the [Resources.md](Resources.md) file

Briefly, at each timestep the system calculates the steady-state value for each
reaction using a Hill function, and then uses an exponential decay calculation
to find how far the system would approach it.

The output value of any reaction depends only on its inputs. It is not affected
by any number of downstream reactions that it may plug into. This differs
fundamentally from chemical reactions. It greatly simplifies design of
models and analysis of signal flow, because all information flow is forward.

All calculations are done using simple Python functions. Despite this, the
basic HillTau algorithm is so efficient that we have not yet felt the need
to build a matching C++ library using pybind 11. This is on the cards and
if anyone feels that this is essential they should let us know.

### HillTau outputs

Output values for all molecules in a HillTau model are stored in the

	model.plotvec

field, which is a list of numpy arrays holding all molecule concs at each 
timestep. Please see *Use of HillTau as a library*, above. One can also sample
from the 

	model.conc

array, at each timestep. This is also documented above.





