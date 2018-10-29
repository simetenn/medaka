# NEURON model for endocrine pituitary cells in Medaka

An implementation of the NEURON models created for the article:
"BK channels have opposite effects on sodium versus calcium mediated action potentials in endocrine pituitary cells.",
and all analysis of the model.


## Docker environment

We have created a [Docker](https://www.docker.com/) environment
with all dependencies installed.
This Docker environment can be started by running the bash script
`run_docker.sh` from within the this directory.
All results have been created in this Docker environment.


## Dependencies

The required dependencies are:

* `numpy`
* `matplotlib`
* `uncertainpy`
* `chaospy`

These can be installed with:

```
pip install numpy
pip install matplotlib
pip install uncertainpy
pip install chaospy
```

Additionaly the [Neuron](https://www.neuron.yale.edu/neuron/download) simulator
with the Python interface is required. NEURON must be manually installed
by the user.

## Running the code

To create Figure 3 and Figure 4 in
"BK channels have opposite effects on sodium versus calcium mediated action potentials in endocrine pituitary cells"
run:

```
python analysis.py
```

This takes around 7 minutes on a workstation computer.

To perform the uncertainty quantification and sensitivity analysis of the model
run:

```
python uq.py
```

This takes around 4 minutes on a workstation computer.


## Content

The content of this folder is:

* `medaka.py` - contains the Medaka model implementation as the functions `medaka_1` for MEDAKA 1 and `medaka_2` for MEDAKA 2.
* `burstines.py` - contains the functions for calculating the burstiness.
* `analysis.py` - contains the analysis of the Medaka model.
* `uq.py` - contains the uncertainty quantification and sensitivity analysis of the Medaka model.


## Platform and package specifications

All results have been generated inside a Docker environment with:

```
Platform: linux
Python: 3.7.0 (default, Jun 28 2018, 13:15:42)
[GCC 7.2.0]
Machine and architecture x86_64 64bit
NumPy: 1.15.1
matplotlib: 3.0.0
Chaospy: 2.3.5
Uncertainpy: 1.1.4
NEURON: 7.6.2-3-g9f36b13
```