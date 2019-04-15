[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1491552.svg)](https://doi.org/10.5281/zenodo.1491552)


# NEURON model for endocrine pituitary cells in medaka

An implementation of the NEURON models created for the article:
"BK channels have opposite effects on sodium versus calcium mediated action potentials in endocrine pituitary cells.",
and all analysis of the models.



## Content

The content of this folder is:

* `medaka.py` - contains the Medaka model.
* `burstines.py` - contains the functions for calculating the burstiness.
* `uq.py` - contains the uncertainty quantification and sensitivity analysis of the models.
* `*.mod` - NEURON files that implements the various ion channels.
* `platform_information.py` - prints platform information.


## Docker environment

We have created a [Docker](https://www.docker.com/) environment
with all dependencies installed.
This Docker environment can be started by running the bash script
`run_docker.sh` from within this directory.
All results have been created in this Docker environment.


## Dependencies

The required dependencies are:

* `numpy`
* `matplotlib`
* `uncertainpy`
* `chaospy`
* `NEURON`

These can be installed with:

```
pip install numpy
pip install matplotlib
pip install uncertainpy
pip install chaospy
```

Additionally, the [Neuron](https://www.neuron.yale.edu/neuron/download) simulator
with the Python interface is required. NEURON must be manually installed
by the user.

## Running the code

To perform the uncertainty quantification and sensitivity analysis of the model
run:

```
python uq.py
```


## Platform and package specifications

The uncertainty quantification and sensitivity analysis results have been generated inside a Docker environment with:

```
Platform: linux
Python: 3.7.1 (default, Dec 14 2018, 19:28:38)
[GCC 7.3.0]
Machine and architecture x86_64 64bit
NumPy: 1.15.2
matplotlib: 3.0.0
Chaospy: 2.3.5
Uncertainpy: 1.1.4
NEURON: 7.6.6
```
