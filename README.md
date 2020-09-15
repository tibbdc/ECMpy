# ECM_eco1515
The  automated process for enzyme-constrained model construction.

## About

The pipeline was written and tested with Python 3.5, though Python versions (Python 3.6) should work. The core libraries essential for the pipeline including: Cobra, Pandas, and related packages. 

## Software

The packages used to run the code in the pipeline was listed in requirements.txt. To install the requirements using pip, run the following at command-line:

```shell
$ pip install -r requirements.txt
```

To create a stand-alone environment named ECM with Python 3.5 and all the reqiured package versions, run the following:

```shell
$ conda create -n ECM python=3.5 
$ source activate ECM
$ pip install -r requirements.txt
```

You can read more about using conda environments in the [Managing Environments](http://conda.pydata.org/docs/using/envs.html) section of the conda documentation. 

## Steps to reproduce the analysis in the publication

 All results can be reproduced by executing the Jupyter Python notebooks:

+ Workflow_construction_of_ECM_eco1515.ipynb
  + the main script of construction enzyme-constrained model
+ Simulation.ipynb
  + the script for simulating part of results in manuscript

