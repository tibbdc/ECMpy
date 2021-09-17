# eciML1515
The process for enzyme-constrained model construction.

## About

The pipeline was written and tested with Python 3.8. The core libraries essential for the pipeline including: cobra, plotly(draw figures), and related packages. 

## Installation

1. create ECMpy environment using conda:

```shell
$ conda create -n ECMpy python=3.8
```

2. install related packages using pip:

```shell 
$ conda activate ECMpy
$ pip install cobra
$ pip install plotly
$ pip install -U kaleido
$ pip install nbformat
$ pip install ipykernel
$ python -m ipykernel install --user --name ECMpy --display-name "ECMpy"
```

## Steps to reproduce the analysis in the publication

Download all data and analysis code from github (directlt download or use git clone). 

 ```shell
$ cd /file path/project save path/
$ git clone https://github.com/tibbdc/ECMpy.git
```

 All results can be reproduced by executing the Jupyter Python notebooks:

+ 01.iML1515_modification_workflow.ipynb
  + get new iML1515 model corrected GPR relationships, reaction direction and EC number.

+ 02.construct_raw_eciML1515.ipynb
  + get raw eciML1515 using machine learning data.

+ 03.construct_final_eciML1515.ipynb
  + modify raw eciML1515 using enzyme usage and C13 data.
  
+ 04.simulation.ipynb
  + simulation results involved in the article.
  
