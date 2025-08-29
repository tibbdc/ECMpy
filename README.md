# ECMpy2.0

Automated construction of enzyme-constrained models using ECMpy workflow.

## 1. Create environment

```shell
$ conda create -n ECMpy2 python=3.8  
$ conda activate ECMpy2
```

## 2. Install the relevant packages

#### Install package

```shell
$ pip install cobra openpyxl requests pebble xlsxwriter Bio Require quest scikit-learn  RDKit seaborn pubchempy torch bioservices==1.10.4 pyprobar xmltodict plotly kaleido nbformat jupyterlab ipykernel
```

## 3. Preprocessing data sources

The "all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration50","atom_dict.pickle", "bond_dict.pickle", "edge_dict.pickle", 'fingerprint_dict.pickle", and "sequence_dict.pickle" files are derived from the DLKcat method, and you can update it from GitHub(https://github.com/SysBioChalmers/DLKcat.git).
The 'bigg_models_metabolites.txt" file is downloaded from BiGG (http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt).
The "brenda_2023_1.txt" file is downloaded from BRENDA ([https://www.brenda-enzymes.org/download.php](https://www.brenda-enzymes.org/download.php)), and "EC_kcat_max.json" is obtained from this file extraction.
The "gene_abundance.csv" file is downloaded and transformed from PaxDB (https://pax-db.org/download).
The "uniprot_data_accession_key.json" is compiled from the UniProt database (only for Swiss-Prot), and we have uploaded to zenodo (https://zenodo.org/record/8119567/files/uniprot_data_accession_key.json?download=1).
The "AutoPACMEN_function.py" file is downloaded and modified from the AutoPACMEN method (https://github.com/klamt-lab/autopacmen.git).

## 4. Documentation

Full documentation is available at https://ecmpy.readthedocs.io/en/latest/.

### Detailed process for constructing enzyme-constrained Models.

+ 00.Model_preview.ipynb
  + Assessment of gene coverage (UniProt ID coverage), reaction coverage (EC number coverage excluding exchange reactions), and metabolite coverage (BiGG ID coverage).
+ 01.get_reactiion_kcat_using_DLKcat.ipynb
  + Using DLKcat for predicting enzyme kinetic parameters directly based on the sequence information of enzymes catalyzing reactions and substrate information.
+ 01.get_reaction_kcat_using_AutoPACMEN.ipynb
  + Employing the AutoPACMEN process for extracting enzyme kinetic parameter information from the BRENDA and SABIO-RK databases.
+ 02.get_ecModel_using_ECMpy.ipynb
  + Using the ECMpy process to construct ecGEM.
+ 03.ecModel_calibration.ipynb
  + An automated parameter calibration process for the ecModel, guided by the principle of enzyme utilization.
+ 04.ecModel_analysis.ipynb
  + Some analysis cases of ecModels.
+ 05.ecModel_ME.ipynbP
  + Predicting metabolic engineering targets using ecModels.
+ 06.One-click_modeling.ipynb
  + Constructing ecGEMs with a one-click approach through the command line.
+ 07.BiGG_to_ecGEM.ipynb
  + Constructing ecGEMs with a one-click approach through the command line for BiGG models.

## 5. Acknowledgement

Here we are deeply grateful to klamt-lab for releasing the code for AutoPACMEN (https://github.com/klamt-lab/autopacmen) and to SysBioChalmers for sharing the code for DLKcat (https://github.com/SysBioChalmers/DLKcat), which enables ECMpy2.0 to rapidly obtain enzyme kinetics parameter information for the corresponding models. We extend our heartfelt thanks to qLSLab for making the code for GPRuler available (https://github.com/qLSLab/GPRuler), as it has inspired ideas for ECMpy2.0 to automatically acquire the subunit composition of proteins.

## 6. How to cite:

Zhitao Mao, Xin Zhao, Xue Yang, Peiji Zhang, Jiawei Du, Qianqian Yuan and Hongwu Ma, ECMpy, a Simplified Workflow for Constructing Enzymatic Constrained Metabolic Network Model,Biomolecules, 2022; [https://doi.org/10.3390/biom12010065](https://doi.org/10.3390/biom12010065)

Zhitao Mao, Jinhui Niu, JianxiaoZhao, Yuanyuan Huang, Ke Wu, Liyuan Yun, Jirun Guan, Qianqian Yuan, Xiaoping Liao, Zhiwen Wang, Hongwu Ma, ECMpy 2.0: A Python package for automated construction and analysis of enzyme-constrained models,Synthetic and Systems Biotechnology, 2024; [https://doi.org/10.1016/j.synbio.2024.04.005](https://doi.org/10.1016/j.synbio.2024.04.005 "Persistent link using digital object identifier")
