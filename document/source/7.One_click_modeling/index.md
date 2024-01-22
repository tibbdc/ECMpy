```python
conda activate ECMpy2
cd userwokplace/ECMpy2.0
```
# 7. One Click Modeling
## Use AutoPACMEN data to reconstruct ecGEM


```python
python ./script/get_ecGEM_onestop.py -m "./data/iML1515_new.xml" -kcat "No" -f "No" -bigg  "./data/bigg_models_metabolites.txt" -gene_abundance "./data/gene_abundance.csv" -tax_id 83333 -org "Escherichia coli" -sigma 1 -ptot 0.56 -kcat_method "AutoPACMEN" -work_folder "./analysis/get_kcat_mw" -brenda "./data/brenda_2023_1.txt" -uniprot "./data/uniprot_data_accession_key.json" -kcat_gap_fill "mean" -r_gap_fill "mean" -ecGEM "./model/ecGEM_autopacmen.json"
```

## Use DLKcat data to reconstruct ecGEM


```python
python ./script/get_ecGEM_onestop.py -m "./data/iML1515_new.xml" -kcat "No" -f "No" -bigg  "./data/bigg_models_metabolites.txt" -gene_abundance "./data/gene_abundance.csv" -tax_id 83333 -sigma 1 -ptot 0.56 -kcat_method "DLKcat" -work_folder "./analysis/get_kcat_mw"  -ecGEM "./model/ecGEM_DLKcat.json"
```

## Use user's data to reconstruct ecGEM


```python
python ./script/get_ecGEM_onestop.py -m "./data/iML1515_new.xml" -kcat "./data/reaction_kcat_MW_template.csv" -f 0.405 -bigg  "./data/bigg_models_metabolites.txt"  -sigma 1 -ptot 0.56 -ecGEM "./model/ecGEM_userdata.json"
```


```python
python ./script/get_ecGEM_onestop.py -m "./data/iML1515_new.xml" -kcat "./data/reaction_kcat_MW_template.csv" -f "No" -bigg  "./data/bigg_models_metabolites.txt"  -gene_abundance "./data/gene_abundance.csv" -tax_id 83333 -sigma 1 -ptot 0.56 -ecGEM "./model/ecGEM_userdata_no_f.json"
```


```python

```
