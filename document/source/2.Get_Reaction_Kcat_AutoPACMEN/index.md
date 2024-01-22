# 2.Get Reaction Kcat Using AutoPACMEN

```python
import cobra
import datetime 
import sys
sys.path.append(r'./script/')
from AutoPACMEN_function import *
from ECMpy_function import *
```



```python
# input files
autopacmen_folder = "./analysis/get_kcat_mw_by_AutoPACMEN/"
kcat_gap_fill= 'mean'#'mean'#'median'
reaction_gap_fill='mean'
sbml_path = "./data/iML1515_new.xml"
organism = "Escherichia coli"
project_name = "iML1515_%s"%kcat_gap_fill
create_file(autopacmen_folder)
protein_kcat_database_path = "none"
bigg_metabolites_file = "./data/bigg_models_metabolites.txt"#date:20230629 http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt
brenda_textfile_path = "./data/brenda_2023_1.txt"#date:20230427 https://www.brenda-enzymes.org/brenda_download/file_download.php
uniprot_data_file='./data/uniprot_data_accession_key.json'#from uniprot

#output files
brenda_json_path = "%skcat_database_brenda.json"%autopacmen_folder
brenda_json_path2 = "%ssa_database_brenda.json"%autopacmen_folder
sabio_rk_json_path = "%skcat_database_sabio_rk.json"%autopacmen_folder
bigg_id_name_mapping_path = "%sbigg_id_name_mapping.json"%autopacmen_folder
brenda_output_json_path = "%skcat_database_brenda_for_model.json"%autopacmen_folder
combined_output_path = "%skcat_database_combined.json"%autopacmen_folder
sub_description_path = '%sget_gene_subunitDescription.csv'%autopacmen_folder
gene_subnum_path = "%sgene_subnum.csv"%autopacmen_folder
reaction_mw_path = "%sreaction_mw.json"%autopacmen_folder
reaction_kcat_mw_path = '%sreaction_kcat_MW.csv'%autopacmen_folder
```

    Path exists
    


## Step 1: get bigg metbolite


```python
starttime=datetime.datetime.now()
# Step 1: get bigg metbolite
print("Starting to deal BIGG metabolites text file...")
parse_bigg_metabolites_file(bigg_metabolites_file, autopacmen_folder)
print("BIGG metabolites text file done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to deal BIGG metabolites text file...
    BIGG metabolites text file done!
    
    0:00:00.073214
    

## Step 2: BRENDA kcat


```python
starttime=datetime.datetime.now()
# Step 2: BRENDA kcat
print("Starting to deal BRENDA textfile...")
parse_brenda_textfile(brenda_textfile_path, autopacmen_folder, brenda_json_path, brenda_json_path2) 
print("BRENDA textfile done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to deal BRENDA textfile...
    WARNING: BRENDA text file line 3.6.1.34 (transferred to  EC 7.1.2.2, H+-transporting two-sector ATPase) is not interpretable!
    WARNING: BRENDA text file line 3.6.3.17 (transferred to subgroup EC 7.5.2.) is not interpretable!
    BRENDA textfile done!
    
    0:00:05.697793
    

## Step 3: Select Brenda kcat for model


```python
starttime=datetime.datetime.now()
# Step 3: Select Brenda kcat for model
print("Starting to deal brenda json for model...")
parse_brenda_json_for_model(sbml_path, brenda_json_path, brenda_output_json_path)
print("BRENDA json for model done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to deal brenda json for model...
    BRENDA json for model done!
    
    0:00:04.686739
    

## Step 4: SABIO-RK kcat for model


```python
starttime=datetime.datetime.now()
# Step 4: SABIO-RK kcat for model
print("Starting EC numbers kcat search in SABIO-RK...")
parse_sabio_rk_for_model_with_sbml(sbml_path, sabio_rk_json_path, bigg_id_name_mapping_path)
print("SABIO-RK done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting EC numbers kcat search in SABIO-RK...
    SABIO-RK done!
    
    0:20:05.752679
    

## Step 5: Brenda and SABIO-RK kcat combined


```python
starttime=datetime.datetime.now()
# Step 5: Brenda and SABIO-RK kcat combined
print("Combining kcat database...")
create_combined_kcat_database(sabio_rk_json_path, brenda_output_json_path, combined_output_path)
print("Combining kcat database done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Combining kcat database...
    Combining kcat database done!
    
    0:00:00.470836
    

## Step 6: subunit number of each reaction


```python
starttime=datetime.datetime.now()
# Step 6: subunit number of each reaction
print("Starting to fetch subunit number of each enzyme")
if re.search('\.xml',sbml_path):
    model = cobra.io.read_sbml_model(sbml_path)
elif re.search('\.json',sbml_path):
    model = cobra.io.json.load_json_model(sbml_path)
get_gene_subunitDescription(sub_description_path,model)#从uniprot的api下载，运行一次就行
subbnumdf = get_subunit_number(sub_description_path,gene_subnum_path)
print("Calculation done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to fetch subunit number of each enzyme
    Calculation done!
    
    0:00:00.132096
    

## Step 7: get mw for model gene (must be uniprot ID)


```python
starttime=datetime.datetime.now()
# Step 7: get mw for model gene (must be uniprot ID)
print("Starting UniProt ID<->Protein mass search using UniProt...")
get_protein_mass_mapping_from_local(sbml_path, autopacmen_folder, project_name, uniprot_data_file)
get_reaction_mw(sbml_path,autopacmen_folder, project_name, reaction_mw_path, gene_subnum_path)
print("Protein ID<->Mass mapping done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting UniProt ID<->Protein mass search using UniProt...
    Protein ID<->Mass mapping done!
    
    0:00:21.431180
    

## Step 8: kcat assignment for model(include sa)


```python
starttime=datetime.datetime.now()
# Step 8: kcat assignment for model(include sa)
print("Starting to assign kcat for model...")
get_reactions_kcat_mapping(sbml_path, autopacmen_folder, project_name, organism, combined_output_path,brenda_json_path2, reaction_mw_path,protein_kcat_database_path,kcat_gap_fill)
print("kcat assignment done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to assign kcat for model...
    kcat assignment done!
    
    0:00:09.861095
    

## Step 9: get_reaction_kcat_mw for model


```python
starttime=datetime.datetime.now()
# Step 9: get_reaction_kcat_mw for model
print("Starting to get reaction kcat_mw for model...")
if re.search('\.xml',sbml_path):
    model = cobra.io.read_sbml_model(sbml_path)
elif re.search('\.json',sbml_path):
    model = cobra.io.json.load_json_model(sbml_path)
get_reaction_kcat_mw(model,autopacmen_folder, project_name, reaction_gap_fill,gene_subnum_path,reaction_kcat_mw_path)       
print("Reaction kcat_mw done!")

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to get reaction kcat_mw for model...
    Default kcat is: 112645.64147536599
    Reaction kcat_mw done!
    0:00:15.806192
    


```python

```
