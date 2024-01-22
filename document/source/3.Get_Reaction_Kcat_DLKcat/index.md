# 3.Get Reaction Kcat Using DLKcat



```python
import cobra
import datetime 
import pandas as pd
import subprocess
# from script.ECMpy_function import *
import sys
sys.path.append(r'./script/')
from ECMpy_function import *
```



```python
dlkcat_folder = "./analysis/get_kcat_mw_by_DLkcat/"
create_file(dlkcat_folder)

# input files
sbml_path = "./data/iML1515_new.xml"

# output files
gene_subnum_path = "%sgene_subnum.csv"%dlkcat_folder
sub_description_path = '%sget_gene_subunitDescription.csv'%dlkcat_folder
inchikey_list_file='%sinchikey_list.csv'%dlkcat_folder
inchikey_list_smilesfile='%sinchikey_list_smiles.csv'%dlkcat_folder
comdf_file= '%scomdf.csv'%dlkcat_folder
DLouputdf_file = '%sDLoutput.tsv'%dlkcat_folder
metdf_outfile='%smetabolites_reactions_gpr_similes_prosequence_mass_dropna.csv'%dlkcat_folder
metabolites_reactions_gpr_file = '%smetabolites_reactions_gpr.csv'%dlkcat_folder
prodf_file = '%sprodf.csv'%dlkcat_folder
DLinput_file= '%sDLinput.tsv'%dlkcat_folder
DL_reaction_kact_mw_file='%sreaction_kcat_MW.csv'%dlkcat_folder
```

    Path exists
    


## Step 0: read GEM


```python
# Step 0: read GEM
if re.search('\.xml',sbml_path):
    model = cobra.io.read_sbml_model(sbml_path)
elif re.search('\.json',sbml_path):
    model = cobra.io.json.load_json_model(sbml_path)
```

## Step 1: subunit number of each reaction


```python
starttime=datetime.datetime.now()
# Step 1: subunit number of each reaction
print("Starting to fetch subunit number of each enzyme")
get_gene_subunitDescription(sub_description_path,model)#Download from the UniProt API, run it once.
subbnumdf = get_subunit_number(sub_description_path,gene_subnum_path)
print("Calculation done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to fetch subunit number of each enzyme
    Start downloading from UniProt...
    [0;38;2;66;227;35m100.00%[0;38;2;186;189;250m|█████████████████████████████|[0;38;2;237;166;178m 0:00:00|1:28:45 [0;38;2;146;52;247m ETC: 07-06 00:02:23[0m[K
    Success downloading! :-)
    Calculation done!
    
    1:28:49.808229
    

## Step 2: convert metbolites bigg id to smiles 


```python
starttime=datetime.datetime.now()
# Step 2: convert metbolites bigg id to smiles 
print("Starting to convert metbolites bigg id to smiles...")
metdf_name = get_met_bigg_id(model)
inchkeydf = convert_bigg_met_to_inchikey(metdf_name['met'],inchikey_list_file)#from BIGG
# inchkeydf = pd.read_csv('./data/inchikey_list.csv')
smilesdf = convert_inchikey_to_smiles(inchkeydf,inchikey_list_smilesfile)#from pubchem
print("Converting done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to convert metbolites bigg id to smiles...
    Converting...
    [0;38;2;66;227;35m100.00%[0;38;2;190;231;233m|█████████████████████████████|[0;38;2;28;97;15m 0:00:00|1:09:21 [0;38;2;146;52;247m ETC: 07-06 10:52:44[0m[K
    [] try again later
    Fail secure!
    Converting done!
    
    2:01:59.345124
    

## Step 3: get protein sequence and mass in model 


```python
starttime=datetime.datetime.now()
# Step 3: get protein sequence and mass in model 
print("Starting to get protein sequence and mass in model...")
subbnumdf = pd.read_csv(gene_subnum_path)
prodf = get_model_protein_sequence_and_mass(model,subbnumdf,prodf_file)
print("Getting done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to get protein sequence and mass in model...
    Getting done!
    
    0:22:14.932906
    

## Step 4: split the substrate of reactions to match the gene


```python
starttime=datetime.datetime.now()
# Step 4: split the substrate of reactions to match the gene
print("Starting to split the substrate of reactions to match the gene...")
spdf = split_substrate_to_match_gene(model,metabolites_reactions_gpr_file)
print("Splitting done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to split the substrate of reactions to match the gene...
    Splitting done!
    
    0:00:10.310155
    

## Step 5: combine the reaction--substrate--gene--protein_sequnce--mass and formate DLKcat input file


```python
starttime=datetime.datetime.now()
# Step 5: combine the reaction--substrate--gene--protein_sequnce--mass and formate DLKcat input file
print("Starting to combine data...")
metdf_name = get_met_bigg_id(model)
smilesdf = pd.read_csv(inchikey_list_smilesfile)
spdf = pd.read_csv(metabolites_reactions_gpr_file)
prodf = pd.read_csv(prodf_file)
comdf = combine_reactions_simles_sequence(spdf,smilesdf,prodf,comdf_file)
DLinputdf = generate_DLKCAT_input(comdf,metdf_name,metdf_outfile,DLinput_file)
print("Combinning done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to combine data...
    DLKCAT input file generated
    Combinning done!
    
    0:00:03.239874
    

## Step 6: use DLKcat calculate kcat


```python
starttime=datetime.datetime.now()
# Step 6: use DLKcat calculate kcat
print("Starting to Use DLKcat calculate kcat...")
cmd_str = "python ./script/prediction_for_input.py ./analysis/get_kcat_mw_by_DLkcat/DLinput.tsv ./analysis/get_kcat_mw_by_DLkcat/DLoutput.tsv"
subprocess.run(cmd_str, shell=True)
print("DLKcat done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to Use DLKcat calculate kcat...
    ./analysis/get_kcat_mw_by_DLkcat/DLinput.tsv
    It's time to start the prediction!
    -----------------------------------

    It takes 1198.0919036865234 seconds to predict Kcat values!
    -----------------------------------
    Prediction success!
    DLKcat done!
    
    0:20:00.334137
    

## Step 7: get the kcat_mw file


```python
starttime=datetime.datetime.now()
# Step 7: get the kcat_mw file
print("Starting to get reaction kcat_mw for model......")
DLouputdf = pd.read_csv(DLouputdf_file, sep='\t')
comdf = pd.read_csv(comdf_file)
DL_reaction_kact_mw = DL_kcat_mw_calculation(DLouputdf, comdf)
DL_reaction_kact_mw.to_csv(DL_reaction_kact_mw_file, index=False)
print("Reaction kcat_mw done!")

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to get reaction kcat_mw for model......
    DL_reaction_kact_mw generated
    Reaction kcat_mw done!
    0:00:00.084022
    
