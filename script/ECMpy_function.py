# -*- coding: utf-8 -*-
# This code is used to introduce enzyme concentration constraint in GEMs
# by COBRApy and to calculate the parameters that need to be entered
# during the construction of the enzyme-constrained model.
#from warnings import warn

import json
import math
import random
import re
import os
import time
import statistics
from typing import Any, Dict, List
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import cobra
import numpy as np
import pandas as pd
from cobra.core import Reaction
from cobra.io.dict import model_to_dict
from cobra.util.solver import set_objective
import pubchempy as pcp
import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import xmltodict
from bioservices import UniProt
from pyprobar import bar, probar
import plotly
from copy import deepcopy
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from optlang.symbolics import Zero, add
#from tqdm import tqdm
import subprocess
from AutoPACMEN_function import *

def create_file(store_path):
    """
    Create a directory at the specified path.

    Args:
        store_path (str): The path of the directory to create.

    """
    if os.path.exists(store_path):
        print("Path exists")
        # Perform any necessary actions if the path already exists
        # For example, remove the existing directory and create a new one
        # shutil.rmtree(store_path)
        # os.makedirs(store_path)
    else:
        os.makedirs(store_path)
        print(store_path)
        
def json_load(path: str) -> Dict[Any, Any]:
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: str ~ The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary

def standardize_folder(folder: str) -> str:
    """Returns for the given folder path is returned in a more standardized way.

    I.e., folder paths with potential \\ are replaced with /. In addition, if
    a path does not end with / will get an added /.

    Argument
    ----------
    * folder: str ~ The folder path that shall be standardized.
    """
    # Standardize for \ or / as path separator character.
    folder = folder.replace("\\", "/")

    # If the last character is not a path separator, it is
    # added so that all standardized folder path strings
    # contain it.
    if folder[-1] != "/":
        folder += "/"

    return folder

def convert_to_irreversible(model):
    """Split reversible reactions into two irreversible reactions

    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    Arguments
    ----------
    * model: cobra.Model ~ A Model object which will be modified in place.

    """
    #warn("deprecated, not applicable for optlang solvers", DeprecationWarning)
    reactions_to_add = []
    coefficients = {}
    for reaction in model.reactions:
        if reaction.lower_bound < 0 and reaction.upper_bound == 0:
            for metabolite in reaction.metabolites:
                original_coefficient = reaction.get_coefficient(metabolite)
                reaction.add_metabolites({metabolite: -2*original_coefficient})
            reaction.id += "_reverse"
            reaction.upper_bound = -reaction.lower_bound
            reaction.lower_bound = 0
        # If a reaction is reverse only, the forward reaction (which
        # will be constrained to 0) will be left in the model.
        if reaction.lower_bound < 0 and reaction.upper_bound > 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
            reverse_reaction.upper_bound = -reaction.lower_bound
            coefficients[
                reverse_reaction] = reaction.objective_coefficient * -1
            reaction.lower_bound = max(0, reaction.lower_bound)
            reaction.upper_bound = max(0, reaction.upper_bound)
            # Make the directions aware of each other
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {k: v * -1
                             for k, v in reaction._metabolites.items()}
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction.subsystem = reaction.subsystem
            reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    model.add_reactions(reactions_to_add)
    set_objective(model, coefficients, additive=True)

    
def get_genes_and_gpr(model,gene_outfile,gpr_outfile):
    """
    Retrieve genes and gene_reaction_rule from a Genome-Scale Metabolic Network Model.

    Args:
        model (cobra.Model): A genome-scale metabolic network model.
        gene_outfile (str): Path to the output file to save the genes.
        gpr_outfile (str): Path to the output file to save the gene_reaction_rule.

    Returns:
        list: A list containing two pandas DataFrames - genes and gpr.

    """
    model_dict = model_to_dict(model, sort=False)
    genes = pd.DataFrame(model_dict['genes']).set_index(['id'])
    genes.to_csv(gene_outfile)
    all_gpr = pd.DataFrame(model_dict['reactions']).set_index(['id'])
    all_gpr.to_csv(gpr_outfile)
    return [genes, all_gpr]

def isoenzyme_split(model):
    """Split isoenzyme reaction to mutiple reaction

    Arguments
    ----------
    * model: cobra.Model.
    
    :return: new cobra.Model.
    """  
    for r in model.reactions:
        if re.search(" or ", r.gene_reaction_rule):
            rea = r.copy()
            gene = r.gene_reaction_rule.split(" or ")
            for index, value in enumerate(gene):
                if index == 0:
                    r.id = r.id + "_num1"
                    r.gene_reaction_rule = value
                else:
                    r_add = rea.copy()
                    r_add.id = rea.id + "_num" + str(index+1)
                    r_add.gene_reaction_rule = value
                    model.add_reaction(r_add)
    for r in model.reactions:
        r.gene_reaction_rule = r.gene_reaction_rule.strip("( )")
    return model
        
def get_reaction_mw(sbml_path,project_folder,project_name,json_output_file,enzyme_unit_number_file):
    """
    Adds proteomic constraints according to sMOMENT to the given stoichiometric model and stores it as SBML.

    Arguments:
    * sbml_path: str - Path to the SBML or JSON file representing the model.
    * project_folder: str - The folder in which the spreadsheets and JSONs with the model's supplemental
      data can be found.
    * project_name: str - The sMOMENTed model creation's name, which will be added at the beginning
      of the created SBML's name.
    * json_output_file: str - Path to the output JSON file to store the reaction molecular weights.
    * enzyme_unit_number_file: str - Path to the enzyme unit number file.

    Output:
    ----------
    Writes the reaction molecular weights to a JSON file.
    """
    if re.search('\.xml',sbml_path):
        model = cobra.io.read_sbml_model(sbml_path)
    elif re.search('\.json',sbml_path):
        model = cobra.io.json.load_json_model(sbml_path)
    basepath: str = project_folder + project_name
    # READ REACTIONS<->KEGG ID XLSX
    protein_id_mass_mapping: Dict[str, float] = json_load(
        basepath + "_protein_id_mass_mapping.json")
    if enzyme_unit_number_file != 'none':
        enzyme_unit_number=pd.read_csv(enzyme_unit_number_file,index_col=0)
    convert_to_irreversible(model)
    #split isoenzyme
    model = isoenzyme_split(model)

    #subunit_num	1 and 1 and 1 
    reaction_mw={}
    for r in model.reactions:
            if re.search('_num',r.id):
                r_id=r.id.split('_num')[0]
            else:
                r_id=r.id            
        #if r_id in reactions_kcat_mapping_database.keys():
            #print(r.id,r.gene_reaction_rule)
            mass_sum = .0
            if re.search(' and ',r.gene_reaction_rule):
                genelist=r.gene_reaction_rule.split(' and ')
                for eachgene in genelist:
                    #enzyme_unit_number=1
                    if eachgene in protein_id_mass_mapping.keys():
                        #mass_sum += protein_id_mass_mapping[eachgene]['mw'] * enzyme_unit_number
                        if enzyme_unit_number_file != 'none' :
                            # Find the row where the ID exists in the 'geneNames' column
                            row = enzyme_unit_number[enzyme_unit_number['geneNames'].str.contains(eachgene, na=False)]
                            if not row.empty:
                                # ID found, output the corresponding subunit number
                                subunit_number = row['subunitnumber'].values[0]
                            else:
                                subunit_number = 1
                            mass_sum += protein_id_mass_mapping[eachgene]['mw'] * int(subunit_number)
                        else:
                            mass_sum += protein_id_mass_mapping[eachgene]['mw']
                #print(mass_sum)
                reaction_mw[r.id]=mass_sum
            else:  # Single enzyme
                eachgene=r.gene_reaction_rule
                #enzyme_unit_number = 1
                if eachgene in protein_id_mass_mapping.keys():
                    #print(protein_id_mass_mapping[eachgene] * enzyme_unit_number)
                    #reaction_mw[r.id]=protein_id_mass_mapping[eachgene]['mw'] * enzyme_unit_number
                    if enzyme_unit_number_file != 'none':
                        row = enzyme_unit_number[enzyme_unit_number['geneNames'].str.contains(eachgene, na=False)]
                        if not row.empty:
                            # ID found, output the corresponding subunit number
                            subunit_number = row['subunitnumber'].values[0]
                        else:
                            subunit_number = 1
                        reaction_mw[r.id]=protein_id_mass_mapping[eachgene]['mw'] * int(subunit_number)
                    else:
                        reaction_mw[r.id]=protein_id_mass_mapping[eachgene]['mw']
    json_write(json_output_file, reaction_mw)        
        
def get_reaction_kcat_mw(model,project_folder, project_name,type_of_default_kcat_selection,enzyme_unit_number_file,json_output_file):
    """Adds proteomic constraints according to sMOMENT to the given stoichiometric model and stores it as SBML.

    Arguments:
    ----------
    * model: cobra.Model ~ A cobra Model representation of the metabolic network. This model will
      be changed using cobrapy functions in order to add the proteomic constraints.
    * project_folder: str ~ The folder in which the spreadsheets and JSONs with the model's supplemental
      data can be found.
    * project_name: str ~ The sMOMENTed model creation's name, which will be added at the beginning
      of the created SBML's name.
    * type_of_default_kcat_selection: str ~ The type of default kcat selection to use. Available options
      are "median", "mean", "max", "random", or "Null".
    * enzyme_unit_number_file: str ~ The file path of the enzyme unit number data.
    * json_output_file: str ~ The file path to save the output as a JSON file.

    Output:
    ----------
    None.
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)

    # Set folder path for newly created SBML and name for the reaction ID addition (added at the end,
    # and used in order to have a programatically convinient way to separate additions such as 'reverse'
    # from the 'actual' reaction ID).
    basepath: str = project_folder + project_name
    id_addition: str = "_num"
    # READ REACTIONS<->KEGG ID XLSX
    protein_id_mass_mapping: Dict[str, float] = json_load(
        basepath + "_protein_id_mass_mapping.json")
    if enzyme_unit_number_file != 'none':
        enzyme_unit_number=pd.read_csv(enzyme_unit_number_file,index_col=0)
    # Make model irreversible, separating all reversible reactions to which a gene rule is given
    # in order to save some reactions.
    convert_to_irreversible(model)
    #split isoenzyme
    model = isoenzyme_split(model)
    
    # Read reaction <-> kcat mapping :-)
    reactions_kcat_mapping_database = json_load(
        basepath + "_reactions_kcat_mapping_combined.json")

    # sMOMENT :D
    # Get all kcats which are not math.nan and calculate the median of them, which will be used as default kcat
    all_kcats = [x["forward"] for x in reactions_kcat_mapping_database.values()] + \
                [x["reverse"] for x in reactions_kcat_mapping_database.values()]
    all_kcats = [x for x in all_kcats if not math.isnan(x)]

    if type_of_default_kcat_selection == "median":
        default_kcat = statistics.median(all_kcats)
    elif type_of_default_kcat_selection == "mean":
        default_kcat = statistics.mean(all_kcats)
    elif type_of_default_kcat_selection == "max":
        default_kcat = np.max(all_kcats)
    elif type_of_default_kcat_selection == "random":
        default_kcat = random.choice(all_kcats)
    else:
        default_kcat = 'Null'

    print(f"Default kcat is: {default_kcat}")

    # Get all reaction IDs of the given model
    model_reaction_ids = [x.id for x in model.reactions]

    # Main loop :D, add enzyme constraints to reactions \o/
    reaction_kcat_mw={}
    for model_reaction_id in model_reaction_ids:
        # Get the reaction and split the ID at the ID addition
        reaction = model.reactions.get_by_id(model_reaction_id)
        splitted_id = reaction.id.split(id_addition)

        # If the reaction has no name, ignore it
        if splitted_id[0] == "":
            continue
        # Take the reaction ID from the first part of the split
        reaction_id = splitted_id[0]
        # Remove GPRSPLIT name addition from reactions with measured protein concentrations
        if "_GPRSPLIT_" in reaction_id:
            reaction_id = reaction_id.split("_GPRSPLIT_")[0]

        # Retrieve the reaction's forward and reverse kcats from the given reaction<->kcat database
        if re.search('_num',reaction_id):
            reaction_id=reaction_id.split('_num')[0]
        if re.search('_reverse',reaction_id):
            reaction_id=reaction_id.split('_reverse')[0]
        if reaction_id not in reactions_kcat_mapping_database.keys():
            continue
        forward_kcat = reactions_kcat_mapping_database[reaction_id]["forward"]
        reverse_kcat = reactions_kcat_mapping_database[reaction_id]["reverse"]

        # If the given reaction<->kcat database contains math.nan as the reaction's kcat,
        # set the default kcat as math.nan means that no kcat could be found.
        if math.isnan(forward_kcat):
            forward_kcat = default_kcat
        if math.isnan(reverse_kcat):
            reverse_kcat = default_kcat

        # Add the given forward or reverse kcat is the reaction was
        # splitted due to its reversibility.
        # If the reaction is not splitted, add the forward kcat (this
        # is the only possible direction for non-splitted=non-reversible
        # reactions)
        if model_reaction_id.endswith(id_addition + "forward"):
            reaction_kcat = forward_kcat
        elif model_reaction_id.endswith(id_addition + "reverse"):
            reaction_kcat = reverse_kcat
        else:
            reaction_kcat = forward_kcat

        reaction_kcat_mw[model_reaction_id]={}
        if reaction_kcat=='Null':
            continue
        reaction_kcat_mw[model_reaction_id]['kcat']=reaction_kcat
        if "data_type" in reactions_kcat_mapping_database[reaction_id].keys():
            reaction_kcat_mw[model_reaction_id]['data_type']=reactions_kcat_mapping_database[reaction_id]["data_type"]
        else:
            reaction_kcat_mw[model_reaction_id]['data_type']='fill'

        #MW
        #subunit_num	1 and 1 and 1 
        reaction_mw={}
        mass_sum = .0
        if re.search(' and ',reaction.gene_reaction_rule):
            genelist=reaction.gene_reaction_rule.split(' and ')
            for eachgene in genelist:
                #enzyme_unit_number=1
                if eachgene in protein_id_mass_mapping.keys():
                    #mass_sum += protein_id_mass_mapping[eachgene]['mw'] * enzyme_unit_number
                    if enzyme_unit_number_file != 'none':
                        # Find the row where the ID exists in the 'geneNames' column
                        row = enzyme_unit_number[enzyme_unit_number['geneNames'].str.contains(eachgene, na=False)]
                        if not row.empty:
                            # ID found, output the corresponding subunit number
                            subunit_number = row['subunitnumber'].values[0]
                        else:
                            subunit_number = 1
                        mass_sum += protein_id_mass_mapping[eachgene]['mw'] * int(subunit_number)
                    else:
                        mass_sum += protein_id_mass_mapping[eachgene]['mw']
            #print(mass_sum)
            if mass_sum>0:
                reaction_mw[reaction.id]=mass_sum
        else:  # Single enzyme
            eachgene=reaction.gene_reaction_rule
            #enzyme_unit_number = 1
            if eachgene in protein_id_mass_mapping.keys():
                #print(protein_id_mass_mapping[eachgene] * enzyme_unit_number)
                #reaction_mw[r.id]=protein_id_mass_mapping[eachgene]['mw'] * enzyme_unit_number
                if enzyme_unit_number_file != 'none':
                    # Find the row where the ID exists in the 'geneNames' column
                    row = enzyme_unit_number[enzyme_unit_number['geneNames'].str.contains(eachgene, na=False)]
                    if not row.empty:
                        # ID found, output the corresponding subunit number
                        subunit_number = row['subunitnumber'].values[0]
                    else:
                        subunit_number = 1
                    reaction_mw[reaction.id]=protein_id_mass_mapping[eachgene]['mw'] * int(subunit_number)
                else:
                    reaction_mw[reaction.id]=protein_id_mass_mapping[eachgene]['mw']
        #print(model_reaction_id,reaction_mw.keys())
        if model_reaction_id in reaction_mw.keys():
            #print(model_reaction_id,reaction_mw[model_reaction_id])
            reaction_kcat_mw[model_reaction_id]['MW']=reaction_mw[model_reaction_id]#Kda
            reaction_kcat_mw[model_reaction_id]['kcat_MW']=reaction_kcat_mw[model_reaction_id]['kcat']*3600000/reaction_mw[model_reaction_id]
    reaction_kcat_mw_df = pd.DataFrame(reaction_kcat_mw)
    reaction_kcat_mw_df_T=reaction_kcat_mw_df.T
    reaction_kcat_mw_df_T_select=reaction_kcat_mw_df_T[abs(reaction_kcat_mw_df_T['kcat_MW'])>0]
    reaction_kcat_mw_df_T_select.to_csv(json_output_file)
    #reaction_kcat_mw_df_T.to_csv(project_folder + 'reaction_kcat_MW_total.csv')
    #return reaction_kcat_mw_df_T
    
def isoenzyme_split(model):
    """Split isoenzyme reaction to mutiple reaction

    Arguments
    ----------
    * model: cobra.Model.
    
    :return: new cobra.Model.
    """  
    for r in model.reactions:
        if re.search(" or ", r.gene_reaction_rule):
            rea = r.copy()
            gene = r.gene_reaction_rule.split(" or ")
            for index, value in enumerate(gene):
                if index == 0:
                    r.id = r.id + "_num1"
                    r.gene_reaction_rule = value
                else:
                    r_add = rea.copy()
                    r_add.id = rea.id + "_num" + str(index+1)
                    r_add.gene_reaction_rule = value
                    model.add_reaction(r_add)
    for r in model.reactions:
        r.gene_reaction_rule = r.gene_reaction_rule.strip("( )")
    return model

def trans_model2enz_json_model_split_isoenzyme(model_file, reaction_kcat_mw_file, f, ptot, sigma, lowerbound, upperbound, json_output_file):
    """Tansform cobra model to json mode with  
    enzyme concentration constraintat.

    Arguments
    ----------
    * model_file:   The path of sbml model
    * reaction_kcat_mw_file: The path of storing kcat/MW value of the enzyme catalyzing each
     reaction in the GEM model
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint. 

    """
    if re.search('\.xml',model_file):
        model = cobra.io.read_sbml_model(model_file)
    elif re.search('\.json',model_file):
        model = cobra.io.json.load_json_model(model_file)
    convert_to_irreversible(model)
    model = isoenzyme_split(model)
    model_name = model_file.split('/')[-1].split('.')[0]
    json_path = "./model/%s_irreversible.json" % model_name
    cobra.io.save_json_model(model, json_path)
    dictionary_model = json_load(json_path)
    dictionary_model['enzyme_constraint'] = {'enzyme_mass_fraction': f, 'total_protein_fraction': ptot,
                                             'average_saturation': sigma, 'lowerbound': lowerbound, 'upperbound': upperbound}
    # Reaction-kcat_mw file.
    # eg. AADDGT,49389.2889,40.6396,1215.299582180927
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)
    for eachreaction in range(len(dictionary_model['reactions'])):
        reaction_id = dictionary_model['reactions'][eachreaction]['id']
        #if re.search('_num',reaction_id):
        #    reaction_id=reaction_id.split('_num')[0]
        if reaction_id in reaction_kcat_mw.index:
            dictionary_model['reactions'][eachreaction]['kcat'] = reaction_kcat_mw.loc[reaction_id, 'kcat']
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
        else:
            dictionary_model['reactions'][eachreaction]['kcat'] = ''
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = ''
    json_write(json_output_file, dictionary_model)
    
def get_enzyme_constraint_model(json_model_file):
    """using enzyme concentration constraint
    json model to create a COBRApy model.

    Arguments
    ----------
    * json_model_file: json Model file.

    :return: Construct an enzyme-constrained model.
    """

    dictionary_model = json_load(json_model_file)
    model = cobra.io.json.load_json_model(json_model_file)

    coefficients = dict()
    for rxn in model.reactions:
        for eachr in dictionary_model['reactions']:
            if rxn.id == eachr['id']:
                if eachr['kcat_MW']:
                    coefficients[rxn.forward_variable] = 1 / float(eachr['kcat_MW'])
                break

    lowerbound = dictionary_model['enzyme_constraint']['lowerbound']
    upperbound = dictionary_model['enzyme_constraint']['upperbound']
    constraint = model.problem.Constraint(0, lb=lowerbound, ub=upperbound)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    return model

def get_enzyme_constraint_model_percent(json_model_file,percent):
    """using enzyme concentration constraint
    json model to create a COBRApy model.

    Arguments
    ----------
    * json_model_file: json Model file.

    :return: Construct an enzyme-constrained model.
    """

    dictionary_model = json_load(json_model_file)
    model = cobra.io.json.load_json_model(json_model_file)

    coefficients = dict()
    for rxn in model.reactions:
        for eachr in dictionary_model['reactions']:
            if rxn.id == eachr['id']:
                if eachr['kcat_MW']:
                    coefficients[rxn.forward_variable] = 1 / float(eachr['kcat_MW'])
                break

    lowerbound = dictionary_model['enzyme_constraint']['lowerbound']
    upperbound = dictionary_model['enzyme_constraint']['upperbound']*percent
    constraint = model.problem.Constraint(0, lb=lowerbound, ub=upperbound)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    return model

def get_fluxes_detail_in_model(model,model_pfba_solution,fluxes_outfile,json_model_file):
    """
    Get the detailed information of each reaction.

    Arguments:
    * model: cobra.Model - The metabolic model.
    * model_pfba_solution: pandas.Series - The pFBA solution containing reaction fluxes.
    * fluxes_outfile: str - Path to the output file for reaction fluxes.
    * json_model_file: str - Path to the JSON model file.

    Returns:
    * model_pfba_solution_detail: pandas.DataFrame - Detailed information of each reaction.
    """

    dictionary_model = json_load(json_model_file)
    model_pfba_solution = model_pfba_solution.to_frame()
    model_pfba_solution_detail = pd.DataFrame()
    for index, row in model_pfba_solution.iterrows():
        reaction_detail = model.reactions.get_by_id(index)
        model_pfba_solution_detail.loc[index, 'fluxes'] = row['fluxes']
        for eachreaction in dictionary_model['reactions']:
            if index ==eachreaction['id']:
                if 'annotation' in eachreaction.keys():
                    if 'ec-code' in eachreaction['annotation'].keys():
                        if isinstance (eachreaction['annotation']['ec-code'],list):
                            model_pfba_solution_detail.loc[index, 'ec-code'] = (',').join(eachreaction['annotation']['ec-code'])
                        else:
                            model_pfba_solution_detail.loc[index, 'ec-code'] = eachreaction['annotation']['ec-code']    
                if 'kcat_MW' in eachreaction.keys():
                    if eachreaction['kcat_MW']:
                        model_pfba_solution_detail.loc[index, 'kcat_MW'] = eachreaction['kcat_MW']
                        model_pfba_solution_detail.loc[index, 'E'] = float(row['fluxes'])/float(eachreaction['kcat_MW'])
                break
        model_pfba_solution_detail.loc[index, 'equ'] = reaction_detail.reaction
    model_pfba_solution_detail.to_csv(fluxes_outfile)
    return model_pfba_solution_detail

def json_write(path, dictionary):
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path:   The path of the JSON file that shall be written
    * dictionary: The dictionary which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(dictionary, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)
        
def GENENAME_2_ACC_from_uniprot(query,outfile):
    '''
    get amin acid sequence and mass by Uniprot id
    
    Arguments
    ----------
    query: Uniprot ID list.
    outfile:  The path of the ACC file that shall be written
    '''
    #print(' '.join(query).replace('511145.',''))
    url = 'https://legacy.uniprot.org/uploadlists/'
    params = {
        'from': 'GENENAME',
        'to': 'ACC',
        'format': 'tab',
        'query': ' '.join(query),
        'columns':'id,entry name,protein names,genes,organism,ec,mass,database(PDB)'
    }
    data = urlencode(params).encode()
    request = Request(url, data)
    # Please set your email address here to help us debug in case of problems.
    contact = ""
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urlopen(request)
    page = response.read()
    outFile = open(outfile,'w') 
    namesRegex = re.compile(r'yourlist:(.*)\n')
    outFile.write(namesRegex.sub('Gene ID\n',page.decode('utf-8')))
    #print(namesRegex.sub('Protein AC\t',page.decode('utf-8')))
    outFile.close()

def calculate_f(model_file, gene_abundance_file, gene_mw_file, gene_mw_colname, gene_abundance_colname):
    """
    Calculate the enzyme mass fraction (f) based on protein abundance and molecular weight.

    Arguments:
    * model_file: str - File path of the SBML model.
    * gene_abundance_file: str - File path of the gene abundance data file.
    * gene_mw_file: str - File path of the gene molecular weight data file.
    * gene_mw_colname: str - Name of the column in the gene molecular weight file containing the molecular weight values.
    * gene_abundance_colname: str - Name of the column in the gene abundance file containing the abundance values.

    Returns:
    * f: float - The enzyme mass fraction.
    """

    if re.search('\.xml',model_file):
        model = cobra.io.read_sbml_model(model_file)
    elif re.search('\.json',model_file):
        model = cobra.io.json.load_json_model(model_file)
    model_gene_list = [eachg.id for eachg in model.genes]
    uni_model_gene_list = list(set(model_gene_list))

    gene_abundance = pd.read_csv(gene_abundance_file, index_col=0)
    gene_list = list(set(gene_abundance.index))
    GENENAME_2_ACC_from_uniprot(gene_list,gene_mw_file)
    gene_mw = pd.read_csv(gene_mw_file, sep='\t', index_col=gene_mw_colname)
    enzy_abundance = 0
    pro_abundance = 0

    for gene_i in gene_abundance.index:
        if gene_i in gene_mw.index:
            mass_value = gene_mw.loc[gene_i, 'Mass']
            if isinstance(mass_value, str):
                abundance = gene_abundance.loc[gene_i, gene_abundance_colname] * int(mass_value.replace(',', ''))
            else:
                abundance = gene_abundance.loc[gene_i, gene_abundance_colname] * int(mass_value[0].replace(',', ''))
            pro_abundance += abundance
            if gene_i in uni_model_gene_list:
                enzy_abundance += abundance

    f = enzy_abundance / pro_abundance
    return f

def calculate_f_v2(model_file, gene_abundance_file, gene_abundance_colname, taxonom_id):
    '''
    Calculate the enzyme mass fraction (f) based on protein abundance.

    Arguments:
    * model_file: str - File path of the SBML model.
    * gene_abundance_file: str - File path of the gene abundance data file.
    * gene_abundance_colname: str - Name of the column in the gene abundance file containing the abundance values.
    * taxonom_id: str - Taxonomy ID of the organism.

    Returns:
    * f: float - The enzyme mass fraction.
    '''
    if re.search('\.xml',model_file):
        model = cobra.io.read_sbml_model(model_file)
    elif re.search('\.json',model_file):
        model = cobra.io.json.load_json_model(model_file)
    model_gene_list = [eachg.id for eachg in model.genes]
    uni_model_gene_list = list(set(model_gene_list))

    gene_abundance = pd.read_csv(gene_abundance_file, index_col=0)
    enzy_abundance = 0
    pro_abundance = 0

    for gene_i in probar(gene_abundance.index):
        uniprot_query_url = f"https://rest.uniprot.org/uniprotkb/search?query={gene_i}+AND+organism_id:{taxonom_id}&format=tsv&fields=accession,mass"
        try:
            uniprot_mass = requests.get(uniprot_query_url).text.split("\n")[1].split("\t")[1]
        except:
            print(gene_i + ' has no mass!')
        else:
            abundance = gene_abundance.loc[gene_i, gene_abundance_colname] * int(uniprot_mass)
            pro_abundance += abundance
            if gene_i in uni_model_gene_list:
                enzy_abundance += abundance

    f = enzy_abundance / pro_abundance
    return f

def get_model_substrate_obj(use_model):
    '''
    change model substrate for single carbon source
    
    Arguments
    ----------
    use_model: cobra.Model ~ A Model object which will be modified in place.
    '''
    
    ATPM='No' 
    substrate_list=[]
    concentration_list=[]
    EX_exclude_reaction_list=['EX_pi_e','EX_h_e','EX_fe3_e','EX_mn2_e','EX_co2_e','EX_fe2_e','EX_h2_e','EX_zn2_e',\
                             'EX_mg2_e','EX_ca2_e','EX_so3_e','EX_ni2_e','EX_no_e','EX_cu2_e','EX_hg2_e','EX_cd2_e',\
                             'EX_h2o2_e','EX_h2o_e','EX_no2_e','EX_nh4_e','EX_so4_e','EX_k_e','EX_na1_e','EX_o2_e',\
                             'EX_o2s_e','EX_ag_e','EX_cu_e','EX_so2_e','EX_cl_e','EX_n2o_e','EX_cs1_e','EX_cobalt2_e']
    EX_exclude_reaction_list=EX_exclude_reaction_list+[i+'_reverse' for i in EX_exclude_reaction_list]
    for r in use_model.reactions:
        if r.objective_coefficient == 1:
            obj=r.id #Product name
        #elif not r.lower_bound==0 and not r.lower_bound==-1000 and not r.lower_bound==-999999 and abs(r.lower_bound)>0.1:#排除很小的值
        elif not r.upper_bound==0 and not r.upper_bound==1000 and not r.upper_bound==999999 and abs(r.upper_bound)>0.1:#排除很小的值
            #print(r.id,r.upper_bound,r.lower_bound)
            if r.id=='ATPM':
                if r.upper_bound>0:
                    ATPM='Yes' #ATP maintenance requirement
            elif r.id not in EX_exclude_reaction_list:
                #print(r.id,r.upper_bound,r.lower_bound)
                #substrate=r.id #Substrate name
                substrate_list.append(r.id)
                #concentration=r.upper_bound #Substrate uptake rate  
                concentration_list.append(r.upper_bound)
    return(obj,substrate_list,concentration_list,ATPM)

def parse_sabio_rk_for_eclist(ec_numbers_list: List[str], json_output_path: str, bigg_id_name_mapping_path: str) -> None:
    """Retrieves kcats from SABIO-RK for the given model and stores it in a JSON for the given model in the given path.

    Algorithm
    ----------
    Using the SABIO-RK REST API (as of 2019/30/04, it is explained under
    http://sabiork.h-its.org/layouts/content/docuRESTfulWeb/RESTWebserviceIntro.gsp),


    Arguments
    ----------
    * eclist: List[str] ~ eclist.
    * json_output_path: str ~ The path of the JSON that shall be created

    Output
    ----------
    * A JSON in the given project folder with the following structure:
    <pre>
        {
            "$EC_NUMBER_OR_KEGG_REACTION_ID": {
                "$SUBSTRATE_WITH_BIGG_ID_1": {
                    "$ORGANISM_1": [
                        $kcat_1,
                        (...)
                        $kcat_n,
                    ]
                },
                (...),
                "REST": {
                    "$ORGANISM_1": [
                        $kcat_1,
                        (...)
                        $kcat_n,
                    ]
                }
            }
            (...),
        }
    </pre>
    'REST' stands for a substrate without found BIGG ID.
    """
    # GET KCATS FOR EC NUMBERS
    ec_number_kcat_mapping = get_ec_number_kcats_wildcard_search(ec_numbers_list, bigg_id_name_mapping_path)
    json_write(json_output_path, ec_number_kcat_mapping)
    
def get_protein_mass_mapping_from_local(sbml_path: str, project_folder: str, project_name: str,uniprot_data_file: str) -> None:
    """Returns a JSON with a mapping of protein IDs as keys, and as values the protein mass in kDa.

    The protein masses are calculated using the amino acid sequence from UniProt (retrieved using
    UniProt's REST API).

    Arguments
    ----------
    * model: cobra.Model ~ The model in the cobrapy format
    * project_folder: str ~ The folder in which the JSON shall be created
    * project_name: str ~ The beginning of the JSON's file name
    * uniprot_data_file: str ~ The gene information obtained from uniprot
    Output
    ----------
    A JSON file with the path project_folder+project_name+'_protein_id_mass_mapping.json'
    and the following structure:
    <pre>
    {
        "$PROTEIN_ID": $PROTEIN_MASS_IN_KDA,
        (...),
    }
    </pre>
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)
    # The beginning of the created JSON's path :D
    basepath: str = project_folder + project_name

    # Read the model from the file
    if re.search('\.xml', sbml_path):
        model = cobra.io.read_sbml_model(sbml_path)
    elif re.search('\.json', sbml_path):
        model = cobra.io.json.load_json_model(sbml_path)
    # GET UNIPROT ID - PROTEIN MAPPING
    uniprot_id_protein_id_mapping: Dict[str, List[str]] = {}
    for gene in model.genes:
        # Without a UniProt ID, no mass mapping can be found
        if "uniprot" not in gene.annotation:
            continue
        uniprot_id = gene.annotation["uniprot"]
        if uniprot_id in uniprot_id_protein_id_mapping.keys():
            uniprot_id_protein_id_mapping[uniprot_id].append(gene.id)
        else:
            uniprot_id_protein_id_mapping[uniprot_id] = [gene.id]

    # GET UNIPROT ID<->PROTEIN MASS MAPPING
    uniprot_id_protein_mass_mapping = json_load(uniprot_data_file)
    
    # Create the final protein ID <-> mass mapping
    protein_id_mass_mapping: Dict[str, float] = {}
    for uniprot_id in list(uniprot_id_protein_mass_mapping.keys()):
        try:
            protein_ids = uniprot_id_protein_id_mapping[uniprot_id]
        except Exception:
            #print(f"No mass found for {uniprot_id}!")
            continue
        for protein_id in protein_ids:
            protein_id_mass_mapping[protein_id] = uniprot_id_protein_mass_mapping[uniprot_id]

    # Write protein mass list JSON :D
    #print("Protein ID<->Mass mapping done!")
    json_write(basepath+"_protein_id_mass_mapping.json", protein_id_mass_mapping)
        
def change_reaction_kcat_by_database(json_model_path,select_reactionlist, EC_max_file, reaction_kcat_mw_file, reaction_kapp_change_file):
    """Use the kcat in database to change reaction kcat in model

    Arguments
    ----------
    * json_model_path: The file storing json model.
    * select_reactionlist: reaction list need to change.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_kapp_change_file: changed file stored reaction kcat/mw.

    :return: a dataframe stored new reaction kcat/mw .
    """
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)
    Brenda_sabio_combined_select = json_load(EC_max_file)
    
    json_model=cobra.io.load_json_model(json_model_path)
    reaction_change_accord_database = []
    for eachreaction in select_reactionlist:
        select_reaction = json_model.reactions.get_by_id(eachreaction)
        if "ec-code" in select_reaction.annotation.keys():
            ec_number = select_reaction.annotation["ec-code"]
            kcat_max_list = []
            if isinstance(ec_number, str):
                if ec_number in Brenda_sabio_combined_select.keys():
                    reaction_kcat_max = Brenda_sabio_combined_select[ec_number]['kcat_max']
                    if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max :
                        reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max #h_1
                        reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                        reaction_change_accord_database.append(eachreaction) 
            else:
                for eachec in ec_number:
                    if eachec in Brenda_sabio_combined_select.keys():
                        kcat_max_list.append(Brenda_sabio_combined_select[eachec]['kcat_max'])
                reaction_kcat_max = np.max(kcat_max_list)     
                if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max:
                    reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max
                    reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                    reaction_change_accord_database.append(eachreaction)    
                
    reaction_kcat_mw.to_csv(reaction_kapp_change_file)
    return(reaction_change_accord_database)

def get_enz_model_use_enz_usage_by_eckcat(enz_ratio,json_model_path, reaction_flux_file, EC_max_file, reaction_kcat_mw_file, f, \
                                          ptot, sigma, lowerbound, upperbound, json_output_file, reaction_mw_outfile):
    """Get new enzyme model using enzyme usage to calibration

    Arguments
    ----------
    * enz_ratio: enzyme ratio which needed change.
    * json_model_path: The file storing json model.
    * reaction_flux_file: reaction-flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_enz_usage_file: enzyme usage of each reaction.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * model_file: cobra model.
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint.  
    * json_output_file: json file store json model
    * reaction_mw_outfile: changed file stored reaction kcat/mw.

    :return: new enzyme model
    """ 
    reaction_fluxes = pd.read_csv(reaction_flux_file, index_col=0)
    reaction_fluxes['enz ratio'] = reaction_fluxes['E']/np.sum(reaction_fluxes['E'])
    reaction_fluxes=reaction_fluxes.sort_values(by="enz ratio", axis=0, ascending=False)
    
    select_reaction = list(
        reaction_fluxes[reaction_fluxes['enz ratio'] > enz_ratio].index)  # more than 1%
    print('need changing reaction: ')
    print(select_reaction)
    change_reaction_list_round1 = change_reaction_kcat_by_database(json_model_path,select_reaction, EC_max_file, reaction_kcat_mw_file, reaction_mw_outfile)
    print('changed reaction: ')
    print(change_reaction_list_round1)

    trans_model2enz_json_model_split_isoenzyme(
        json_model_path, reaction_mw_outfile, f, ptot, sigma, lowerbound, upperbound, json_output_file)

    enz_model = get_enzyme_constraint_model(json_output_file)
    return enz_model

def adj_reaction_kcat_by_database(json_model_path,select_reactionlist, need_change_reaction_list, changed_reaction_list,EC_max_file, reaction_kcat_mw):
    """Use the kcat in database to change reaction kcat in model

    Arguments
    ----------
    * json_model_path: The file storing json model.
    * select_reactionlist: reaction list need to change.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_kapp_change_file: changed file stored reaction kcat/mw.

    :return: a dataframe stored new reaction kcat/mw .
    """
    Brenda_sabio_combined_select = json_load(EC_max_file)
    
    json_model=cobra.io.load_json_model(json_model_path)
    for eachreaction in select_reactionlist:
        need_change_reaction_list.append(eachreaction)
        if eachreaction in list(reaction_kcat_mw.index):
            select_reaction = json_model.reactions.get_by_id(eachreaction)
            if "ec-code" in select_reaction.annotation.keys():
                ec_number = select_reaction.annotation["ec-code"]
                kcat_max_list = []
                if isinstance(ec_number, str):
                    if ec_number in Brenda_sabio_combined_select.keys():
                        reaction_kcat_max = Brenda_sabio_combined_select[ec_number]['kcat_max']
                        if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max:
                            reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max#h_1
                            reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                            changed_reaction_list.append(eachreaction) 
                else:
                    for eachec in ec_number:
                        if eachec in Brenda_sabio_combined_select.keys():
                            kcat_max_list.append(Brenda_sabio_combined_select[eachec]['kcat_max'])
                    if len(kcat_max_list)>0:        
                        reaction_kcat_max = np.max(kcat_max_list)  
                    else:
                        reaction_kcat_max = 0  
                    if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max:
                        reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max
                        reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                        changed_reaction_list.append(eachreaction)    
                
    return(need_change_reaction_list,changed_reaction_list,reaction_kcat_mw)

def adj_trans_model2enz_model(model_file, reaction_kcat_mw, f, ptot, sigma, lowerbound, upperbound, json_output_file):
    """Tansform cobra model to json mode with  
    enzyme concentration constraintat.

    Arguments
    ----------
    * model_file:   The path of sbml model
    * reaction_kcat_mw_file: The path of storing kcat/MW value of the enzyme catalyzing each
     reaction in the GEM model
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint. 

    """
    if re.search('\.xml',model_file):
        model = cobra.io.read_sbml_model(model_file)
    elif re.search('\.json',model_file):
        model = cobra.io.json.load_json_model(model_file)
    convert_to_irreversible(model)
    model = isoenzyme_split(model)
    model_name = model_file.split('/')[-1].split('.')[0]
    json_path = "./model/%s_irreversible.json" % model_name
    cobra.io.save_json_model(model, json_path)
    dictionary_model = json_load(json_path)
    dictionary_model['enzyme_constraint'] = {'enzyme_mass_fraction': f, 'total_protein_fraction': ptot,
                                             'average_saturation': sigma, 'lowerbound': lowerbound, 'upperbound': upperbound}
    # Reaction-kcat_mw file.
    # eg. AADDGT,49389.2889,40.6396,1215.299582180927
    for eachreaction in range(len(dictionary_model['reactions'])):
        reaction_id = dictionary_model['reactions'][eachreaction]['id']
        #if re.search('_num',reaction_id):
        #    reaction_id=reaction_id.split('_num')[0]
        if reaction_id in reaction_kcat_mw.index:
            dictionary_model['reactions'][eachreaction]['kcat'] = reaction_kcat_mw.loc[reaction_id, 'kcat']
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
        else:
            dictionary_model['reactions'][eachreaction]['kcat'] = ''
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = ''
    json_write(json_output_file, dictionary_model)
    
def change_enz_model_by_enz_usage(json_model_path, reaction_flux_file, EC_max_file, reaction_kcat_mw, need_change_reaction_list, changed_reaction_list,f, ptot, sigma, lowerbound, upperbound, json_output_file):
    """Get new enzyme model using enzyme usage to calibration

    Arguments
    ----------
    * json_model_path: The file storing json model.
    * reaction_flux_file: reaction-flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_enz_usage_file: enzyme usage of each reaction.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * model_file: cobra model.
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint.  
    * json_output_file: json file store json model
    * reaction_mw_outfile: changed file stored reaction kcat/mw.

    :return: new enzyme model
    """ 
    reaction_fluxes = pd.read_csv(reaction_flux_file, index_col=0)
    reaction_fluxes['enz ratio'] = reaction_fluxes['E']/np.sum(reaction_fluxes['E'])
    reaction_fluxes=reaction_fluxes.sort_values(by="enz ratio", axis=0, ascending=False)

    i=0
    select_reaction = reaction_fluxes.index[0]
    while ((select_reaction in need_change_reaction_list) and (i<len(list(reaction_fluxes.index)))):
        i=i+1
        #print(i)
        select_reaction = reaction_fluxes.index[i+1]
        
    print('Need changing reaction: ')
    print(select_reaction)
    [need_change_reaction_list, changed_reaction_list,reaction_kcat_mw] = adj_reaction_kcat_by_database(json_model_path,[select_reaction], need_change_reaction_list, changed_reaction_list, EC_max_file, reaction_kcat_mw)
    print('Changed reaction: ')
    print(changed_reaction_list)

    adj_trans_model2enz_model(json_model_path, reaction_kcat_mw, f, ptot, sigma, lowerbound, upperbound, json_output_file)

    enz_model = get_enzyme_constraint_model(json_output_file)
    return (enz_model,reaction_kcat_mw,need_change_reaction_list, changed_reaction_list)

def get_met_bigg_id(model):
    '''
    This function is used to get the bigg id of metabolites in the model.

    Arguments
    ----------
    * model: cobra.Model ~ A Model object which will be modified in place.

    Returns
    ----------
    * metdf_unique: pandas.DataFrame ~ A DataFrame containing unique metabolite IDs and names.
    '''

    metlist = []
    metnamelist = []
    for met in model.metabolites:
        if '_c' in str(met):
            metlist.append(str(met.id).split('_c')[0])
            metnamelist.append(str(met.name))
        elif '_e' in str(met):
            metlist.append(str(met.id).split('_e')[0])
            metnamelist.append(str(met.name))
        elif '_p' in str(met):
            metlist.append(str(met.id).split('_p')[0])
            metnamelist.append(str(met.name))
        else:
            metlist.append(str(met.id))
            metnamelist.append(str(met.name))

    metdf = pd.DataFrame()
    metdf_unique = pd.DataFrame()

    metdf['met'] = metlist
    metdf['name'] = metnamelist
    metdf.set_index('met', inplace=True)

    metlist = list(set(metlist))
    metdf_unique['met'] = metlist
    metdf_unique.set_index('met', inplace=True)

    for row in metdf_unique.itertuples():
        if type(metdf.loc[str(row.Index), 'name']) == str:
            metdf_unique.loc[str(row.Index), 'name'] = metdf.loc[str(row.Index), 'name']
        else:
            for i in metdf.loc[str(row.Index), 'name']:
                metdf_unique.loc[str(row.Index), 'name'] = i
    metdf_unique.reset_index(inplace=True)
    # metdf_unique.to_csv('./metbigg_name_unique.csv')
    return metdf_unique

def convert_bigg_met_to_inchikey(metlist, inchikey_list_file):
    '''
    This function is used to convert the bigg id of metabolites in the model to inchikey.

    Arguments
    ----------
    * metlist: list ~ List of metabolite IDs.
    * inchikey_list_file: str ~ File path to save the resulting DataFrame as a CSV file.

    Returns
    ----------
    * inchkeydf: pandas.DataFrame ~ DataFrame containing metabolite names and InChIKeys.
    '''

    inchikeylist = []
    for met in probar(metlist):
        time.sleep(2.0)
        try:
            url = f'http://bigg.ucsd.edu/api/v2/universal/metabolites/{met}'
            response = requests.get(url, headers={"Accept": "application/json"})
            jsonData = json.loads(response.text)
            inchikeylist.append(jsonData['database_links']['InChI Key'][0]['id'])
        except:
            inchikeylist.append('NA')

    inchkeydf = pd.DataFrame({'metname': metlist, 'inchikey': inchikeylist})
    inchkeydf.to_csv(inchikey_list_file, index=False)

    return inchkeydf

def convert_inchikey_to_smiles(inchkeydf, inchikey_list_smilesfile):
    '''
    This function is used to convert the InChIKey of metabolites in the model to SMILES.

    Arguments
    ----------
    * inchkeydf: pandas.DataFrame ~ DataFrame object containing the InChIKey of metabolites.
    * inchikey_list_smilesfile: str ~ File path to save the resulting DataFrame as a CSV file.

    Returns
    ----------
    * inchkeydf: pandas.DataFrame ~ DataFrame with added 'Canonical_SMILES' column containing SMILES.
    '''

    print('Converting...')
    inchkeydf.set_index('metname', inplace=True)
    inchkeydf['Canonical_SMILES'] = 'noinchikey'  # Initialize all values as 'noinchikey'
    valid_inchikeys = inchkeydf[inchkeydf['inchikey'].notna()].index  # Find non-null InChIKeys
    faillist = []

    for eachc in probar(valid_inchikeys):
        try:
            clist = pcp.get_compounds(inchkeydf.loc[eachc, 'inchikey'], 'inchikey')
            if len(clist) == 0:
                inchkeydf.loc[eachc, 'Canonical_SMILES'] = 'noinchikey_inpubchem'
            else:
                compound = clist[0]
                try:
                    inchkeydf.loc[eachc, 'Canonical_SMILES'] = compound.canonical_smiles
                except:
                    inchkeydf.loc[eachc, 'Canonical_SMILES'] = 'nosmiles'
        except:
            faillist.append(eachc)
        time.sleep(2.0)

    print(f'{str(faillist)} try again later')
    print('Fail secure!')

    while True:
        successlist = []
        for eachc in faillist:
            try:
                clist = pcp.get_compounds(inchkeydf.loc[eachc, 'inchikey'], 'inchikey')
                time.sleep(2.0)
                if len(clist) == 0:
                    inchkeydf.loc[eachc, 'Canonical_SMILES'] = 'noinchikey_inpubchem'
                else:
                    compound = clist[0]
                    try:
                        inchkeydf.loc[eachc, 'Canonical_SMILES'] = compound.canonical_smiles
                    except:
                        inchkeydf.loc[eachc, 'Canonical_SMILES'] = 'nosmiles'
                successlist.append(eachc)
                if len(faillist) - len(successlist) > 0:
                    print(f'still have {str(len(faillist) - len(successlist))} fail.')
            except:
                continue
        if len(faillist) == len(successlist):
            break
        else:
            faillist = list(set(faillist) - set(successlist))
            print(f'{str(faillist)} try again later.')

    inchkeydf.reset_index(inplace=True)
    inchkeydf.to_csv(inchikey_list_smilesfile, index=False)

    return inchkeydf

def get_model_protein_sequence_and_mass(model, subbnumdf, prodf_file):
    '''
    This function is used to get the protein sequence and mass of the model.
    
    Arguments
    ----------
    * model: cobra.Model ~ A Model object which will be modified in place.
    * subbnumdf: pandas.DataFrame ~ DataFrame containing enzyme subunit information.
    * prodf_file: str ~ File path to save the resulting DataFrame as a CSV file.

    Returns
    ----------
    * prodf: pandas.DataFrame ~ DataFrame with protein information.
    '''
    genelist = []
    prolist = []
    enzyme_unit_number = subbnumdf
    for gene in model.genes:
        if 'uniprot' in gene.annotation:
            pro = gene.annotation['uniprot']
            genelist.append(gene.id)
            prolist.append(str(pro))
    prodf = pd.DataFrame(columns=['geneid', 'pro'])
    prodf['geneid'] = genelist
    prodf['pro'] = prolist

    def fetch_protein_data(row):
        try:
            query = str(row)
            uniprot_query_url = f"https://rest.uniprot.org/uniprotkb/search?query=accession:{query}&format=tsv&fields=accession,sequence"
            uniprot_data = requests.get(uniprot_query_url).text.split("\n")[1].split("\t")[1]
            return uniprot_data
        except:
            return None

    prodf['aaseq'] = prodf['pro'].apply(fetch_protein_data)
    prodf.dropna(subset=['aaseq'], inplace=True)
    prodf['mass'] = prodf['aaseq'].apply(lambda seq: ProteinAnalysis(seq, monoisotopic=False).molecular_weight())
    prodf['subunitmass'] = prodf.apply(lambda row: row['mass'] * int(enzyme_unit_number.loc[row['pro'], 'subunitnumber'])
                                      if row['pro'] in enzyme_unit_number.index else row['mass'], axis=1)

    prodf.to_csv(prodf_file, index=False)
    return prodf

def get_currency_metabolites():
    '''
    Returns a set of currency metabolites.
    '''
    currencylist1 = ['coa_c', 'co2_c', 'co2_e', 'co2_p', 'cobalt2_c', 'cobalt2_e', 'cobalt2_p', 'h_c', 'h_e', 'h_p', 'h2_c', 'h2_e', 'h2_p', 'h2o_c', 'h2o_e', 'h2o_p', 'h2o2_c', 'h2o2_e', 'h2o2_p', 'nh4_c', 'nh4_e', 'nh4_p', 'o2_c', 'o2_e', 'o2_p', 'pi_c', 'pi_e', 'pi_p', 'ppi_c', 'pppi_c', 'q8_c', 'q8h2_c', 'no_p', 'no3_p', 'no3_e', 'no_c', 'no2_c', 'no_e', 'no3_c', 'no2_e', 'no2_p', 'n2o_c', 'n2o_e', 'n2o_p', 'h2s_c', 'so3_c', 'so3_p', 'o2s_c', 'h2s_p', 'so2_e', 'so4_e', 'h2s_e', 'o2s_p', 'o2s_e', 'so4_c', 'so4_p', 'so3_e', 'so2_c', 'so2_p', 'ag_c', 'ag_e', 'ag_p','na1_c','na1_e','na1_p','ca2_c', 'ca2_e', 'ca2_p', 'cl_c', 'cl_e', 'cl_p', 'cd2_c', 'cd2_e', 'cd2_p', 'cu_c', 'cu_e', 'cu_p', 'cu2_c', 'cu2_e', 'cu2_p', 'fe2_c', 'fe2_e', 'fe2_p', 'fe3_c', 'fe3_e', 'fe3_p', 'hg2_c', 'hg2_e', 'hg2_p', 'k_c', 'k_e', 'k_p', 'mg2_c', 'mg2_e', 'mg2_p', 'mn2_c', 'mn2_e', 'mn2_p', 'zn2_c', 'zn2_e', 'zn2_p','nh3']
    currencylist2 = ['amp_c', 'amp_e', 'amp_p', 'adp_c', 'adp_e', 'adp_p', 'atp_c', 'atp_e', 'atp_p', 'cmp_c', 'cmp_e', 'cmp_p', 'cdp_c', 'cdp_e', 'cdp_p', 'ctp_c', 'ctp_e', 'ctp_p', 'gmp_c', 'gmp_e', 'gmp_p', 'gdp_c', 'gdp_e', 'gdp_p', 'gtp_c', 'gtp_e', 'gtp_p', 'imp_c', 'imp_e', 'imp_p', 'idp_c', 'idp_e', 'idp_p', 'itp_c', 'itp_e', 'itp_p', 'ump_c', 'ump_e', 'ump_p', 'udp_c', 'udp_e', 'udp_p', 'utp_c', 'utp_e', 'utp_p', 'xmp_e', 'xmp_c', 'xmp_p', 'xdp_c', 'xdp_e', 'xdp_p', 'xtp_c', 'xtp_e', 'xtp_p', 'damp_c', 'damp_e', 'damp_p', 'dadp_c', 'dadp_e', 'dadp_p', 'datp_c', 'datp_e', 'datp_p', 'dcmp_c', 'dcmp_e', 'dcmp_p', 'dcdp_c', 'dcdp_e', 'dcdp_p', 'dctp_c', 'dctp_e', 'dctp_p', 'dgmp_c', 'dgmp_e', 'dgmp_p', 'dgdp_c', 'dgdp_e', 'dgdp_p', 'dgtp_c', 'dgtp_e', 'dgtp_p', 'dimp_c', 'dimp_e', 'dimp_p', 'didp_c', 'didp_e', 'didp_p', 'ditp_c', 'ditp_e', 'ditp_p', 'dump_c', 'dump_e', 'dump_p', 'dudp_c', 'dudp_e', 'dudp_p', 'dutp_c', 'dutp_e', 'dutp_p', 'dtmp_c', 'dtmp_e', 'dtmp_p', 'dtdp_c', 'dtdp_e', 'dtdp_p', 'dttp_c', 'dttp_e', 'dttp_c', 'fad_c', 'fad_p', 'fad_e', 'fadh2_c', 'fadh2_e', 'fadh2_p', 'nad_c', 'nad_e', 'nad_p', 'nadh_c', 'nadh_e', 'nadh_p', 'nadp_c', 'nadp_e', 'nadp_p', 'nadph_c', 'nadph_e', 'nadph_p']
    currencylist3 = ['cdp', 'ag', 'dctp', 'dutp', 'ctp', 'gdp', 'gtp', 'ump', 'ca2', 'h2o', 'datp', 'co2', 'no2', 'no', 'k', 'zn2', 'no3', 'o2', 'cl', 'udp', 'damp', 'ditp', 'dump', 'q8h2', 'pppi', 'idp', 'dimp', 'pi', 'dttp', 'so4', 'adp', 'xtp', 'dgtp', 'dadp', 'coa', 'ppi', 'h2', 'cmp', 'fe2', 'o2s', 'h', 'gmp', 'itp', 'q8', 'cobalt2', 'n2o', 'xmp', 'xdp', 'nadph', 'cu', 'cu2', 'atp', 'dgmp', 'imp', 'h2s', 'utp', 'dtmp', 'fadh2', 'so3', 'fad', 'cd2', 'dgdp', 'nad', 'nadh', 'hg2', 'dcmp', 'dudp', 'dtdp', 'didp', 'mn2', 'dcdp', 'nh4', 'amp', 'fe3', 'nadp', 'so2', 'h2o2', 'mg2']
    return set(currencylist1 + currencylist2 + currencylist3)

def extract_metabolite_id(met_id):
    '''
    Extracts the metabolite ID without compartment information.
    '''
    if '_c' in met_id:
        return met_id.split('_c')[0]
    elif '_p' in met_id:
        return met_id.split('_p')[0]
    elif '_e' in met_id:
        return met_id.split('_e')[0]
    return met_id

def preprocess_metabolites(model, currency_metabolites):
    '''
    Preprocesses the metabolites and creates a DataFrame.

    Arguments
    ----------
    * model: cobra.Model ~ A Model object to extract metabolite information from.
    * currency_metabolites: set ~ A set of currency metabolite IDs.

    Returns
    ----------
    * metdf: pandas.DataFrame ~ DataFrame with preprocessed metabolite information.
    '''
    rlist, sublist, subtotal, gprlist, genetotal = [], [], [], [], []

    for reaction in model.reactions:
        if reaction.gene_reaction_rule and len(reaction.reactants) > 1:
            for metabolite in reaction.reactants:
                if metabolite.id not in currency_metabolites:
                    for gene in reaction.genes:
                        rlist.append(reaction.id)
                        sublist.append(metabolite.id)
                        gprlist.append(reaction.gene_reaction_rule)
                        genetotal.append(gene.id)
                        subtotal.append(extract_metabolite_id(metabolite.id))

    metdf = pd.DataFrame({
        'reactions': rlist,
        'metabolites': sublist,
        'metabolitestotal': subtotal,
        'gpr': gprlist,
        'genes': genetotal,
    })

    return metdf

def split_substrate_to_match_gene(model, metabolites_reactions_gpr_file):
    '''
    Splits the substrate of reactions to match the gene of reactions and saves the results to a file.

    Arguments
    ----------
    * model: cobra.Model ~ A Model object which will be modified in place.
    * metabolites_reactions_gpr_file: str ~ File path to save the resulting DataFrame as a CSV file.

    Returns
    ----------
    * metdf: pandas.DataFrame ~ DataFrame with metabolite information.
    '''
    currency_metabolites = get_currency_metabolites()
    
    # Perform necessary modifications to the model
    convert_to_irreversible(model)
    model = isoenzyme_split(model)

    # Preprocess the metabolites and create the DataFrame
    metdf = preprocess_metabolites(model, currency_metabolites)
    
    # Save the DataFrame to a file
    metdf.to_csv(metabolites_reactions_gpr_file, index=False)
    return metdf

def combine_reactions_simles_sequence(metdf, smilesdf, prodf, comdf_file):
    '''
    This function is used to combine the reaction--substrate--gene--protein_sequnce--mass.
    
    Arguments:
    * metdf: metabolites_reactions_gpr
    * similesdf: inchkeydf
    * prodf: prodf
    '''
    # Create a dictionary for mapping metname to Canonical_SMILES
    metname_to_smiles = dict(zip(smilesdf['metname'], smilesdf['Canonical_SMILES']))

    # Set the index of prodf to 'geneid'
    prodf.set_index('geneid', inplace=True)
    
    for index, row in metdf.iterrows():
        metname = row['metabolitestotal']
        genes = row['genes']
        gpr = row['gpr']

        if metname in metname_to_smiles:
            metdf.loc[index, 'similes'] = metname_to_smiles[metname]
        
        if genes in prodf.index:
            metdf.loc[index, 'prosequence'] = prodf.loc[genes, 'aaseq']
            metdf.loc[index, 'mass'] = prodf.loc[genes, 'mass']
        
        totalmass = 0.0
        genelist = gpr.split(' and ')
        for gene in genelist:
            if gene in prodf.index:
                totalmass += prodf.loc[gene, 'subunitmass']
        
        metdf.loc[index, 'totalmass'] = totalmass
    
    # Reset the index of prodf
    prodf.reset_index(inplace=True)
    
    metdf.to_csv(comdf_file, index=False)
    return metdf

def generate_DLKCAT_input(metdf, metdf_name, metdf_outfile, DLinputdf_file):
    '''
    This function is used to generate the DLKCAT input file.
    
    Arguments:
    * metdf: pandas.DataFrame - DataFrame containing reaction, substrate, gene, protein sequence, and mass information.
    * metdf_name: pandas.DataFrame - DataFrame containing metabolite names.
    * metdf_outfile: str - File path to save the modified metdf DataFrame as a CSV file.
    * DLinputdf_file: str - File path to save the DLKCAT input DataFrame as a CSV file.

    Returns:
    * DLinputdf: pandas.DataFrame - DataFrame containing DLKCAT input data.
    '''
    # Generate the DLKCAT input file
    metdf_name.index = metdf_name['met']
    metdf['metname'] = metdf_name.loc[metdf['metabolitestotal'], 'name'].values

    # Drop rows with missing prosequence
    metdf.dropna(subset=['prosequence'], inplace=True)

    # Replace specific values with None
    metdf['similes'].replace(['noinchikey_inpubchem', 'noinchikey'], 'None', inplace=True)

    # Remove trailing decimals from similes
    metdf['similes'] = metdf['similes'].astype(str).str.split('.').str[0]

    metdf.reset_index(drop=True, inplace=True)
    metdf.to_csv(metdf_outfile, index=False)

    DLinputdf = metdf[['metname', 'similes', 'prosequence']].copy()
    DLinputdf.dropna(subset=['prosequence'], inplace=True)
    DLinputdf.reset_index(drop=True, inplace=True)
    DLinputdf.to_csv(DLinputdf_file, sep='\t', index=False)

    print('DLKCAT input file generated')
    return DLinputdf

def DL_kcat_mw_calculation(DLouputdf, metdf):
    '''
    Make the kcat_mw input file for ecGEM construction.

    Arguments:
    * DLouputdf: pandas.DataFrame - DLinputdf from DLKCAT.
    * metdf: pandas.DataFrame - DataFrame containing metabolite, reaction, gene, similes, prosequence, mass information.
    '''
    DLouputdf.reset_index(drop=True, inplace=True)
    metdf.reset_index(drop=True, inplace=True)

    # Combine DLouputdf and metdf based on index
    DLoutputdf_rex = pd.concat([metdf, DLouputdf['Kcat value (1/s)']], axis=1)

    # Remove rows with 'None' in Kcat value
    DLoutputdf_rex = DLoutputdf_rex[DLoutputdf_rex['Kcat value (1/s)'] != 'None']

    # Convert Kcat value to float
    DLoutputdf_rex['Kcat value (1/s)'] = DLoutputdf_rex['Kcat value (1/s)'].astype(float)

    # Sort by Kcat value and keep only the first occurrence of each reaction
    DLoutputdf_rex = DLoutputdf_rex.sort_values('Kcat value (1/s)', ascending=False).drop_duplicates(subset=['reactions'], keep='first')

    # Calculate kcat_mw
    DLoutputdf_rex['kcat_mw'] = DLoutputdf_rex['Kcat value (1/s)'] * 3600 * 1000 / DLoutputdf_rex['totalmass']

    # Prepare DL_reaction_kact_mw DataFrame
    DL_reaction_kact_mw = pd.DataFrame()
    DL_reaction_kact_mw['reactions'] = DLoutputdf_rex['reactions']
    DL_reaction_kact_mw['data_type'] = 'DLkcat'
    DL_reaction_kact_mw['kcat'] = DLoutputdf_rex['Kcat value (1/s)']
    DL_reaction_kact_mw['MW'] = DLoutputdf_rex['mass']
    DL_reaction_kact_mw['kcat_MW'] = DLoutputdf_rex['kcat_mw']
    DL_reaction_kact_mw.reset_index(drop=True, inplace=True)

    print('DL_reaction_kact_mw generated')
    return DL_reaction_kact_mw

def get_gene_subunitDescription(sub_description_path, model):
    '''
    Get gene subunit descriptions from UniProt for a given model.

    Arguments:
    * sub_description_path: str - File path to save the resulting DataFrame as a CSV file.
    * model: cobra.Model - A Model object.

    Returns:
    None
    '''
    print('Start downloading from UniProt...')
    start = time.time()
    
    gprlist = []
    for eachg in model.genes:
        if 'uniprot' in eachg.annotation.keys():
            gprlist.append(eachg.annotation['uniprot'])
    
    df = pd.DataFrame()
    u = UniProt(verbose=False)
    lALLUniprotIds = []
    lAllProteinNames = []
    lAllGeneNames = []
    lsubunits = []
    
    for el in probar(gprlist):
        res = u.search("%s" % el, frmt="xml")
        time.sleep(3.0)
        geneNames = []
        proteinNames = []
        
        if res != '':
            try:
                dRes = xmltodict.parse(res)
                if 'entry' in dRes['uniprot'].keys():
                    try:
                        for k, val in dRes['uniprot']['entry']['protein']['recommendedName'].items():
                            if k == 'fullName' or k == 'shortName':
                                if type(val) == str:
                                    proteinNames.append(val)
                                else:
                                    proteinNames.append(val['#text'])
                    except:
                        pass

                    try:
                        for diz in dRes['uniprot']['entry']['protein']['alternativeName']:
                            for k, val in diz.items():
                                if k == 'fullName' or k == 'shortName':
                                    if type(val) == str:
                                        proteinNames.append(val)
                                    else:
                                        proteinNames.append(val['#text'])
                    except:
                        pass

                    try:
                        for k, val in dRes['uniprot']['entry']['protein']['submittedName'].items():
                            if k == 'fullName' or k == 'shortName':
                                if type(val) == str:
                                    proteinNames.append(val)
                                else:
                                    proteinNames.append(val['#text'])
                    except:
                        pass

                    try:
                        for diz in dRes['uniprot']['entry']['gene']['name']:
                            for k, val in diz.items():
                                if k == '#text':
                                    geneNames.append(val)
                    except:
                        pass

                    subunitDescription = ''

                    if 'entry' in dRes['uniprot'].keys():
                        if 'comment' in dRes['uniprot']['entry'].keys() and type(dRes['uniprot']['entry']['comment']) == list:
                            try:
                                for diz in dRes['uniprot']['entry']['comment']:
                                    if diz['@type'] == 'subunit':
                                        for k, val in diz.items():
                                            if k == 'text':
                                                subunitDescription = val['#text']
                            except:
                                for ell in dRes['uniprot']['entry']['comment']:
                                    if ell['@type'] == 'subunit':
                                        subunitDescription = ell['text']

                    else:
                        print(el + ' no entry.')
            except:
                pass
            lALLUniprotIds.append(el)
            lAllGeneNames.append(geneNames)
            lsubunits.append(subunitDescription)
            lAllProteinNames.append(proteinNames)
            
    df['UniProtID'] = lALLUniprotIds
    df['proteinNames'] = lAllProteinNames
    df['txt_subunit'] = lsubunits
    df['geneNames'] = lAllGeneNames

    end = time.time()
    run_time = end - start
    
    hour = run_time // 3600
    minute = (run_time - 3600 * hour) // 60
    second = run_time - 3600 * hour - 60 * minute
    
    df = df.fillna(" ")
    df.to_csv(sub_description_path)
    
    print('Success downloading! :-)')

def get_subunit_number(sub_description_path, gene_subnum_path):
    '''
    Get the number of subunits for each gene from the subunit description DataFrame.

    Arguments:
    * sub_description_path: str - File path to the subunit description DataFrame (CSV file).
    * gene_subnum_path: str - File path to save the resulting DataFrame with subunit numbers as a CSV file.

    Returns:
    * df: pd.DataFrame - DataFrame with subunit numbers for each gene.
    '''
    # Load subunit description DataFrame
    df = pd.read_csv(sub_description_path)
    
    # Drop rows with missing values in the 'txt_subunit' column
    df = df.dropna(subset=['txt_subunit'])
    
    # Initialize a list to store subunit numbers
    lSubunitNum = []
    
    for row in df.itertuples():
        lPossibleSubunitsFromName = []
        lAllNames = eval(row.proteinNames) + eval(row.geneNames)
        lnameEnzyme = []
        possibleWords = ['subunit', 'subunits', 'component', 'alpha chain', 'beta chain', 'gamma chain',
                         '30S ribosomal protein', '50S ribosomal protein', 'binding protein', 'large chain',
                         'small chain', 'permease protein', 'insertase', 'translocase protein', 'accessory protein',
                         'UvrABC system protein', 'Chaperonin', 'Co-chaperonin', 'assembly factor', 'recombinase',
                         'flavoprotein']
        namesWithSubunitsIndication = []
        
        for p in possibleWords:
            namesWithSubunitsIndication = [el for el in lAllNames if p in el]
            
            if re.search('RNA-binding protein|methyltransferase|Nucleotide-binding protein|Phosphotransferase|\
                          Transport|AMP-forming|Two-component|SsrA|Ribonuclease', str(namesWithSubunitsIndication)):
                continue
            
            if len(namesWithSubunitsIndication) != 0:
                for n in namesWithSubunitsIndication:
                    if n.endswith(p) or p + ',' in n:
                        try:
                            lPossibleSubunitsFromName.append(n.split(p)[0].split()[-1].strip(','))
                            lnameEnzyme.append(' '.join(n.split(p)[0].split()[:-1]).lower().strip())
                            break
                        except:
                            print(n)
                    elif 'alpha-ketoacid' in namesWithSubunitsIndication:
                        try:
                            lPossibleSubunitsFromName.append(n.split(p)[0].split()[-1].strip(','))
                            lnameEnzyme.append(' '.join(n.split(p)[0].split()[:-1]).lower().strip())
                            break
                        except:
                            print(n)
                    else:
                        try:
                            lPossibleSubunitsFromName.append(n.split(p)[1].split()[0].strip(','))
                            lnameEnzyme.append(n.split(p)[0].lower().strip())
                            break
                        except:
                            print(n)
        
        sub_dict = {
            'Monomer.': '1', 'Homodimer.': '2', 'homodimers': '2', 'homodimer around DNA': '2', 'Homotrimer': '3',
            'Homotetramer.': '4', 'Homopentamer.': '5', 'Homohexamer.': '6', 'Homohexameric': '6', 'Homohexamer': '6',
            '7 subunits': '7', 'ClpP subunits': '7', 'Heterooligomer': '7', 'Homooctamer.': '8', 'Homodecamer.': '10',
            'Homodecamer;': '10', 'homodecamer,': '10', 'dodecamer': '12', 'cylinder': '7', 'Heterotrimer': '1',
            'Associates with': '1', 'Dimer of': '1', 'cytochrome bc1': '1', 'Heterodimer': '1', 'heterodimer': '1',
            'Tat system': '1', 'Heterotetramer': '2', 'heterotetramer': '2', 'Tetramer of': '2', 'Pup': '2',
            'Composed of two chains': '2', 'UreD, UreF and UreG': '2', 'Heterooctamer': '4', '24-polypeptide': '24',
            '24 subunits': '24'
        }
        
        for eachkey in sub_dict.keys():
            sub_num = 'manual'
            
            if re.search(eachkey, row.txt_subunit):
                sub_num = sub_dict[eachkey]
                break
            
            if len(row.txt_subunit) == 1:
                sub_num = 'vacuum'
                break
            
            # Extract subunit number from complex txt_subunit sentences
            if sub_num == 'manual':
                possibleWords = list(set(lPossibleSubunitsFromName))
                
                if len(possibleWords) > 0:
                    for p in possibleWords:
                        namesWithSubunitsIndication = []
                        
                        try:
                            namesWithSubunitsIndication = row.txt_subunit[0]
                        except:
                            namesWithSubunitsIndication = row.txt_subunit
                        
                        if re.search(p, str(namesWithSubunitsIndication)):
                            if re.search('FGAM|RNAP', namesWithSubunitsIndication):
                                try:
                                    sub_num = namesWithSubunitsIndication.split(p)[0].split()[-1].strip(',')
                                except:
                                    sub_num == 'manual'
                            elif 'CF' in namesWithSubunitsIndication:
                                sub_num = namesWithSubunitsIndication.split(p)[1].split()[0].strip(',')
                                if p == 'a':
                                    sub_num = '1'
                                    break
                            else:
                                sub_num = '1'
                            
                            if re.search('UreD', p):
                                sub_num = '1'
                        else:
                            sub_num = '1'
                else:
                    sub_num = '1'
        
        namesWithSubunitsIndication = row.txt_subunit
        
        if re.search('ATP-binding proteins', namesWithSubunitsIndication):
            sub_num = '2'
        if re.search('transmembrane proteins', namesWithSubunitsIndication):
            sub_num = '2'
        if re.search('solute-binding protein', namesWithSubunitsIndication):
            sub_num = '1'
        if re.search('ribosom', namesWithSubunitsIndication):
            sub_num = '1'
        
        sub_num = sub_num.strip('.').strip('()')
        lSubunitNum.append(sub_num)
    
    # Add subunit numbers to the DataFrame
    df['subunitnumber'] = lSubunitNum
    df.set_index('UniProtID', inplace=True)
    
    # Save the resulting DataFrame to a CSV file
    df.to_csv(gene_subnum_path)
    
    return df

def get_PhPP_data(model_file, model_type, obj, substrate_bound, o2_bound, number, outputfile, x_id, z_id):
    '''
    Perform PhPP analysis on a model using different substrate and oxygen bounds.

    Arguments:
    * model_file: str - File path to the model file (JSON or SBML format).
    * model_type: str - Type of the model ('GEM' for Genome-Scale Model, other types for enzyme-constraint models).
    * obj: str - ID of the objective reaction in the model.
    * substrate_bound: float - Upper bound for the substrate constraint.
    * o2_bound: float - Upper bound for the oxygen constraint.
    * number: int - Number of intervals to divide the bounds.
    * outputfile: str - File path to save the resulting DataFrame as a CSV file.
    * x_id: str - ID of the first reaction to constrain.
    * z_id: str - ID of the second reaction to constrain.

    Returns:
    * df1: pd.DataFrame - DataFrame containing the results of PhPP analysis.
    '''

    # Generate equally spaced lists
    exlistn = np.linspace(1, substrate_bound, number+1).round(4).tolist()
    exlistmo2 = np.linspace(1, o2_bound, number+1).round(4).tolist()

    df1 = pd.DataFrame(index=range(number+1), columns=range(number+1))

    for v, i in enumerate(exlistmo2[1:], start=1):
        condi = i
        for k, j in enumerate(exlistn[1:], start=1):
            condj = j

            # Load the model
            if model_type == 'GEM':
                model = cobra.io.json.load_json_model(model_file)
            else:
                model = get_enzyme_constraint_model(model_file)

            try:
                model.reactions.get_by_id(z_id).bounds = (0, 0)
            except:
                print(z_id + ' not in the model')
            else:
                try:
                    model.reactions.get_by_id(z_id+'_reverse').bounds = (0, condi)
                except:
                    model.reactions.get_by_id(z_id).bounds = (0, condi)

            try:
                model.reactions.get_by_id(x_id).bounds = (0, 0)
            except:
                print(x_id + ' not in the model')
            else:
                try:
                    model.reactions.get_by_id(x_id+'_reverse').bounds = (0, condj)
                except:
                    model.reactions.get_by_id(x_id).bounds = (0, condj)

            model.objective = obj
            enz_model_pfba_solution = cobra.flux_analysis.pfba(model)
            df1.iloc[k, v] = enz_model_pfba_solution.fluxes[obj]

    df1.to_csv(outputfile, index=False)
    return df1

def drawphpp(GEM_glc_o2_df, ecGEM_glc_o2_df, substrate_bound, o2_bound, obj_bound, number, PhPP_output_fig_file):
    '''
    Generate a 3D surface plot of PhPP results for GEM and ecGEM.

    Arguments:
    * GEM_glc_o2_df: pd.DataFrame - DataFrame containing GEM growth rates for different glucose and oxygen uptake rates.
    * ecGEM_glc_o2_df: pd.DataFrame - DataFrame containing ecGEM growth rates for different glucose and oxygen uptake rates.
    * substrate_bound: float - Upper bound for the substrate uptake rates.
    * o2_bound: float - Upper bound for the oxygen uptake rates.
    * obj_bound: float - Upper bound for the growth rates.
    * number: int - Number of intervals to divide the bounds.
    * PhPP_output_fig_file: str - File path to save the generated figure.

    Returns:
    * fig: go.Figure - Generated 3D surface plot figure.
    '''

    fig = make_subplots(rows=1, cols=2, column_widths=[0.5, 0.5], specs=[[{"type": "surface"}, {"type": "surface"}]])

    x_values = np.linspace(1, o2_bound, number)
    y_values = np.linspace(1, substrate_bound, number)

    fig.add_trace(
        go.Surface(
            y=y_values,
            x=x_values,
            z=GEM_glc_o2_df.values,
            coloraxis="coloraxis"
        ),
        row=1,
        col=1
    )

    fig.add_trace(
        go.Surface(
            y=y_values,
            x=x_values,
            z=ecGEM_glc_o2_df.values,
            coloraxis="coloraxis"
        ),
        row=1,
        col=2
    )

    fig.update_layout(
        scene=dict(
            xaxis=dict(
                range=[0, o2_bound],
                tickfont=dict(size=13, family='Times New Roman'),
                backgroundcolor="lightgrey",
                title=dict(text="<b>O2 uptake rates<br>(mmol/gDW/h)</b>", font=dict(size=15, family='Times New Roman'))
            ),
            yaxis=dict(
                range=[0, o2_bound],
                tickfont=dict(size=13, family='Times New Roman'),
                backgroundcolor="lightgrey",
                title=dict(text="<b>Glucose uptake rates<br>(mmol/gDW/h)</b>", font=dict(size=15, family='Times New Roman'))
            ),
            zaxis=dict(
                range=[0, obj_bound],
                tickfont=dict(size=13, family='Times New Roman'),
                backgroundcolor="grey",
                gridcolor="white",
                title=dict(text="<b>GEM Growth rates (1/h)</b>", font=dict(size=15, family='Times New Roman'))
            )
        ),
        scene2=dict(
            xaxis=dict(
                range=[0, o2_bound],
                tickfont=dict(size=13, family='Times New Roman'),
                backgroundcolor="lightgrey",
                title=dict(text="<b>O2 uptake rates<br>(mmol/gDW/h)</b>", font=dict(size=15, family='Times New Roman'))
            ),
            yaxis=dict(
                range=[0, o2_bound],
                tickfont=dict(size=13, family='Times New Roman'),
                backgroundcolor="lightgrey",
                title=dict(text="<b>Glucose uptake rates<br>(mmol/gDW/h)</b>", font=dict(size=15, family='Times New Roman'))
            ),
            zaxis=dict(
                range=[0, obj_bound],
                tickfont=dict(size=13, family='Times New Roman'),
                backgroundcolor="grey",
                gridcolor="white",
                title=dict(text="<b>ecGEM Growth rates (1/h)</b>", font=dict(size=15, family='Times New Roman'))
            )
        ),
        autosize=False,
        width=1150,
        height=550,
        margin=dict(l=1, r=3, b=10, t=20)
    )

    fig.update_layout(
        coloraxis={'colorscale': [[0, 'green'], [1, 'red']]}
    )

    fig.write_image(PhPP_output_fig_file)

    return fig

def draw_cdf_fig(data_cdf_data, output_file, x_name, y_name, y_index, nticks):
    '''
    Generate a cumulative distribution function (CDF) plot.

    Arguments:
    * data_cdf_data: list - Data for the CDF plot.
    * output_file: str - File path to save the generated figure.
    * x_name: str - Name of the x-axis.
    * y_name: str - Name of the y-axis.
    * y_index: list - Index values for the y-axis.
    * nticks: int - Number of ticks for the x-axis.

    Returns:
    * fig: go.Figure - Generated CDF plot figure.
    '''

    trace0 = go.Scatter(
        x=data_cdf_data,
        y=y_index,
        mode='lines',
        marker={'color': 'blue'}
    )
    trace1 = go.Scatter(
        x=data_cdf_data,
        y=y_index,
        mode='lines',
        marker={'color': 'blue'},
        line={'color': 'blue', 'width': 3},
        xaxis='x2',
        yaxis="y2"
    )
    data1 = [trace0, trace1]

    layout = go.Layout(
        plot_bgcolor='lightgrey',
        xaxis=dict(
            title=dict(
                text=x_name,
                font=dict(size=20, family='Times New Roman')
            ),
            type="log",
            rangemode="tozero",
            tickfont=dict(color='black', size=20, family='Times New Roman'),
            linecolor='black',
            ticks='inside',
            tickcolor='black',
            zeroline=False,
            showexponent='all',
            exponentformat="power",
            gridcolor="yellow"
        ),
        xaxis2=dict(
            linecolor='black',
            showticklabels=False,
            type="log",
            tickfont=dict(color='black', size=20, family='Times New Roman'),
            rangemode="tozero",
            overlaying='x',
            side='top',
            nticks=nticks,
            zeroline=False,
            showexponent='all',
            exponentformat="power",
            gridcolor="white"
        ),
        yaxis=dict(
            title=dict(
                text=y_name,
                font=dict(size=20, family='Times New Roman')
            ),
            range=[0, 1],
            showgrid=False,
            zeroline=False,
            rangemode="tozero",
            tickfont=dict(color='black', size=20, family='Times New Roman'),
            ticks='inside',
            tickcolor='black',
            linecolor='black'
        ),
        yaxis2=dict(
            range=[0, 1],
            linecolor='black',
            showgrid=False,
            zeroline=False,
            tickfont=dict(color='black', size=20, family='Times New Roman'),
            showticklabels=False,
            overlaying='y',
            side='right'
        ),
        showlegend=False,
        height=450,
        width=750,
        margin=go.layout.Margin(l=10, r=10, b=10, t=10)
    )

    fig = go.Figure(data1, layout=layout)
    fig.add_hline(y=0.5, line_width=2, line_color="orange")
    fig.write_image(output_file)

    return fig

def draw_3d_rbas(z_data, substrate_bound, o2_bound, obj_bound, number, PhPP_output_fig_file):
    '''
    Generate a 3D plot of reaction bounds analysis (RBA) results.

    Arguments:
    * z_data: pd.DataFrame - Dataframe containing the z-axis values for the plot.
    * substrate_bound: float - Upper bound of the substrate.
    * o2_bound: float - Upper bound of O2 uptake rates.
    * obj_bound: float - Upper bound of growth rates.
    * number: int - Number of data points.
    * PhPP_output_fig_file: str - File path to save the generated figure.

    Returns:
    * fig: go.Figure - Generated 3D plot figure.
    '''

    x_values = np.linspace(1, o2_bound, number)
    y_values = np.linspace(1, substrate_bound, number)

    layout = go.Layout(
        template="none",
        plot_bgcolor='lightgrey',
        scene=dict(
            xaxis=dict(
                range=[0, o2_bound],
                tickfont=dict(size=13, family='Times New Roman'),
                backgroundcolor="lightgrey",
                title=dict(text="<b>O2 uptake rates<br>(mmol/gDW/h)</b>",
                           font=dict(size=18, family='Times New Roman'))
            ),
            yaxis=dict(
                range=[0, substrate_bound],
                tickfont=dict(size=13, family='Times New Roman'),
                backgroundcolor="lightgrey",
                title=dict(text="<b>Glucose uptake rates<br>(mmol/gDW/h)</b>",
                           font=dict(size=18, family='Times New Roman'))
            ),
            zaxis=dict(
                range=[0, obj_bound],
                tickfont=dict(size=13, family='Times New Roman'),
                backgroundcolor="grey",
                gridcolor="white",
                title=dict(text="<b>Growth rates<br>(1/h)</b>",
                           font=dict(size=18, family='Times New Roman'))
            )
        ),
        autosize=False,
        scene_camera_eye=dict(x=-0.8, y=-2.1, z=0.3),
        width=850,
        height=850,
        margin=dict(l=20, r=20, b=20, t=20)
    )

    fig = go.Figure(data=[go.Surface(x=x_values, y=y_values, z=z_data.values)], layout=layout)

    fig.update_traces(contours_z=dict(usecolormap=True, highlightcolor="mistyrose", project_z=True))
    fig.update_scenes(yaxis_tickangle=0, xaxis_tickangle=0)

    fig.write_image(PhPP_output_fig_file)
    return fig

def calculate_yield(model, use_substrate, substrate_name,glc_concentration,columns, correspond_rxn,yield_list):
    '''
    Calculate yield using the provided model and glucose concentration.

    Arguments:
    * model: cobra.Model - COBRApy model object.
    * use_substrate: str - substrate ID in model.
    * substrate_name: str - substrate name.
    * glc_concentration: float - substrate concentration.
    * columns: list - List of reactions corresponding to the columns in yield_list
    * yield_list: pd.DataFrame - DataFrame to store the calculated yields.

    Returns:
    * yield_list: pd.DataFrame - DataFrame with updated yield values.
    '''
    try:
        model.reactions.get_by_id(use_substrate).bounds = (-glc_concentration, 0)
    except:
        print(use_substrate + ' not in the model')
    else:
        try:
            model.reactions.get_by_id(use_substrate+'_reverse').bounds = (0, 0)
        except:
            pass

    pfba_solution = cobra.flux_analysis.pfba(model)
    for column in range(len(columns)):
        yield_list.loc[glc_concentration, columns[column]] = pfba_solution.fluxes[correspond_rxn[column]]
    yield_list.loc[glc_concentration, substrate_name] = glc_concentration
    return yield_list

def create_trace(x, y, name, color, size=10, symbol=None, width=None, xaxis=None, yaxis=None):
    '''
    Create a scatter trace for a plot.

    Arguments:
    * x: list or array-like - x-axis data.
    * y: list or array-like - y-axis data.
    * name: str - Name of the trace.
    * color: str - Color of the markers and line.
    * size: int, optional - Size of the markers (default: 10).
    * symbol: str, optional - Symbol for the markers (default: None).
    * width: int or None, optional - Width of the line (default: None).
    * xaxis: str or None, optional - Name of the x-axis (default: None).
    * yaxis: str or None, optional - Name of the y-axis (default: None).

    Returns:
    * trace: go.Scatter - Scatter trace object.
    '''

    # Set marker properties
    marker = {'color': color, 'size': size}

    # Set line properties if width is provided
    line = {'color': color, 'width': width} if width else None

    # Create Scatter trace object
    trace = go.Scatter(
        x=x,
        y=y,
        mode='lines+markers',
        name=name,
        marker=marker,
        line=line,
        xaxis=xaxis,
        yaxis=yaxis
    )

    return trace

def draw_overfolw_fig(GEMyield_list, ecGEMyield_list, column_list, y_axis_loc_list, color_list, substrate_name, substrate_bound, obj_bound, secrate_bound, trade_off_biomass_yield_figfile):
    '''
    Create a figure showing trade-off between biomass yield and substrate/secretion rates.

    Arguments:
    * GEMyield_list: pandas DataFrame - Data for GEM yields.
    * ecGEMyield_list: pandas DataFrame - Data for ecGEM yields.
    * column_list: list - List of column names.
    * y_axis_loc_list: list - List of y-axis locations ('left' or 'right') for each column.
    * color_list: list - List of colors for the traces.
    * substrate_bound: float - Upper limit for substrate uptake rate.
    * obj_bound: float - Upper limit for growth rate.
    * secrate_bound: float - Upper limit for secretion rate.
    * trade_off_biomass_yield_figfile: str - Output file name for the figure.

    Returns:
    * fig: go.Figure - Plotly figure object.
    '''

    data1 = []

    for column in range(len(column_list)):
        if y_axis_loc_list[column] == 'left':
            trace1 = create_trace(GEMyield_list[substrate_name], GEMyield_list[column_list[column]], '%s_GEM' % column_list[column], color_list[column], symbol=5, size=10, width=3)
            data1.append(trace1)
            trace2 = create_trace(ecGEMyield_list[substrate_name], ecGEMyield_list[column_list[column]], '%s_ecGEM' % column_list[column], color_list[len(column_list) + column], symbol=5, size=10, width=3)
            data1.append(trace2)
        else:
            trace1 = create_trace(GEMyield_list[substrate_name], GEMyield_list[column_list[column]], '%s_GEM' % column_list[column], color_list[column], size=10, symbol=4, width=3, xaxis='x2', yaxis='y2')
            data1.append(trace1)
            trace2 = create_trace(ecGEMyield_list[substrate_name], ecGEMyield_list[column_list[column]], '%s_ecGEM' % column_list[column], color_list[len(column_list) + column], size=10, symbol=4, width=3, xaxis='x2', yaxis='y2')
            data1.append(trace2)

    layout = go.Layout(
        plot_bgcolor='white',
        xaxis=dict(
            title=dict(text="<b>Substrate uptake rate (mmol/gDW/h)</b>",
                       font=dict(size=20, family='Times New Roman')),
            tickfont=dict(color='black', size=15, family='Times New Roman'),
            range=[1, substrate_bound],
            linecolor='black',
            ticks='inside',
            tickcolor='black'
        ),
        xaxis2=dict(
            linecolor='black',
            showticklabels=False,
            overlaying='x',
            side='top',
            range=[1, substrate_bound]
        ),
        yaxis=dict(
            title=dict(text="<b>Growth rate (h<sup>-1</sup>)</b>",
                       font=dict(size=20, family='Times New Roman')),
            tickfont=dict(color='black', size=15, family='Times New Roman'),
            ticks='inside',
            tickcolor='black',
            linecolor='black',
            range=[0, obj_bound]
        ),
        yaxis2=dict(
            title=dict(text="<b>Secrete rate (mmol/gDW/h)</b>",
                       font=dict(size=20, family='Times New Roman')),
            overlaying='y',
            tickfont=dict(color='black', size=15, family='Times New Roman'),
            ticks='inside',
            tickcolor='black',
            linecolor='black',
            side="right",
            range=[0, secrate_bound]
        ),
        legend=dict(x=0.05, y=0.98, font=dict(size=20, color="black", family='Times New Roman')),
        width=800,
        height=600
    )

    fig = go.Figure(data=data1, layout=layout)
    fig.write_image(trade_off_biomass_yield_figfile)

    return fig

def generate_random_colors(num_colors):
    '''
    Generate a list of random colors in hexadecimal format.

    Arguments:
    * num_colors: int - Number of colors to generate.

    Returns:
    * colors: list - List of random colors in hexadecimal format.
    '''

    colors = []

    for _ in range(num_colors):
        # Generate random RGB values
        red = random.randint(0, 255)
        green = random.randint(0, 255)
        blue = random.randint(0, 255)

        # Convert RGB values to hexadecimal color code
        color = "#{:02x}{:02x}{:02x}".format(red, green, blue)
        colors.append(color)

    return colors

def get_min_enzyme_cost(model, dict_coeff):
    """Get model flux using Minimum enzyme cost algorithm

    Arguments
    ----------
    * model: cobra model.
    * dict_coeff: {reaction ID: coeffient}.
    
    :return: cobra solution.
    """
    with model:
        bounds = (model.slim_optimize(), model.slim_optimize())
        cons_obj = model.problem.Constraint(
            model.objective.expression,
            lb=min(bounds), ub=max(bounds))
        model.add_cons_vars(cons_obj)

        dict_obj = dict()
        for r in model.reactions:
            if r.id in list(dict_coeff.index):
                #print(dict_coeff.loc[r.id,'kcat_MW'])
                dict_obj[r.forward_variable] = 1 / dict_coeff.loc[r.id,'kcat_MW']

        model_obj = model.problem.Objective(Zero, direction="min", sloppy=True)
        model.objective = model_obj
        model.objective.set_linear_coefficients(dict_obj)

        solution = model.optimize()
    return solution

def get_yield_cost_efficiency(enz_model, glc_concentration_list, use_substrate, obj, reaction_kcat_MW, efficiency_file):
    '''
    Calculate yield, cost, and efficiency for different glucose concentrations.

    Arguments:
    * enz_model: cobra.Model - Enzyme-constrained metabolic model.
    * glc_concentration_list: list - List of glucose concentrations to evaluate.
    * use_substrate: str - ID of the substrate reaction in the model.
    * obj: str - ID of the objective reaction in the model.
    * reaction_kcat_MW: dict - Dictionary mapping reaction IDs to kcat/MW values.
    * efficiency_file: str - File path to save the efficiency data.

    Returns:
    * yield_cost_efficiency_df: pd.DataFrame - DataFrame containing yield, cost, and efficiency values.
    '''

    yield_cost_efficiency_df = pd.DataFrame()

    for glc_concentration in glc_concentration_list:
        with enz_model as growth_model:
            # Set the bounds and optimize the model
            try:
                growth_model.reactions.get_by_id('%s_reverse' % use_substrate).bounds = (0.0, 0.0)
            except:
                print(use_substrate + ' not in the model')
            else:
                pass

            growth_model.reactions.get_by_id('ATPM').bounds = (0.0, 1000)
            growth_model.reactions.get_by_id(use_substrate).bounds = (-glc_concentration, 0.0)
            pfba_solution = cobra.flux_analysis.pfba(growth_model)

            # Calculate enzyme cost
            enz_solution = get_min_enzyme_cost(growth_model, reaction_kcat_MW)

            # Store the results in the DataFrame
            yield_cost_efficiency_df.loc[glc_concentration, 'biomass'] = pfba_solution.fluxes[obj]
            yield_cost_efficiency_df.loc[glc_concentration, 'glucose_simu'] = -pfba_solution.fluxes[use_substrate]
            yield_cost_efficiency_df.loc[glc_concentration, 'glucose_set'] = glc_concentration
            yield_cost_efficiency_df.loc[glc_concentration, 'biomass yield'] = -pfba_solution.fluxes[obj] / (
                    pfba_solution.fluxes[use_substrate] * 0.18)
            yield_cost_efficiency_df.loc[glc_concentration, 'min enzyme cost'] = enz_solution.objective_value
            yield_cost_efficiency_df.loc[glc_concentration, 'enzyme efficiency'] = pfba_solution.fluxes[obj] / enz_solution.objective_value

    yield_cost_efficiency_df.to_csv(efficiency_file)
    return yield_cost_efficiency_df

def draw_trade_off(efficiency_pfba, trade_off_enzyme_efficiency_figfile):
    '''
    Draw a trade-off plot showing the relationship between biomass yield and enzyme efficiency.

    Arguments:
    * efficiency_pfba: pd.DataFrame - DataFrame containing efficiency data.
    * trade_off_enzyme_efficiency_figfile: str - File path to save the trade-off plot.

    Returns:
    * fig: go.Figure - Plotly Figure object.
    '''

    # Create traces
    trace0 = go.Scatter(
        x=efficiency_pfba['glucose_set'],
        y=efficiency_pfba['biomass yield'],
        mode='lines',
        line={'dash': 'dash', 'color': 'red', 'width': 3},
        name='biomass yield',
        marker={'color': 'red', 'size': 10}
    )
    trace1 = go.Scatter(
        x=efficiency_pfba['glucose_set'],
        y=efficiency_pfba['enzyme efficiency'],
        mode='lines',
        name='enzyme efficiency',
        marker={'color': 'blue', 'symbol': 5, 'size': 10},
        line={'color': 'blue', 'width': 3},
        xaxis='x2',
        yaxis='y2'
    )
    data1 = [trace0, trace1]

    layout = go.Layout(
        plot_bgcolor='white',
        xaxis=dict(
            title=dict(text="<b>Substrate uptake rate (mmol/gDW/h)<b>", font=dict(size=20, family='Times New Roman')),
            tickfont=dict(color='black', size=15, family='Times New Roman'),
            linecolor='black',
            ticks='inside',
            tickcolor='black',
            range=[1, np.max(efficiency_pfba['glucose_simu']) + 0.1]
        ),
        xaxis2=dict(
            linecolor='black',
            showticklabels=False,
            overlaying='x',
            side='top',
            range=[1, np.max(efficiency_pfba['glucose_simu']) + 0.1]
        ),
        yaxis=dict(
            title=dict(text="<b>Biomass yield (gDW/g glucose)<b>", font=dict(size=20, family='Times New Roman')),
            range=[np.min(efficiency_pfba['biomass yield']) - 0.1, np.max(efficiency_pfba['biomass yield']) + 0.1],
            tickfont=dict(color='black', size=15, family='Times New Roman'),
            ticks='inside',
            tickcolor='black',
            linecolor='black'
        ),
        yaxis2=dict(
            title=dict(text="<b>Enzyme efficiency (gDW/g enzyme)<b>", font=dict(size=20, family='Times New Roman')),
            range=[np.min(efficiency_pfba['enzyme efficiency']) - 0.1, np.max(efficiency_pfba['enzyme efficiency']) + 0.1],
            overlaying='y',
            tickfont=dict(color='black', size=15, family='Times New Roman'),
            ticks='inside',
            tickcolor='black',
            linecolor='black',
            side="right"
        ),
        legend=dict(x=0.07, y=0.95, font=dict(size=20, color="black", family='Times New Roman')),
        width=800,
        height=600
    )

    fig = go.Figure(data=data1, layout=layout)
    fig.write_image(trade_off_enzyme_efficiency_figfile)

    return fig

def get_enz_foldchange(ecModel_file, obj, substrate, substrate_con, biomass_id, biomass_min, biomass_max, FC_threshold,fluxes_outfile, enzcost_diff_file):
    '''
    Calculate the enzyme fold change between the minimum and maximum biomass conditions.

    Arguments:
    * ecModel_file: str - File path to the enzyme-constrained model.
    * obj: str - ID of the objective reaction in the model.
    * substrate: str - ID of the substrate reaction.
    * substrate_con: float - Substrate concentration.
    * biomass_id: str - ID of the biomass reaction.
    * biomass_min: float - Minimum biomass value.
    * biomass_max: float - Maximum biomass value.
    * fluxes_outfile: str - File path to save the fluxes data.
    * enzcost_diff_file: str - File path to save the enzyme cost difference.

    Returns:
    * enzcost_diff: pd.DataFrame - DataFrame containing enzyme cost difference data.
    '''

    enz_model = get_enzyme_constraint_model(ecModel_file)
    enzcost_diff = pd.DataFrame()

    with enz_model as tmp_model:
        tmp_model.reactions.get_by_id(substrate).bounds = (-substrate_con, 0)
        try:
            tmp_model.reactions.get_by_id('%s_reverse' % substrate).bounds = (0, 0)
        except:
            print(substrate + ' not in the model')
        else:
            pass

        tmp_model.reactions.get_by_id(biomass_id).bounds = (biomass_min, biomass_min)
        tmp_model.objective = obj
        enz_model_pfba_solution = cobra.flux_analysis.pfba(tmp_model)
        enzcost_min = get_fluxes_detail_in_model(enz_model, enz_model_pfba_solution, fluxes_outfile, ecModel_file)

    with enz_model as tmp_model:
        tmp_model.reactions.get_by_id(substrate).bounds = (-substrate_con, 0)
        try:
            tmp_model.reactions.get_by_id('%s_reverse' % substrate).bounds = (0, 0)
        except:
            print(substrate + ' not in the model')
        else:
            pass

        tmp_model.reactions.get_by_id(biomass_id).bounds = (biomass_max, biomass_max)
        tmp_model.objective = obj
        enz_model_pfba_solution = cobra.flux_analysis.pfba(tmp_model)
        enzcost_max = get_fluxes_detail_in_model(enz_model, enz_model_pfba_solution, fluxes_outfile, ecModel_file)

    enzcost_diff['E_μ%.2g' % biomass_min] = enzcost_min['E']
    enzcost_diff['E_μ%.2g' % biomass_max] = enzcost_max['E']

    non_zero_min = enzcost_min['E'] != 0
    non_zero_max = enzcost_max['E'] != 0
    valid_values = non_zero_min & non_zero_max

    enzcost_diff['log2_foldchange(max/min)'] = np.nan
    enzcost_diff['log2_foldchange(min/max)'] = np.nan
    enzcost_diff.loc[valid_values, 'log2_foldchange(max/min)'] = np.log2(enzcost_max.loc[valid_values, 'E'] / enzcost_min.loc[valid_values, 'E'])
    enzcost_diff.loc[valid_values, 'log2_foldchange(min/max)'] = np.log2(enzcost_min.loc[valid_values, 'E'] / enzcost_max.loc[valid_values, 'E'])

    weaken_target = enzcost_diff['log2_foldchange(max/min)'] > FC_threshold
    enhance_target = enzcost_diff['log2_foldchange(min/max)'] > FC_threshold
    enzcost_diff['type'] = np.select([enhance_target, weaken_target], ['enhance_target', 'weaken_target'], 'normal')

    enzcost_diff['foldchange(max/min)'] = enzcost_max['E'] / enzcost_min['E']
    enzcost_diff['equ'] = enzcost_min['equ']
    enzcost_diff['ec-code'] = enzcost_min['ec-code']

    filtered_df = enzcost_diff[(enzcost_diff['type'] == 'enhance_target') | (enzcost_diff['type'] == 'weaken_target')]
    filtered_df = filtered_df.sort_values('type', ascending=False)

    filtered_df.to_csv(enzcost_diff_file)
    return filtered_df

def run_FSEOF(model, substrate, substrate_con, biomass_id, obj,FSEOF_file):
    '''
    Perform Flux Scanning-based Enzyme-Feasible Flux Optimization (FSEOF) analysis.

    Arguments:
    * model: cobra.Model - COBRA model object.
    * substrate: str - ID of the substrate reaction.
    * substrate_con: float - Substrate concentration.
    * biomass_id: str - ID of the biomass reaction.
    * obj: str - ID of the objective reaction.

    Returns:
    * FSEOFdf_f: pd.DataFrame - DataFrame containing FSEOF analysis results.
    * FESEOF_gene: pd.DataFrame - DataFrame containing genes and their regulations.
    '''

    model.reactions.get_by_id(substrate).bounds = (-substrate_con, 0)
    try:
        model.reactions.get_by_id(f"{substrate}_reverse").bounds = (0, 0)
    except:
        print(substrate + ' not in the model')
    else:
        pass

    model.objective = biomass_id
    model_pfba_solution = cobra.flux_analysis.pfba(model)
    WTbiomass = model_pfba_solution.fluxes[biomass_id]
    model_solution_frame_wt = model_pfba_solution.to_frame()

    exlist = np.linspace(0.5 * WTbiomass, WTbiomass * 0.9, 10)
    exlistn = [format(i, '.2f') for i in exlist]
    FSEOFdf = pd.DataFrame()

    for cond in exlistn:
        cond = float(cond)
        with model as tmp_model:
            tmp_model.reactions.get_by_id(substrate).bounds = (-substrate_con, 0)
            try:
                tmp_model.reactions.get_by_id(f"{substrate}_reverse").bounds = (0, 0)
            except:
                print(substrate + ' not in the model')
            else:
                pass

            tmp_model.reactions.get_by_id(biomass_id).bounds = (cond, cond)
            tmp_model.objective = obj
            model_pfba_solution = cobra.flux_analysis.pfba(tmp_model)
        model_solution_frame = model_pfba_solution.to_frame()
        FSEOFdf['id'] = model_solution_frame.index
        FSEOFdf.set_index('id', inplace=True)
        FSEOFdf[f"growth = {cond}"] = model_solution_frame['fluxes']
        FSEOFdf[f"FC growth = {cond}"] = FSEOFdf[f"growth = {cond}"] / model_solution_frame_wt['fluxes']

    FSEOFdf['WT'] = model_solution_frame_wt['fluxes']
    FSEOFdf['GPR'] = FSEOFdf.apply(lambda row: model.reactions.get_by_id(row.name).gene_reaction_rule, axis=1)

    FSEOFdf = FSEOFdf.applymap(lambda x: np.where(x != '', x, None))
    FSEOFdf_drop = FSEOFdf.dropna(subset=['GPR'])
    FSEOFdf_drop = FSEOFdf_drop[FSEOFdf_drop.iloc[:, 3:14].notnull().any(axis=1)]
    FSEOFdf_drop.iloc[:, 3:14] = FSEOFdf_drop.iloc[:, 3:14].fillna(1)
    FSEOFdf_drop.replace([np.inf, -np.inf], 1000, inplace=True)
    FSEOFdf_drop['FC_mean'] = FSEOFdf_drop.iloc[:, 3:14].mean(axis=1)
    FSEOFdf_drop.sort_values(by='FC_mean', ascending=False, inplace=True)

    unchanged = FSEOFdf_drop[(FSEOFdf_drop.iloc[:, 3:13] <= 1).all(axis=1) & (FSEOFdf_drop.iloc[:, 3:13] >= 0.95).all(axis=1)]
    unchanged['regulation'] = 'unchanged'

    always_down = FSEOFdf_drop[(FSEOFdf_drop.iloc[:, 3:13] < 0.95).all(axis=1)]
    always_down['regulation'] = 'down'

    always_up = FSEOFdf_drop[(FSEOFdf_drop.iloc[:, 3:13] > 1).all(axis=1)]
    always_up['regulation'] = 'up'

    FSEOFdf_done = pd.concat([always_up, unchanged, always_down])
    FSEOFdf_done.loc[FSEOFdf_done['FC_mean'] == 0, 'regulation'] = 'knockout'

    reaction_ids = FSEOFdf_done.index.tolist()
    unique_reaction_ids = list(set(reaction_ids))

    reactions_equ = {}
    for r_id in unique_reaction_ids:
        reactions_equ[r_id] = model.reactions.get_by_id(r_id).reaction

    FSEOFdf_done['reactions'] = FSEOFdf_done.index.map(reactions_equ)   

    FSEOFdf_done.to_csv(FSEOF_file)

    return FSEOFdf_done

def Determine_suitable_ecGEM(model_file, bigg_met_file):
    '''
    Determine if the provided model is suitable for constructing an enzyme-constrained model.

    Arguments:
    * model_file: str - File path of the model in SBML or JSON format.
    * bigg_met_file: str - File path of the Bigg metabolite data file in tab-separated format.

    Returns:
    * result: str or list - If the model is suitable, returns "Suitable for constructing enzyme-bound models."
                           If there are errors, returns a list of error messages.
    '''

    error_list = []

    # Read the model from the file
    if re.search('\.xml', model_file):
        model = cobra.io.read_sbml_model(model_file)
    elif re.search('\.json', model_file):
        model = cobra.io.json.load_json_model(model_file)

    # Load the Bigg metabolite data
    bigg_met_df = pd.read_csv(bigg_met_file, sep='\t')

    # Check metabolite coverage
    model_met_list = [met.id for met in model.metabolites]
    model_met_in_bigg = [met.id for met in model.metabolites if met.id in bigg_met_df['bigg_id'].tolist()]
    met_coverage = len(model_met_in_bigg) / len(model_met_list)
    if met_coverage < 0.33:
        met_error = f"The coverage of metabolites is too low ({met_coverage*100:.1f}%), and it is not recommended to construct an enzyme-constrained model."
        error_list.append(met_error)

    # Check gene coverage
    model_gene_list = [gene.id for gene in model.genes]
    model_gene_in_uniprot = [gene.id for gene in model.genes if 'uniprot' in gene.annotation]
    gene_coverage = len(model_gene_in_uniprot) / len(model_gene_list)
    if gene_coverage < 0.33:
        gene_error = f"The coverage of genes is too low ({gene_coverage*100:.1f}%), and it is not recommended to construct an enzyme-constrained model."
        error_list.append(gene_error)

    # Check reaction coverage
    model_reaction_list = [reaction.id for reaction in model.reactions if not reaction.id.startswith('EX_')]
    model_reaction_with_EC = [reaction.id for reaction in model.reactions if 'ec-code' in reaction.annotation]
    reaction_coverage = len(model_reaction_with_EC) / len(model_reaction_list)
    if reaction_coverage < 0.33:
        reaction_error = f"The coverage of reactions is too low ({reaction_coverage*100:.1f}%), and it is not recommended to use DLKcat to obtain enzyme kinetic data."
        error_list.append(reaction_error)

    sui_or_not='Yes'
    if len(error_list) > 0:
        sui_or_not='No'
    else:
        error_list.append("Suitable for constructing enzyme-bound models.")

    return (sui_or_not,error_list)

def get_reaction_kcatmw_onestop_by_AutoPACMEN(autopacmen_folder,sbml_path,bigg_metabolites_file,brenda_textfile_path,project_name,uniprot_data_file,organism,protein_kcat_database_path,kcat_gap_fill,reaction_gap_fill):
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

    # Step 1: get bigg metbolite
    print("Starting to deal BIGG metabolites text file...")
    parse_bigg_metabolites_file(bigg_metabolites_file, autopacmen_folder)
    print("BIGG metabolites text file done!")

    # Step 2: BRENDA kcat
    print("Starting to deal BRENDA textfile...")
    parse_brenda_textfile(brenda_textfile_path, autopacmen_folder, brenda_json_path, brenda_json_path2) 
    print("BRENDA textfile done!")

    # Step 3: Select Brenda kcat for model
    print("Starting to deal brenda json for model...")
    parse_brenda_json_for_model(sbml_path, brenda_json_path, brenda_output_json_path)
    print("BRENDA json for model done!")

    # Step 4: SABIO-RK kcat for model
    print("Starting EC numbers kcat search in SABIO-RK...")
    parse_sabio_rk_for_model_with_sbml(sbml_path, sabio_rk_json_path, bigg_id_name_mapping_path)
    print("SABIO-RK done!")

    # Step 5: Brenda and SABIO-RK kcat combined
    print("Combining kcat database...")
    create_combined_kcat_database(sabio_rk_json_path, brenda_output_json_path, combined_output_path)
    print("Combining kcat database done!")

    # Step 6: subunit number of each reaction
    print("Starting to fetch subunit number of each enzyme")
    #model=cobra.io.read_sbml_model(sbml_path)
    if re.search('\.xml',sbml_path):
        model = cobra.io.read_sbml_model(sbml_path)
    elif re.search('\.json',sbml_path):
        model = cobra.io.json.load_json_model(sbml_path)   
    get_gene_subunitDescription(sub_description_path,model)#从uniprot的api下载，运行一次就行
    subbnumdf = get_subunit_number(sub_description_path,gene_subnum_path)
    print("Calculation done!")

    # Step 7: get mw for model gene (must be uniprot ID)
    print("Starting UniProt ID<->Protein mass search using UniProt...")
    get_protein_mass_mapping_from_local(sbml_path, autopacmen_folder, project_name, uniprot_data_file)
    get_reaction_mw(sbml_path,autopacmen_folder, project_name, reaction_mw_path, gene_subnum_path)
    print("Protein ID<->Mass mapping done!")

    # Step 8: kcat assignment for model(include sa)
    print("Starting to assign kcat for model...")
    get_reactions_kcat_mapping(sbml_path, autopacmen_folder, project_name, organism, combined_output_path,brenda_json_path2, reaction_mw_path,protein_kcat_database_path,kcat_gap_fill)
    print("kcat assignment done!")

    # Step 9: get_reaction_kcat_mw for model
    print("Starting to get reaction kcat_mw for model...")
    if re.search('\.xml',sbml_path):
        model = cobra.io.read_sbml_model(sbml_path)
    elif re.search('\.json',sbml_path):
        model = cobra.io.json.load_json_model(sbml_path)   
    get_reaction_kcat_mw(model,autopacmen_folder, project_name, reaction_gap_fill,gene_subnum_path,reaction_kcat_mw_path)       
    print("Reaction kcat_mw done!")
    return reaction_kcat_mw_path

def get_reaction_kcatmw_onestop_by_DLKcat(dlkcat_folder,sbml_path):
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
    # Step 0: read GEM
    if re.search('\.xml',sbml_path):
        model = cobra.io.read_sbml_model(sbml_path)
    elif re.search('\.json',sbml_path):
        model = cobra.io.json.load_json_model(sbml_path)

    # Step 1: subunit number of each reaction
    print("Starting to fetch subunit number of each enzyme")
    get_gene_subunitDescription(sub_description_path,model)#Download from the UniProt API, run it once.
    subbnumdf = get_subunit_number(sub_description_path,gene_subnum_path)
    print("Calculation done!")

    # Step 2: convert metbolites bigg id to smiles 
    print("Starting to convert metbolites bigg id to smiles...")
    metdf_name = get_met_bigg_id(model)
    inchkeydf = convert_bigg_met_to_inchikey(metdf_name['met'],inchikey_list_file)#from BIGG
    inchkeydf = pd.read_csv('./data/inchikey_list.csv')
    smilesdf = convert_inchikey_to_smiles(inchkeydf,inchikey_list_smilesfile)#from pubchem
    print("Converting done!")
    
    # Step 3: get protein sequence and mass in model 
    print("Starting to get protein sequence and mass in model...")
    subbnumdf = pd.read_csv(gene_subnum_path,index_col=0)
    prodf = get_model_protein_sequence_and_mass(model,subbnumdf,prodf_file)
    print("Getting done!")

    # Step 4: split the substrate of reactions to match the gene
    print("Starting to split the substrate of reactions to match the gene...")
    spdf = split_substrate_to_match_gene(model,metabolites_reactions_gpr_file)
    print("Splitting done!")

    # Step 5: combine the reaction--substrate--gene--protein_sequnce--mass and formate DLKcat input file
    print("Starting to combine data...")
    comdf = combine_reactions_simles_sequence(spdf,smilesdf,prodf,comdf_file)
    DLinputdf = generate_DLKCAT_input(comdf,metdf_name,metdf_outfile,DLinput_file)
    print("Combinning done!")

    # Step 6: use DLKcat calculate kcat
    print("Starting to Use DLKcat calculate kcat...")
    cmd_str = "python ./script/prediction_for_input.py %s/DLinput.tsv %s/DLoutput.tsv"%(dlkcat_folder,dlkcat_folder)
    subprocess.run(cmd_str, shell=True)
    print("DLKcat done!")

    # Step 7: get the kcat_mw file
    print("Starting to get reaction kcat_mw for model......")
    DLouputdf = pd.read_csv(DLouputdf_file, sep='\t')
    DL_reaction_kact_mw = DL_kcat_mw_calculation(DLouputdf, comdf)
    DL_reaction_kact_mw.to_csv(DL_reaction_kact_mw_file, index=False)
    print("Reaction kcat_mw done!")    
    return DL_reaction_kact_mw_file

def get_reaction_kcatmw_onestop(work_folder, kcat_method, model_file, bigg_metabolites_file, brenda_textfile_path, uniprot_data_file, organism, kcat_gap_fill, reaction_gap_fill):
    if kcat_method == 'AutoPACMEN':
        autopacmen_folder = work_folder+'_by_AutoPACMEN/'
        create_file(autopacmen_folder)
        project_name = "model_%s"%kcat_gap_fill
        protein_kcat_database_path = "none"
        reaction_kcat_MW_file = get_reaction_kcatmw_onestop_by_AutoPACMEN(autopacmen_folder,model_file,bigg_metabolites_file,brenda_textfile_path,project_name,uniprot_data_file,organism,protein_kcat_database_path,kcat_gap_fill,reaction_gap_fill)
    elif kcat_method == 'DLKcat':
        dlkcat_folder = work_folder+'_by_DLKcat/'
        create_file(dlkcat_folder)
        reaction_kcat_MW_file = get_reaction_kcatmw_onestop_by_DLKcat(dlkcat_folder,model_file)
    else:
        print('The method you provided is not supported, please check the spelling!')
        
    return reaction_kcat_MW_file