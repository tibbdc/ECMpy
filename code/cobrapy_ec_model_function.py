# -*- coding: utf-8 -*-
# This code is used to introduce enzyme concentration constraint in GEMs
# by COBRApy and to calculate the parameters that need to be entered
# during the construction of the enzyme-constrained model.
#from warnings import warn

import pandas as pd
import numpy as np
import json
import cobra
import math
import re
import random
import statistics
import os
import shutil
from cobra.core import Reaction
from cobra.io.dict import model_to_dict
from cobra.util.solver import set_objective
from xml.dom import minidom
from optlang.symbolics import Zero, add

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
            reverse_reaction._gene_reaction_rule = reaction._gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    model.add_reactions(reactions_to_add)
    set_objective(model, coefficients, additive=True)


def get_genes_and_gpr(model,gene_outfile,gpr_outfile):
    """Retrieving genes and gene_reaction_rule from GEM.

    Arguments
    ----------
    * model: cobra.Model ~ A genome scale metabolic network model for
        constructing the enzyme-constrained model.

    :return: all genes and gpr in model.
    """
    model_dict = model_to_dict(model, sort=False)
    genes = pd.DataFrame(model_dict['genes']).set_index(['id'])
    genes.to_csv(gene_outfile)
    all_gpr = pd.DataFrame(model_dict['reactions']).set_index(['id'])
    all_gpr.to_csv(gpr_outfile)
    return [genes, all_gpr]

def get_reaction_gene_subunit_MW(reaction_gene_subunit_file,gene_molecular_weight_file,save_file):
    """Retrieving genes,subunits and MW, and split 'or' type of reaction

    Arguments
    ----------
    * reaction_gene_subunit_file: gene-molecular_weight file eg. b3500,48771.94
    * gene_molecular_weight_file: manually get the subunit of each protein from EcoCy  

    :return: all reaction with gene, subunit, and MW.
    """
    reaction_gene_subunit_MW_new = pd.DataFrame()
    reaction_gene_subunit = pd.read_csv(reaction_gene_subunit_file, index_col=0)
    protein_mw=pd.read_csv(gene_molecular_weight_file, index_col=0)
    for reaction, data in reaction_gene_subunit.iterrows():
        if re.search(" or ", data['gene_reaction_rule']):
            gene = enumerate(data['gene_reaction_rule'].split(" or "))
            subunit_num = data['subunit_num'].split(" or ")
            for index, value in gene:
                if index == 0:
                    reaction_new = reaction + "_num1"
                    reaction_gene_subunit_MW_new.loc[reaction_new,'name'] = data['name']
                    reaction_gene_subunit_MW_new.loc[reaction_new, 'gene_reaction_rule'] = value
                    reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_num'] = subunit_num[index]
                    if re.search(" and ", value):
                        reaction_gene_subunit_MW = []
                        gene2 = enumerate(value.replace('(', '').replace(")", '').replace(" ", '').split('and'))
                        for index2, value2 in gene2:
                            reaction_gene_subunit_MW.append(str(protein_mw.loc[value2,'mw']))
                        reaction_gene_subunit_MW=' and '.join(reaction_gene_subunit_MW)
                        if re.search('\(',value):
                            reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = '( '+reaction_gene_subunit_MW+ ' )'
                        else:
                            reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = reaction_gene_subunit_MW
                    else:
                        reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = protein_mw.loc[value,'mw']
                else:
                    reaction_new = reaction + "_num" + str(index+1)
                    reaction_gene_subunit_MW_new.loc[reaction_new, 'name'] = data['name']
                    reaction_gene_subunit_MW_new.loc[reaction_new, 'gene_reaction_rule'] = value
                    reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_num'] = subunit_num[index]
                    if re.search(" and ", value):
                        reaction_gene_subunit_MW = []
                        gene3 = enumerate(value.replace('(', '').replace(")", '').replace(" ", '').split('and'))
                        for index3, value3 in gene3:
                            reaction_gene_subunit_MW.append(str(protein_mw.loc[value3,'mw']))
                        reaction_gene_subunit_MW=' and '.join(reaction_gene_subunit_MW)
                        if re.search('\(',value3):
                            reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = '( '+reaction_gene_subunit_MW+ ' )'
                        else:
                            reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = reaction_gene_subunit_MW
                    else:
                        reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = protein_mw.loc[value,'mw']                                     
        elif re.search(" and ", data['gene_reaction_rule']):
            reaction_gene_subunit_MW = []
            gene4 = enumerate(data['gene_reaction_rule'].replace('(', '').replace(")", '').replace(" ", '').split('and'))
            reaction_gene_subunit_MW_new.loc[reaction, 'name'] = data['name']
            reaction_gene_subunit_MW_new.loc[reaction, 'gene_reaction_rule'] = data['gene_reaction_rule']
            reaction_gene_subunit_MW_new.loc[reaction,'subunit_num'] = data['subunit_num']
            for index4, value4 in gene4:
                reaction_gene_subunit_MW.append(str(protein_mw.loc[value4,'mw']))
            reaction_gene_subunit_MW=' and '.join(reaction_gene_subunit_MW)
            if re.search('\(',value4):
                reaction_gene_subunit_MW_new.loc[reaction,'subunit_mw'] = '( '+reaction_gene_subunit_MW+ ' )'
            else:
                reaction_gene_subunit_MW_new.loc[reaction,'subunit_mw'] = reaction_gene_subunit_MW
        else:
            reaction_gene_subunit_MW_new.loc[reaction, 'name'] = data['name']
            reaction_gene_subunit_MW_new.loc[reaction,
                                             'gene_reaction_rule'] = data['gene_reaction_rule']
            reaction_gene_subunit_MW_new.loc[reaction,
                                             'subunit_mw'] = protein_mw.loc[data['gene_reaction_rule'],'mw']
            reaction_gene_subunit_MW_new.loc[reaction,
                                             'subunit_num'] = data['subunit_num']
    reaction_gene_subunit_MW_new.to_csv(save_file)
    return reaction_gene_subunit_MW_new

def calculate_reaction_mw(reaction_gene_subunit_MW,reaction_mw_outfile):
    """Calculate the molecular weight of the enzyme that catalyzes each
    reaction in GEM based on the number of subunits and
    molecular weight of each gene.

    Arguments
    ----------
    * reaction_gene_subunit_MW: A CSV file contains the GPR relationship
     for each reaction in the GEM model,the number of subunit components 
     of each gene expressed protein, and the molecular weight of each 
     gene expressed protein.

    :return: The molecular weight of the enzyme that catalyzes each reaction
     in the GEM model.
    """
    reaction_gene_subunit_MW = pd.read_csv(
        reaction_gene_subunit_MW, index_col=0)
    reaction_mw = pd.DataFrame()
    for reaction_id in reaction_gene_subunit_MW.index:
        subunit_mw_list = reaction_gene_subunit_MW.loc[reaction_id, 'subunit_mw'].\
            replace('(', '').replace(")", '').replace(" ", '').split('or')
        subunit_num_list = reaction_gene_subunit_MW.loc[reaction_id, 'subunit_num'].\
            replace('(', '').replace(")", '').replace(" ", '').split('or')

        mw_s = ''
        for mw_i in range(0, len(subunit_mw_list)):
            mw_list = np.array(subunit_mw_list[mw_i].split('and'))
            num_list = np.array(subunit_num_list[mw_i].split('and'))
            mw_list = list(map(float, mw_list))
            num_list = list(map(float, num_list))
            mw_s = mw_s + \
                str(round(np.sum(np.multiply(mw_list, num_list)), 4)) + ' or '

        reaction_mw.loc[reaction_id, 'MW'] = mw_s.rstrip(' or ')
    reaction_mw.to_csv(reaction_mw_outfile)
    return reaction_mw


def calculate_reaction_mw_not_consider_subunit(reaction_gene_subunit_MW, save_file):
    """Calculate the molecular weight of the enzyme that catalyzes each
    reaction in GEM based on the number of subunits and
    molecular weight of each gene.

    Arguments
    ----------
    * reaction_gene_subunit_MW: A CSV file contains the GPR relationship
     for each reaction in the GEM model,the number of subunit components 
     of each gene expressed protein, and the molecular weight of each 
     gene expressed protein.

    :return: The molecular weight of the enzyme that catalyzes each reaction
     in the GEM model.
    """
    reaction_gene_subunit_MW = pd.read_csv(
        reaction_gene_subunit_MW, index_col=0)
    reaction_mw = pd.DataFrame()
    for reaction_id in reaction_gene_subunit_MW.index:
        subunit_mw_list = reaction_gene_subunit_MW.loc[reaction_id, 'subunit_mw'].\
            replace('(', '').replace(")", '').replace(" ", '').split('or')
        subunit_num_list = reaction_gene_subunit_MW.loc[reaction_id, 'subunit_num'].\
            replace('(', '').replace(")", '').replace(" ", '').split('or')

        mw_s = ''
        for mw_i in range(0, len(subunit_mw_list)):
            mw_list = np.array(subunit_mw_list[mw_i].split('and'))
            num_list = np.array(subunit_num_list[mw_i].split('and'))
            mw_list = list(map(float, mw_list))
            mw_s = mw_s + \
                str(round(np.sum(np.multiply(mw_list, 1)), 4)) + ' or '

        mw_s = mw_s.rstrip(' or ')
        reaction_mw.loc[reaction_id, 'MW'] = mw_s
    reaction_mw.to_csv(save_file)
    return reaction_mw


def calculate_reaction_kcat_mw(reaction_kcat_file, reaction_MW_file, save_file, select_key):
    """Calculating kcat/MW

    Arguments
    ----------
    * reaction_kcat_file: A CSV file contains the kcat values for each
    reaction in the model.
    * reaction_mw: The molecular weight of the enzyme that catalyzes
     each reaction in the GEM model.

    :return: The kcat/MW value of the enzyme catalyzing each reaction
     in the GEM model.
    """
    reaction_kcat = json_load(reaction_kcat_file)
    reaction_kcat_mw = pd.DataFrame()
    reaction_MW = pd.read_csv(reaction_MW_file, index_col=0)
    for reaction_idmw in reaction_MW.index:
        if re.search('_num', reaction_idmw):
            reaction_id = reaction_idmw.split('_num')[0]
        else:
            reaction_id = reaction_idmw
        for key,value in reaction_kcat.items():
            if re.search('_b',key):
                reaction_kcat_id = key.split('_b')[0]+'_reverse'
                if reaction_id == reaction_kcat_id and str(value[select_key]) != 'nan':
                    reaction_kcat_mw.loc[reaction_idmw, 'MW'] = reaction_MW.loc[reaction_idmw, 'MW']
                    reaction_kcat_mw.loc[reaction_idmw,'kcat'] = value[select_key]*3600
                    kcat_MW = value[select_key]*3600/reaction_MW.loc[reaction_idmw, 'MW']
                    reaction_kcat_mw.loc[reaction_idmw, 'kcat_MW'] = kcat_MW
            elif re.search('_f',key):
                reaction_kcat_id = key.split('_f')[0]
                if reaction_id == reaction_kcat_id and str(value[select_key]) != 'nan':
                    reaction_kcat_mw.loc[reaction_idmw, 'MW'] = reaction_MW.loc[reaction_idmw, 'MW']
                    reaction_kcat_mw.loc[reaction_idmw,'kcat'] = value[select_key]*3600
                    kcat_MW = value[select_key]*3600/reaction_MW.loc[reaction_idmw, 'MW']
                    reaction_kcat_mw.loc[reaction_idmw, 'kcat_MW'] = kcat_MW
            elif reaction_id == key and str(value[select_key]) != 'nan':
                reaction_kcat_mw.loc[reaction_idmw, 'MW'] = reaction_MW.loc[reaction_idmw, 'MW']
                reaction_kcat_mw.loc[reaction_idmw,'kcat'] = value[select_key]*3600
                kcat_MW = value[select_key]*3600/reaction_MW.loc[reaction_idmw, 'MW']
                reaction_kcat_mw.loc[reaction_idmw, 'kcat_MW'] = kcat_MW                
    reaction_kcat_mw.to_csv(save_file)
    return reaction_kcat_mw

def calculate_f(genes, gene_abundance_file, subunit_molecular_weight_file):
    """Calculating f (the mass fraction of enzymes that are accounted
    in the model out of all proteins) based on the protein abundance
    which can be obtained from PAXdb database.

    Arguments
    ----------
    * genes: All the genes in the model.
    * gene_abundance_file: The protein abundance of each gene
     in the E. coli genome.
    * subunit_molecular_weight_file: The molecular weight of the
     protein subunit expressed by each gene.

    :return: The enzyme mass fraction f.
    """
    gene_abundance = pd.read_csv(gene_abundance_file, index_col=0)
    subunit_molecular_weight = pd.read_csv(
        subunit_molecular_weight_file, index_col=0)
    enzy_abundance = 0
    pro_abundance = 0
    for gene_i in gene_abundance.index:
        abundance = gene_abundance.loc[gene_i, 'abundance'] * \
            subunit_molecular_weight.loc[gene_i, 'mw']
        pro_abundance += abundance
        if gene_i in genes.index:
            enzy_abundance += abundance
    f = enzy_abundance/pro_abundance
    return f


def set_enzyme_constraint(model, reaction_kcat_mw, lowerbound, upperbound):
    """Introducing enzyme concentration constraint
    by COBRApy using the calculated parameters.

    Arguments
    ----------
    * model: cobra.Model ~ A genome scale metabolic network model for
        constructing the enzyme-constrained model.
    * reaction_kcat_mw: The kcat/MW value of the enzyme catalyzing each
     reaction in the GEM model.
    * lowerbound: The lower bound of enzyme concentration constraint in
     the enzyme-constrained model.
    * upperbound: The upper bound of enzyme concentration constraint in
     the enzyme-constrained model.

    :return: Construct an enzyme-constrained model.
    """
    coefficients = dict()
    for rxn in model.reactions:
        if rxn.id in reaction_kcat_mw.index:
            coefficients[rxn.forward_variable] = 1 / \
                float(reaction_kcat_mw.loc[rxn.id, 'kcat_MW'])
    constraint = model.problem.Constraint(0, lb=lowerbound, ub=upperbound)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    return model

def create_file(store_path):
    """Create flod.

    Arguments
    ----------
    * store_path: The path of need create file
    """
    if os.path.exists(store_path):
        print("path exists")
        #shutil.rmtree(store_path)
        #os.makedirs(store_path)
    else:      
        os.makedirs(store_path)
        print(store_path) 
        
def json_load(path):
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary


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
    model = cobra.io.read_sbml_model(model_file)
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
    reaction_kcay_mw_dict = {}
    for eachreaction in range(len(dictionary_model['reactions'])):
        reaction_id = dictionary_model['reactions'][eachreaction]['id']
        if reaction_id in reaction_kcat_mw.index:
            dictionary_model['reactions'][eachreaction]['kcat'] = reaction_kcat_mw.loc[reaction_id, 'kcat']
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
            reaction_kcay_mw_dict[reaction_id] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
        else:
            dictionary_model['reactions'][eachreaction]['kcat'] = ''
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = ''
    dictionary_model['enzyme_constraint']['kcat_MW'] = reaction_kcay_mw_dict
    json_write(json_output_file, dictionary_model)


def trans_model2enz_json_model(model_file, reaction_kcat_mw_file, f, ptot, sigma, lowerbound, upperbound, json_output_file):
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
    * json_output_file: json file store json model
    """
    model = cobra.io.read_sbml_model(model_file)
    convert_to_irreversible(model)
    model_name = model_file.split('/')[-1].split('.')[0]
    json_path = "./model/%s_irreversible.json" % model_name
    cobra.io.save_json_model(model, json_path)
    dictionary_model = json_load(json_path)
    dictionary_model['enzyme_constraint'] = {'enzyme_mass_fraction': f, 'total_protein_fraction': ptot,
                                             'average_saturation': sigma, 'lowerbound': lowerbound, 'upperbound': upperbound}
    # Reaction-kcat_mw file.
    # eg. AADDGT,49389.2889,40.6396,1215.299582180927
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)
    reaction_kcay_mw_dict = {}
    for eachreaction in range(len(dictionary_model['reactions'])):
        reaction_id = dictionary_model['reactions'][eachreaction]['id']
        if reaction_id in reaction_kcat_mw.index:
            dictionary_model['reactions'][eachreaction]['kcat'] = reaction_kcat_mw.loc[reaction_id, 'kcat']
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
            reaction_kcay_mw_dict[reaction_id] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
        else:
            dictionary_model['reactions'][eachreaction]['kcat'] = ''
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = ''
    dictionary_model['enzyme_constraint']['kcat_MW'] = reaction_kcay_mw_dict
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
        if rxn.id in dictionary_model['enzyme_constraint']['kcat_MW'].keys():
            coefficients[rxn.forward_variable] = 1 / \
                float(dictionary_model['enzyme_constraint']['kcat_MW'][rxn.id])

    lowerbound = dictionary_model['enzyme_constraint']['lowerbound']
    upperbound = dictionary_model['enzyme_constraint']['upperbound']
    constraint = model.problem.Constraint(0, lb=lowerbound, ub=upperbound)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    return model


def get_enzyme_usage(enz_total, reaction_flux_file, reaction_kcat_mw_file, reaction_enz_usage_file):
    """Get the enzyme usage of each reaction

    Arguments
    ----------
    * enz_total: total enzyme amount(e.g. 0.227).
    * reaction_flux_file: reaction flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_enz_usage_file: the output file contain reaction and enzyme usage of each reaction.

    :return: reaction_enz_usage_file.
    """
    reaction_fluxes = pd.read_csv(reaction_flux_file, index_col=0)
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)

    reaction_enz_usage_df = pd.DataFrame()
    for index, row in reaction_kcat_mw.iterrows():
        if index in reaction_fluxes.index:
            reaction_enz_usage_df.loc[index, 'kcat_mw'] = row['kcat_MW']
            reaction_enz_usage_df.loc[index,
                                      'flux'] = reaction_fluxes.loc[index, 'fluxes']
            reaction_enz_usage_df.loc[index,
                                      'enz useage'] = reaction_fluxes.loc[index, 'fluxes']/row['kcat_MW']
            reaction_enz_usage_df.loc[index, 'enz ratio'] = reaction_fluxes.loc[index,
                                                                                'fluxes']/row['kcat_MW']/enz_total

    reaction_enz_usage_df = reaction_enz_usage_df.sort_values(
        by="enz ratio", axis=0, ascending=False)
    reaction_enz_usage_df.to_csv(reaction_enz_usage_file)
    return reaction_enz_usage_df

    
def get_fluxes_detail_in_model(model, fluxes_outfile, reaction_kcat_mw_file):
    """Get the detailed information of each reaction

    Arguments
    ----------
    * model: cobra.Model.
    * fluxes_outfile: reaction flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.

    :return: fluxes, kcat, MW and kcat_MW in dataframe.
    """
    model_pfba_solution = cobra.flux_analysis.pfba(model)
    model_pfba_solution = model_pfba_solution.to_frame()
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)
    model_pfba_solution_detail = pd.DataFrame()
    for index, row in model_pfba_solution.iterrows():
        reaction_detail = model.reactions.get_by_id(index)
        model_pfba_solution_detail.loc[index, 'fluxes'] = row['fluxes']
        if index in reaction_kcat_mw.index:
            model_pfba_solution_detail.loc[index,
                                           'kcat'] = reaction_kcat_mw.loc[index, 'kcat']
            model_pfba_solution_detail.loc[index,
                                           'MW'] = reaction_kcat_mw.loc[index, 'MW']
            model_pfba_solution_detail.loc[index,
                                           'kcat_MW'] = reaction_kcat_mw.loc[index, 'kcat_MW']
            if 'source' in reaction_kcat_mw.columns:
                model_pfba_solution_detail.loc[index,
                                               'source'] = reaction_kcat_mw.loc[index, 'source']
        model_pfba_solution_detail.loc[index, 'equ'] = reaction_detail.reaction
    model_pfba_solution_detail.to_csv(fluxes_outfile)
    return model_pfba_solution


def change_reaction_kcat_by_database(json_model_path,select_reactionlist, kcat_database_combined_file, reaction_kcat_mw_file, reaction_kapp_change_file):
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
    Brenda_sabio_combined_select = json_load(kcat_database_combined_file)
    json_model=cobra.io.load_json_model(json_model_path)
    reaction_change_accord_database = []
    for eachreaction in select_reactionlist:
        select_reaction = json_model.reactions.get_by_id(eachreaction)
        if "ec-code" in select_reaction.annotation.keys():
            ec_number = select_reaction.annotation["ec-code"]
            kcat_max_list = []
            if isinstance(ec_number, str):
                if ec_number in Brenda_sabio_combined_select.keys():
                    reaction_kcat_max = Brenda_sabio_combined_select[ec_number]['kcat']
                    if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max * 3600:
                        reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max * 3600
                        reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600/reaction_kcat_mw.loc[eachreaction, 'MW']
                        reaction_change_accord_database.append(eachreaction) 
            else:
                for eachec in ec_number:
                    if eachec in Brenda_sabio_combined_select.keys():
                        kcat_max_list.append(Brenda_sabio_combined_select[eachec]['kcat'])
                reaction_kcat_max = np.max(kcat_max_list)     
                if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max * 3600:
                    reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max * 3600
                    reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600/reaction_kcat_mw.loc[eachreaction, 'MW']
                    reaction_change_accord_database.append(eachreaction)              
    reaction_kcat_mw.to_csv(reaction_kapp_change_file)
    return(reaction_change_accord_database)


def get_enz_model_use_enz_usage_by_eckcat(enz_ratio,json_model_path, reaction_flux_file, reaction_kcat_mw_file, reaction_enz_usage_file, kcat_database_combined_file, model_file, f, ptot, sigma, lowerbound, upperbound, json_output_file, reaction_mw_outfile):

    """Get new enzyme model using enzyme usage to calibration

    Arguments
    ----------
    * enz_ratio: enzyme ratio which needed change.
    * json_model_path: The file storing json model.
    * reaction_flux_file: reaction-flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_enz_usage_file： enzyme usage of each reaction.
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
    reaction_enz_usage_df = get_enzyme_usage(
        upperbound, reaction_flux_file, reaction_kcat_mw_file, reaction_enz_usage_file)

    select_reaction = list(
        reaction_enz_usage_df[reaction_enz_usage_df['enz ratio'] > enz_ratio].index)  # more than 1%
    print('need changing reaction: ')
    print(select_reaction)
    change_reaction_list_round1 = change_reaction_kcat_by_database(json_model_path,select_reaction, kcat_database_combined_file, reaction_kcat_mw_file, reaction_mw_outfile)
    print('changed reaction: ')
    print(change_reaction_list_round1)

    trans_model2enz_json_model_split_isoenzyme(
        model_file, reaction_mw_outfile, f, ptot, sigma, lowerbound, upperbound, json_output_file)

    enz_model = get_enzyme_constraint_model(json_output_file)
    return enz_model


def select_calibration_reaction_by_c13(reaction_kcat_mw_file, c13reaction_file, enzyme_amount, percentage, sigma):
    """Get reaction list need change kcat using c13 data

    Arguments
    ----------
    * reaction_kcat_mw_file: original file stored kcat/mw of reaction.
    * c13reaction_file: The file contained reaction list whcih had c13 flux.
    * enzyme_amount:  total enzyme amount(e.g. 0.227).
    * percentage:  percentage which needed change.
    * sigma: The approximated average saturation of enzyme. 
    
    :return: fluxes, kcat, MW and kcat_MW in dataframe.
    """    
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)
    c13reaction = pd.read_csv(c13reaction_file, index_col=0)
    c13reaction_selecet = []
    for index, row in c13reaction.iterrows():
        if index in reaction_kcat_mw.index:
            ECMpy_c13_reaction_flux = reaction_kcat_mw.loc[index,
                                                           'kcat_MW']*enzyme_amount*percentage*sigma
            if ECMpy_c13_reaction_flux < row['Flux norm']:
                c13reaction_selecet.append(index)
    return(c13reaction_selecet)


def reaction_gene_subunit_MW_split(reaction_gene_subunit_MW,save_file):
    """Split reaction MW accord subunit

    Arguments
    ----------
    * reaction_gene_subunit_MW: reaction gene-mw file contained subunit information.
    * save_file: new reaction gene-mw file. 
    
    """  
    reaction_gene_subunit_MW_new = pd.DataFrame()
    for reaction, data in reaction_gene_subunit_MW.iterrows():
        if re.search(" or ", data['gene_reaction_rule']):
            gene = enumerate(data['gene_reaction_rule'].split(" or "))
            subunit_mw = data['subunit_mw'].split(" or ")
            subunit_num = data['subunit_num'].split(" or ")
            for index, value in gene:
                if index == 0:
                    reaction_new = reaction + "_num1"
                    reaction_gene_subunit_MW_new.loc[reaction_new,
                                                     'name'] = data['name']
                    reaction_gene_subunit_MW_new.loc[reaction_new,
                                                     'gene_reaction_rule'] = value
                    reaction_gene_subunit_MW_new.loc[reaction_new,
                                                     'subunit_mw'] = subunit_mw[index]
                    reaction_gene_subunit_MW_new.loc[reaction_new,
                                                     'subunit_num'] = subunit_num[index]
                else:
                    reaction_new = reaction + "_num" + str(index+1)
                    reaction_gene_subunit_MW_new.loc[reaction_new,
                                                     'name'] = data['name']
                    reaction_gene_subunit_MW_new.loc[reaction_new,
                                                     'gene_reaction_rule'] = value
                    reaction_gene_subunit_MW_new.loc[reaction_new,
                                                     'subunit_mw'] = subunit_mw[index]
                    reaction_gene_subunit_MW_new.loc[reaction_new,
                                                     'subunit_num'] = subunit_num[index]
        else:
            reaction_gene_subunit_MW_new.loc[reaction, 'name'] = data['name']
            reaction_gene_subunit_MW_new.loc[reaction,
                                             'gene_reaction_rule'] = data['gene_reaction_rule']
            reaction_gene_subunit_MW_new.loc[reaction,
                                             'subunit_mw'] = data['subunit_mw']
            reaction_gene_subunit_MW_new.loc[reaction,
                                             'subunit_num'] = data['subunit_num']
    reaction_gene_subunit_MW_new.to_csv(save_file)


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


def get_diff_reaction_use_c13(c13reaction_file, model_fluxes):
    """Split isoenzyme reaction to mutiple reaction

    Arguments
    ----------
    * c13reaction_file: The file contained reaction list whcih had c13 flux.
    * model_fluxes: calulated fluxes
    
    :return: reaction list.
    """   
    c13reaction = pd.read_csv(c13reaction_file, index_col=0)
    c13reaction = list(c13reaction.index)
    enz_model_pfba_solution_select = model_fluxes[model_fluxes['fluxes'] > 0]
    enz_model_pfba_solution_select_id = []
    for eachreaction in enz_model_pfba_solution_select.index:
        if re.search('_num', eachreaction):
            enz_model_pfba_solution_select_id.append(
                eachreaction.split('_num')[0])
        else:
            enz_model_pfba_solution_select_id.append(eachreaction)
    c13reaction_2_enz_model_diff = list(
        set(c13reaction).difference(set(enz_model_pfba_solution_select_id)))
    return(c13reaction_2_enz_model_diff)


def get_enz_model_use_c13(reaction_kcat_mw_file,json_model_path, c13reaction_file, percentage, kcat_database_combined_file, model_file, f, ptot, sigma, lowerbound, upperbound, json_output_file, reaction_mw_outfile):
    """Get new enzyme model using C13 reaction to calibration

    Arguments
    ----------
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * json_model_path: The file storing json model.
    * c13reaction_file: The file contained reaction list whcih had c13 flux.
    * percentage:  percentage which needed change.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.  
    * model_file: cobra model.
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint.  
    * json_output_file: json file store json model
    * reaction_mw_outfile: changed file stored reaction kcat/mw.

    :return: new enzyme model.
    """
    c13reaction_selecet = select_calibration_reaction_by_c13(
        reaction_kcat_mw_file, c13reaction_file, upperbound, percentage, sigma)
    print('need changing reaction: ')
    print(c13reaction_selecet)

    #if isinstance(df_reaction_select, pd.DataFrame):
    #    reaction_kcat_mw_file = "./analysis/reaction_change_by_biomass.csv"

    reaction_kapp_change_file = reaction_mw_outfile
    #c13reaction_selecet=['CS','ACONTa','ACONTb','ICDHyr','MALS', 'MDH', 'ICL', 'SUCOAS_reverse', 'SUCDi', 'AKGDH']
    change_reaction_list_round1 = change_reaction_kcat_by_database(json_model_path,
        c13reaction_selecet, kcat_database_combined_file, reaction_kcat_mw_file, reaction_kapp_change_file)
    print('changed reaction: ')
    print(change_reaction_list_round1)

    reaction_kcat_mw_file = reaction_mw_outfile
    trans_model2enz_json_model_split_isoenzyme(
        model_file, reaction_kcat_mw_file, f, ptot, sigma, lowerbound, upperbound, json_output_file)
    enz_model = get_enzyme_constraint_model(json_output_file)
    return enz_model


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
            if r.id in dict_coeff:
                dict_obj[r.forward_variable] = 1 / dict_coeff[r.id]

        model_obj = model.problem.Objective(Zero, direction="min", sloppy=True)
        model.objective = model_obj
        model.objective.set_linear_coefficients(dict_obj)

        solution = model.optimize()
    return solution