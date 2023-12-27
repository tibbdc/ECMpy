import json
import os
import cobra
import copy
import csv
import io
import requests
import time
import math
import random
import pickle
import sys
import statistics
from typing import Any, Dict, List
from Bio import Entrez
from Bio.SeqUtils.ProtParam import ProteinAnalysis  
import re#Mr.Mao

def create_file(store_path):
    if os.path.exists(store_path):
        print("path exists")
        #shutil.rmtree(store_path)
        #os.makedirs(store_path)
    else:      
        os.makedirs(store_path)
        print(store_path) 
        
def pickle_write(path: str, pickled_object: Any) -> None:
    """Writes the given object as pickled file with the given path

    Arguments
    ----------
    * path: str ~ The path of the pickled file that shall be created
    * pickled_object: Any ~ The object which shall be saved in the pickle file
    """
    pickle_file = open(path, 'wb')
    pickle.dump(pickled_object, pickle_file)
    pickle_file.close()
    
def pickle_load(path: str) -> Any:
    """Returns the value of the given pickle file.

    Arguments
    ----------
    * path: str ~ The path to the pickle file.
    """
    pickle_file = open(path, 'rb')
    pickled_object = pickle.load(pickle_file)
    pickle_file.close()
    return pickled_object    

def json_write(path: str, dictionary: Dict[Any, Any]) -> None:
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path: str ~  The path of the JSON file that shall be written
    * dictionary: Dict[Any, Any] ~ The dictionary which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(dictionary, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)

def parse_bigg_metabolites_file(bigg_metabolites_file_path: str, json_output_folder: str) -> None:
    """Parses a BIGG metabolites text file and returns a dictionary for this file.

    As of 29/04/2019, a BIGG metabolites list of all BIGG-included metabolites
    is retrievable under http://bigg.ucsd.edu/data_access

    Arguments
    ----------
    * bigg_metabolites_file_path: str ~ The file path to the BIGG metabolites file.
      The usual file name (which has to be included too in this argument) is
      bigg_models_metabolites.txt
    * output_folder: str ~ The folder in which the JSON including the parsed BIGG
      metabolites file data is stored with the name 'bigg_id_name_mapping.json'

    Output
    ----------
    * A JSON file with the name 'bigg_id_name_mapping.json' in the given output folder,
      with the following structure:
    <pre>
     {
         "$BIGG_ID": "$CHEMICAL_OR_USUAL_NAME",
         (...),
         "$BIGG_ID": "$BIGG_ID",
         (...),
     }
    </pre>
    The BIGG ID <-> BIGG ID mapping is done for models which already use the BIGG IDs.
    """
    # Standardize output folder
    json_output_folder = standardize_folder(json_output_folder)

    # Open the BIGG metabolites file as string list, and remove all newlines
    with open(bigg_metabolites_file_path, "r") as f:
        lines = f.readlines()
    lines = [x.replace("\n", "") for x in lines if len(x) > 0]

    # Mapping variable which will store the BIGG ID<->
    bigg_id_name_mapping = {}
    # Go through each BIGG metabolites file line (which is a tab-separated file)
    # and retrieve the BIGG ID and the name (if there is a name for the given BIGG
    # ID)
    for line in lines:
        bigg_id = line.split("\t")[1]
        # Exception to check if there is no name :O
        try:
            name = line.split("\t")[2].lower()
        except Exception:
            continue

        bigg_id_name_mapping[name] = bigg_id
        bigg_id_name_mapping[bigg_id] = bigg_id

    # Write the JSON in the given folder :D
    json_write(json_output_folder+"bigg_id_name_mapping.json",
               bigg_id_name_mapping)

def json_load(path: str) -> Dict[Any, Any]:
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: str ~ The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary

def parse_brenda_textfile(brenda_textfile_path: str, bigg_metabolites_json_folder: str,
                          json_output_path: str,json_output_path2: str) -> None:
    """Goes through a BRENDA database textfile and converts it into a machine-readable JSON.

    The JSON includes kcats for found organisms and substrates.
    As of 29/04/2019, the BRENDA database can be downloaded as textfile under
    https://www.brenda-enzymes.org/download_brenda_without_registration.php

    The BRENDA database is not in a completely standardized format, so that this functions
    contains many convoluted checks and circumventions of non-standardized data.

    kcats from mutated enzymes are excluded.

    Arguments
    ----------
    * brenda_textfile_path: str ~ The BRENDA database text file path
    * bigg_metabolites_json_folder: str ~ The folder in which the BIGG metabolites
      database is stored (it has to have the name 'bigg_id_name_mapping.json').
    * json_output_path: str ~ The path of the JSON that shall be created

    Output
    ----------
    * A JSON containing the BRENDA textfile kcat data in a machine-readable format:
    <pre>
        {
            "$EC_NUMBER": {
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
    # Standardize output folder
    bigg_metabolites_json_folder = standardize_folder(
        bigg_metabolites_json_folder)

    # Load BIGG ID <-> metabolite name mapping :D
    bigg_id_name_mapping: Dict[str, str] = json_load(
        bigg_metabolites_json_folder+"bigg_id_name_mapping.json")

    # Load BRENDA textfile as list of strings without newlines :D
    with open(brenda_textfile_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    lines = [x.replace("\n", "") for x in lines]

    # Go through each line and collect the organism lines and kcat lines for each EC number
    in_turnover_numbers = False
    in_specific_activity = False#Mr.Mao
    in_organism_reference = False
    ec_number_kcat_lines_mapping: Dict[str, List[str]] = {}
    ec_number_sa_lines_mapping: Dict[str, List[str]] = {}#Mr.Mao
    ec_number_organsism_lines_mapping: Dict[str, List[str]] = {}
    current_ec_number = ""
    organism_lines: List[str] = []
    kcat_lines: List[str] = []
    sa_lines: List[str] = []#Mr.Mao
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("ID\t"):
            if current_ec_number != "":
                ec_number_organsism_lines_mapping[current_ec_number] = organism_lines
                ec_number_kcat_lines_mapping[current_ec_number] = kcat_lines
                ec_number_sa_lines_mapping[current_ec_number] = sa_lines
            current_ec_number = line.replace("ID\t", "").replace(" ()", "")
            organism_lines = []
            kcat_lines = []
            sa_lines = []

        if len(line) == 0:
            in_turnover_numbers = False
            in_organism_reference = False
        elif line.startswith("PROTEIN"):
            in_organism_reference = True
            i += 1
            line = lines[i]
        elif line.startswith("TURNOVER_NUMBER"):
            in_turnover_numbers = True
            i += 1
            line = lines[i]
            
        ##Mr.MAO    
        elif line.startswith("SPECIFIC_ACTIVITY"):
            in_specific_activity = True
            i += 1
            line = lines[i]
            
        if in_organism_reference:
            if line.startswith("PR"):
                organism_lines.append("")
            if len(organism_lines[-1]) > 0:
                organism_lines[-1] += " "
            organism_lines[-1] += " " + line

        elif in_turnover_numbers:
            if line.startswith("TN"):
                kcat_lines.append("")
            if len(kcat_lines[-1]) > 0:
                kcat_lines[-1] += " "
            kcat_lines[-1] += line
        
        #SPECIFIC_ACTIVITY
        ##Mr.MAO 
        elif in_specific_activity:
            if line.startswith("SA"):
                sa_lines.append("")
            if len(sa_lines[-1]) > 0:
                sa_lines[-1] += " "
            sa_lines[-1] += line
            
        if len(line) == 0:
            in_turnover_numbers = False
            in_organism_reference = False
            in_specific_activity = False#Mr.Mao
            
        i += 1

    # Create the BRENDA database dictionary using the collected kcat and organism lines
    # of each EC number :D
    ec_numbers = list(ec_number_kcat_lines_mapping.keys())
    brenda_kcat_database: Dict[str, Any] = {}
    brenda_sa_database: Dict[str, Any] = {}#Mr.Mao
    for ec_number in ec_numbers:
        if "(transferred to " in ec_number:
            actual_ec_number = ec_number.split(" (transferred")[0]
            try:
                brenda_kcat_database[actual_ec_number] = {}
                brenda_kcat_database[actual_ec_number]["TRANSFER"] = \
                    ec_number.lower().replace("  ", " ").split(
                        "(transferred to ec")[1].replace(")", "").lstrip()
            except Exception:
                # Some transfers go to general subgroups instead of single EC numbers so that
                # no kcat database can be built from it D:
                print("WARNING: BRENDA text file line " + ec_number + " is not interpretable!")
            continue

        brenda_kcat_database[ec_number] = {}
        brenda_sa_database[ec_number] = {}
        
        reference_number_organism_mapping = {}
        organism_lines = ec_number_organsism_lines_mapping[ec_number]
        for organism_line in organism_lines:
            reference_number = organism_line.split("#")[1]
            organism_line_split_first_part = organism_line.split("# ")[1]
            organism_line_split = organism_line_split_first_part.split(" ")
            organism_line_split = [
                x for x in organism_line_split if len(x) > 0]

            end = 1
            for part in organism_line_split:
                # Some organism names contain their SwissProt or UniProt ID,
                # since we don't nned them they are excluded
                if ("swissprot" in part.lower()) or \
                    (part.lower() == "and") or \
                    ("uniprot" in part.lower()) or \
                    ("genbank" in part.lower()) or \
                        ("trembl" in part.lower()):
                    end -= 2
                    break

                if ("<" in part) or ("(" in part):
                    end -= 1
                    break

                end += 1
            organism_name = " ".join(organism_line_split[:end])
            reference_number_organism_mapping[reference_number] = organism_name

        kcat_lines = ec_number_kcat_lines_mapping[ec_number]
        for kcat_line in kcat_lines:
            kcat_line = kcat_line
            # Exclude kcats of mutated/changed proteins since
            # they may not have a biological relevance
            if ("mutant" in kcat_line.lower()) or ("mutated" in kcat_line.lower()):
                continue
            reference_number = kcat_line.split("#")[1].split(",")[0]
            organism = reference_number_organism_mapping[reference_number]
            kcat_str = "".join(kcat_line.split("#")[2]).split("{")[
                0].lstrip().rstrip()
            kcat = max([float(x) for x in kcat_str.split("-") if len(x) > 0])
            substrate = "".join(kcat_line.split("{")[1]).split("}")[0]

            substrate = substrate.lower()
            if substrate in bigg_id_name_mapping.keys():
                substrate = bigg_id_name_mapping[substrate]
            else:
                substrate = "REST"

            if substrate not in brenda_kcat_database[ec_number].keys():
                brenda_kcat_database[ec_number][substrate] = {}
            if organism not in brenda_kcat_database[ec_number][substrate].keys():
                brenda_kcat_database[ec_number][substrate][organism] = []
            brenda_kcat_database[ec_number][substrate][organism].append(kcat)
        
        #Mr.Mao
        sa_lines = ec_number_sa_lines_mapping[ec_number]
        for sa_line in sa_lines:
            sa_line = sa_line
            # Exclude kcats of mutated/changed proteins since
            # they may not have a biological relevance
            if ("mutant" in sa_line.lower()) or ("mutated" in sa_line.lower()):
                continue
            reference_number = sa_line.split("#")[1].split(",")[0]
            organism = reference_number_organism_mapping[reference_number]
            if re.search('\(',sa_line):
                sa_str = "".join(sa_line.split("#")[2].split("(")[
                0]).split("{")[0].lstrip().rstrip()
            elif re.search('<',sa_line):
                sa_str = "".join(sa_line.split("#")[2].split("<")[
                0]).split("{")[0].lstrip().rstrip()                
            else:
                sa_str = "".join(sa_line.split("#")[2]).split("{")[
                0].lstrip().rstrip()
            if re.search('e-',sa_str):
                sa = max([float(x) for x in sa_str.split(" ") if len(x) > 0]) 
            else:
                sa = max([float(x) for x in sa_str.split("-") if len(x) > 0])#为什么是-？是不是有1-2的值？，那如果是-e呢？
            if organism not in brenda_sa_database[ec_number].keys():
                brenda_sa_database[ec_number][organism] = []
            brenda_sa_database[ec_number][organism].append(sa)
    # Write final BRENDA kcat database :D
    json_write(json_output_path, brenda_kcat_database)
    json_write(json_output_path2, brenda_sa_database)
    
def is_fitting_ec_numbers(ec_number_one: str, ec_number_two: str, wildcard_level: int) -> bool:
    """Check whether the EC numbers are the same under the used wildcard level.

    Arguments
    ----------
    * ec_number_one: str ~ The first given EC number.
    * ec_number_two: str ~ The second given EC number.
    * wildcard_level: int ~ The wildcard level.
    """
    if wildcard_level == 0:
        ec_number_one_full_numbers = ec_number_one.split(".")
        ec_number_two_full_numbers = ec_number_two.split(".")
    else:
        ec_number_one_full_numbers = ec_number_one.split(".")[:-wildcard_level]
        ec_number_two_full_numbers = ec_number_two.split(".")[:-wildcard_level]

    if ec_number_one_full_numbers == ec_number_two_full_numbers:
        return True
    else:
        return False
    
def _get_transfer_ec_number_entry(ec_number_entry_key: str, brenda_kcat_database_original: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Returns the new EC number to which the given EC number was transferred.

    This is indicated in the given dictionary by the 'TRANSFER' key.
    Especially since the EC class 7's (translocases) introduction in 2018, many EC numbers are being transferred to new ones.

    Arguments
    ----------
    *ec_number_entry_key: str ~ The EC number for which its newly assigned one shall be searched.
    *brenda_kcat_database_original: Dict[str, Dict[str, Any]] ~ The BRENDA database dictionary with TRANSFER data.

    Output
    ----------
    Either...
    * {'ERROR': None} if the transferred new EC number is transferred to the old one (a futile cycle :O), or if the new EC number
      is not in BRENDA's database :O, or if the new EC number does not have any kcat entries :O
    * the dictionary containing the substrates, organisms and associated kcats of the new EC number
    """
    ec_number_entry = brenda_kcat_database_original[ec_number_entry_key]
    while "TRANSFER" in ec_number_entry.keys():
        new_ec_number = ec_number_entry["TRANSFER"]
        if new_ec_number == ec_number_entry_key:
            return {"ERROR": None}
        if new_ec_number not in brenda_kcat_database_original.keys():
            return {"ERROR": None}
        ec_number_entry = brenda_kcat_database_original[new_ec_number]

    if ec_number_entry == {}:
        return {"ERROR": None}

    return copy.deepcopy(ec_number_entry)

def parse_brenda_json_for_model(sbml_path: str, brenda_json_path: str, output_json_path: str) -> None:
    """Reads out a BRENDA JSON file created with parse_brenda_textfile and creates a model-specific JSON.

    Arguments
    ----------
    * sbml_path: str ~ The path of the SBML model of which a specific BRENDA JSON kcat database
      shall be created
    * brenda_json_path: str ~ The full path to the BRENDA JSON created with parse_brenda_textfile.
    * output_json_path: str ~ The full path to the newly created JSON.

    Output
    ----------
    A JSON in the given folder and the name 'kcat_database_brenda.json', and with the following structure:
    <pre>
    {
        '$EC_NUMBER': {
            '$BIGG_ID_METABOLITE': {
                '$ORGANISM': [
                    kcat_list: float
                ],
                (...)
            },
            (...)
        },
        (...)
    }
    </pre>
    """
    #model: cobra.Model = cobra.io.read_sbml_model(sbml_path)
    if re.search('\.xml',sbml_path):
        model = cobra.io.read_sbml_model(sbml_path)
    elif re.search('\.json',sbml_path):
        model = cobra.io.json.load_json_model(sbml_path)
    # Get EC numbers of the model's reactions
    ec_numbers_of_model: List[str] = []
    for reaction in model.reactions:
        if "ec-code" not in reaction.annotation.keys():
            continue

        ec_numbers_of_reaction = reaction.annotation["ec-code"]
        if type(ec_numbers_of_reaction) is str:
            ec_numbers_of_reaction = [ec_numbers_of_reaction]
        ec_numbers_of_model += ec_numbers_of_reaction
    ec_numbers_of_model = list(set(ec_numbers_of_model))

    # Get EC number entries for each EC number of the model
    brenda_kcat_database_original = json_load(brenda_json_path)
    brenda_kcat_database_for_model = {}
    for ec_number in ec_numbers_of_model:
        entry_error = False
        if ec_number in brenda_kcat_database_original.keys():
            ec_number_entry = _get_transfer_ec_number_entry(
                ec_number, brenda_kcat_database_original)
            if "ERROR" in ec_number_entry.keys():
                entry_error = True
            else:
                ec_number_entry["WILDCARD"] = False
                brenda_kcat_database_for_model[ec_number] = ec_number_entry

        if (ec_number not in brenda_kcat_database_original.keys()) or entry_error:
            eligible_ec_number_entries: List[Dict[str, Any]] = []
            for wildcard_level in range(1, 5):
                for database_ec_number in list(brenda_kcat_database_original.keys()):
                    if is_fitting_ec_numbers(ec_number, database_ec_number, wildcard_level):
                        database_ec_number_entry = _get_transfer_ec_number_entry(
                            database_ec_number, brenda_kcat_database_original)
                        if "ERROR" not in database_ec_number_entry.keys():
                            eligible_ec_number_entries.append(
                                database_ec_number_entry)
                if len(eligible_ec_number_entries) > 0:
                    break
            ec_number_entry = {}
            for eligible_ec_number_entry in eligible_ec_number_entries:
                for metabolite_key in eligible_ec_number_entry.keys():
                    metabolite_entry = eligible_ec_number_entry[metabolite_key]
                    if metabolite_key not in ec_number_entry.keys():
                        ec_number_entry[metabolite_key] = metabolite_entry
                    else:
                        ec_number_entry[metabolite_key] = {
                            **ec_number_entry[metabolite_key], **metabolite_entry}
            ec_number_entry["WILDCARD"] = True
            brenda_kcat_database_for_model[ec_number] = ec_number_entry

    json_write(output_json_path, brenda_kcat_database_for_model)
    
# SCRIPT-WIDE CONSTANTS
# URL for SABIO-RK's kcat REST API
QUERY_URL = "http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv"
# Time in seconds to wait between SABIO-RK API calls
WAIT_TIME = 1.5
# Constant unit multipliers since SABIO-RK contains kcats in different units :O
# The multiplication numbers show with which number a kcat of the corresponding unit
# has to be multiplied in order to get the kcat in s^(-1), the kcat unit that is
# used all over AutoPACMEN.
UNIT_MULTIPLIER: Dict[str, float] = {
    "s^(-1)": 1.0,
    "min^(-1)": 1/60,
    "h^(-1)": 1/(60*60),
}
    
def ensure_folder_existence(folder: str) -> None:
    """Checks if the given folder exists. If not, the folder is created.

    Argument
    ----------
    * folder: str ~ The folder whose existence shall be enforced.
    """
    if os.path.isdir(folder):
        return
    os.makedirs(folder)

def get_files(path: str) -> List[str]:
    """Returns the names of the files in the given folder as a list of strings.

    Arguments
    ----------
    * path: str ~ The path to the folder of which the file names shall be returned
    """
    files: List[str] = []
    for (_, _, filenames) in os.walk(path):
        files.extend(filenames)
    return files
    
def _add_wildcard_to_ec_number(ec_number: str, level: int) -> str:
    """Adds asterisk wildcards to the given EC number.

    Input
    ----------
    * ec_number: str ~ The EC number which shall be wildcarded
    * level: int ~ The wildcard level. It has to be at least 0 and
      maximally 5


    Example
    ----------
    <pre>
    >>> _add_wildcard_to_ec_number("1.1.1.1", 1)
    1.1.1.*
    </pre>

    Output
    ----------
    * The wildcarded EC number as string
    """
    ec_number_list = ec_number.split(".")[::-1]
    #print(ec_number_list,level)
    i = 1
    while i <= level:
        if i < level:
        #对于枯草，有些酶ec长度不对
        #if i < level and len(ec_number_list)==4:#
            ec_number_list[i] = "*"
        i += 1
    wildcarded_ec_number = ".".join(ec_number_list[::-1])
    return wildcarded_ec_number

def _extract_kcat_lines(result: csv.DictReader) -> List[str]:
    """Converts a SABIO-RK csv.DictReader line into a list of strings with the content of these lines.

    Input
    ----------
    * result: csv.DictReader ~ The csv.DictReader instance of a SABIO-RK kcat saerch output line

    Output
    ----------
    A list of all kcat lines as strings
    """
    kcat_lines: List[Any] = []
    for row in result:
        if (row["parameter.type"] == "kcat") and (row["parameter.startValue"] != ""):
            kcat_lines.append(row)
    return kcat_lines

def _get_species_results(result: csv.DictReader) -> List[str]:
    """Returns the organism data from the given SABIO-RK API csv.DictReader

    Arguments
    ----------
    * result: csv.DictReader ~ The csv.DictReader instance of a SABIO-RK kcat saerch output line

    Output
    ----------
    The organism data from the SABIO-RK API result
    """
    species_results: List[str] = []
    for row in result:
        species = row["Organism"]
        if species not in species_results:
            species_results.append(species)
    return species_results

def sabio_rk_query_with_string(query_string: str) -> str:
    """Call SABIO-RK API with the given query string.

    The query string is the actual term that is searched

    Arguments
    ----------
    * query_string: str ~ The query string

    Output
    ----------
    SABIO-RK's request result as str:

    If results could be found, it is a csv-like structure that can
    be futher processed e.g. with the Python standard library function csv.DictReader.
    The returned result fields are the EC number, the KEGG reaction ID, the organism
    and the 'parameter', i.e. the category of the given result. A possible category
    is 'kcat'.
    If no results could be found, it is the str "NO_RESULT".
    """
    # Build-up SABIO-RK query with a dictionary format that can be understood by the
    # requests library (this is also the way to retrieve potential huge amounts of
    # information as shown in SABIO-RK's API documentation).
    query = {"fields[]": ["ECNumber", "KeggReactionID", "Organism", "Parameter", "Substrate"], "q": query_string}

    # Send the request to SABIO-RK :D
    request = requests.post(QUERY_URL, params=query)

    # Error check whether the API call was successful or
    # not. 'Not successful' means that no search result with the given query
    # could be found. In this case, "NO_RESULT" is returned.
    try:
        request.raise_for_status()
    except Exception:
        #print("SABIO-RK API error with query:")
        #print(query_string)
        time.sleep(WAIT_TIME)
        return "NO_RESULT"

    # Wait time in order to not overload SABIO_RK's server
    time.sleep(WAIT_TIME)

    # Return successful search result in the given csv-like structure format.
    return request.text

def sabio_rk_query_with_query_dicts(query_dicts: List[Any]) -> str:
    """Performs a SABIO-RK query with given query dicts and converts them into a string

    Input
    ----------
    * query_dicts: List[Any] ~ A list of query dicts in the form of {id_name: ec_number, "Parametertype": "kcat", "EnzymeType": "wildtype"}

    Output
    ----------
    * As in sabio_rk_query_with_string()
    """
    i = 0
    # Add 'AND' for all variables of a query dict
    for i in range(len(query_dicts)):
        query_dicts[i] = " AND ".join([f"{k}:{v}" for k, v in query_dicts[i].items()])
    # Add 'OR' between each individual query dict
    query_string = "(" + ") OR (".join(query_dicts) + ")"
    # Add parentheses for URL conformity
    query_string = "(" + query_string + ")"

    # Perform the API call :D
    return sabio_rk_query_with_string(query_string)

def sabio_rk_query_get_csv_lines(query_dicts: List[Any]) -> Any:
    """Returns the result of a SABIo-RK API call with the given query dicts as csvReader lines

    Input
    ----------
    * query_dicts: List[Any] ~ A list of query dicts in the form of {id_name: ec_number, "Parametertype": "kcat", "EnzymeType": "wildtype"}

    Output
    ----------
    * The sabio_rk_query_with_string() result in a different form:
      Either "NO_RESULT" (a str) if no search esult was found,
      or a list of strings which includes the lines of the CSV that
      is returned by SABIO-RK
    """
    result = sabio_rk_query_with_query_dicts(query_dicts)
    if result == "NO_RESULT":
        return "NO_RESULT"
    else:
        return list(csv.DictReader(io.StringIO(result), delimiter="\t"))

def get_id_associated_kcats(searched_ids: List[str], id_type: str,
                            bigg_id_name_mapping_path: str, batch_size: int = 5) -> Dict[str, Any]:
    """Returns a dictionary with SABIO-RK kcat data for the given EC numbers or KEGG IDs.

    This function calls the SABIO-RK API.

    Input
    ----------
    * searched_ids: List[str] ~ The list of searched IDs
    * id_type: str ~ Must be either 'EC' or 'KEGG', depending on whether you are looking for kcats for EC numbers
      or KEGG IDs.
    * batch_size: int = 5 ~ The SABIO-RK API search batching number (i.e., with satch_size=5 five IDs are searched at once)

    Output
    ----------
    A dictionary with the following content:
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
    # Set-up the cache if it does not exist yet \o/
    cache_basepath = "./_cache/sabio_rk_total/"
    ensure_folder_existence("./_cache/")
    ensure_folder_existence(cache_basepath)
    cache_files = get_files(cache_basepath)
    # Load the given BIGG ID<->metabolite common name mapping
    bigg_id_name_mapping = json_load(bigg_id_name_mapping_path)
    # In order to save search time, use the seat (i.e., a list where
    # every member occurs only once) of the given searched IDs
    searched_ids = list(set(searched_ids))
    
    # Set the given ID name to the name which SABIO-RK uses for them
    if id_type == "EC":
        id_name = "ECNumber"
    elif id_type == "KEGG":
        id_name = "KeggReactionID"

    # Depending on the wildcard level which is serched, either
    # the output or the wildcard output will be used as output
    # These central dictionaries will contain the ID<->kcat mapping
    output = {}
    wildcard_output = {}
    # We use batched searched in order to save search time :D
    batch_start = 0
    # Loop while not all IDs were searched \o/
    while batch_start < len(searched_ids):
        # Get the batch for the search :-)
        batch = searched_ids[batch_start: batch_start + batch_size]
        # The query dicts contain a list of dictionaries which contain
        # the data for a SABIO-RK search entry
        query_dicts: List[Dict[str, str]] = []
        # Go through each single EC number in the search bath
        for ec_number in batch:
            # Create the cache filename
            cache_filename = ec_number.replace(".", "_").replace("*", "W") + ".json"
            # If the EC number is already searched, i.e. it can be found in the cache,
            # take the results from there in order to save much search time :D
            if cache_filename in cache_files:
                cache_filepath = cache_basepath + cache_filename
                output[ec_number] = json_load(cache_filepath)
                #print(f"Loading {cache_filename}...")
            # Otherwise, create an actual SABIO-RK API search query
            else:
                query_dicts.append({id_name: ec_number, "Parametertype": "kcat", "EnzymeType": "wildtype"})
        # If not all of the searched IDs are present in the cache...
        if len(query_dicts) > 0:
            # ...use SABIO-RK's API :D
            #print(f"Performing query {query_dicts}...")
            result = sabio_rk_query_get_csv_lines(query_dicts)

            # If there was an error with the SABIO-RK result (i.e., no result found or an invalid given ID),
            # continue with the next batch
            if result == "NO_RESULT":
                batch_start += batch_size
                continue
        # ...otherwise set the query result to nothing
        else:
            result = []

        # Loop through every SABIO-RK API query call result :D
        temp_ec_numbers_found_in_search = []
        result = _extract_kcat_lines(result)
        for row in result:
            # Get the unit of the parameter
            unit = row["parameter.unit"]
            # If it is a weird unusable unit, do not use this result and continue with the next result \o/
            if unit not in list(UNIT_MULTIPLIER.keys()):  # e.g. (s^-1)*(mg^-1)
                continue

            # Get the serached ID
            ec_number = row[id_name]
            # Generate a lowercarse and semicolon seperated list of substrates
            substrates_names = row["Substrate"]
            substrates_list = [x.lower() for x in substrates_names.replace("+", "").split(";")]
            substrates_list = sorted(substrates_list)
            # Convert the substrates name list into a BIGG ID list (only works
            # if there is a name<->BIGG ID mapping present for each substrate)
            bigg_ig_substrates_list = []
            for substrate in substrates_list:
                if substrate in bigg_id_name_mapping.keys():
                    bigg_id = bigg_id_name_mapping[substrate]
                    bigg_ig_substrates_list.append(bigg_id)
                # If one of the substrates cannot be found, use the pseudometabolite "REST"
                # and break :O
                else:
                    bigg_ig_substrates_list = ["REST"]
                    break
            # Set the substrate list to a semicolon-connected string
            substrate = ";".join(bigg_ig_substrates_list)
            # Get the result's organism :D
            species = row["Organism"]
            # Get the kcat and set
            # it to 1/s for consistent behaviour :D
            raw_kcat = float(row["parameter.startValue"])  # Without unit correction
            kcat = raw_kcat * UNIT_MULTIPLIER[unit]  # With unit correction 🎉

            # Add the result to the output for the given EC number, sustrate and species
            if ec_number not in output.keys():
                output[ec_number] = {}
            if substrate not in output[ec_number].keys():
                output[ec_number][substrate] = {}
            if species not in output[ec_number][substrate].keys():
                output[ec_number][substrate][species] = []
            output[ec_number][substrate][species].append(kcat)

            # Since we found a result, add the EC number :D
            temp_ec_numbers_found_in_search.append(ec_number)

        # Create cache files for all newly found EC numbers which were not present
        # in the cache
        temp_ec_numbers_found_in_search = list(set(temp_ec_numbers_found_in_search))
        for ec_number in temp_ec_numbers_found_in_search:
            cache_filename = ec_number.replace(".", "_") + ".json"
            if cache_filename not in cache_files:
                json_write(cache_basepath + cache_filename, output[ec_number])

        # Get all wildcarded searched EC numbers...
        wildcarded_searched_ec_numbers = [x for x in batch if "*" in x]
        # ...and loop through them in order to create a result for the EC numbers
        # which fit into the wildcard (i.e 1.1.1.123 in 1.1.1.*) :D
        for wildcarded_ec_number in wildcarded_searched_ec_numbers:
            # Ste the cache name for the wildcarded EC number
            cache_filename = wildcarded_ec_number.replace(".", "_").replace("*", "W") + ".json"
            # If the wildcarded EC number cannot be found in the cache, search for
            # fitting EC numbers, and combine their entries into a huge entry for the
            # wildcarded EC number
            if cache_filename not in cache_files:
                fitting_ec_numbers = []
                for found_ec_number in temp_ec_numbers_found_in_search:
                    if is_fitting_ec_numbers(wildcarded_ec_number, found_ec_number, wildcarded_ec_number.count("*")):
                        fitting_ec_numbers.append(found_ec_number)

                # Combine the EC number entries of fitting EC numbers :D
                wildcarded_ec_number_dict: Dict[str, Any] = {}
                for fitting_ec_number in fitting_ec_numbers:
                    fitting_ec_number_result = output[fitting_ec_number]
                    for metabolite_key in fitting_ec_number_result.keys():
                        if metabolite_key not in wildcarded_ec_number_dict.keys():
                            wildcarded_ec_number_dict[metabolite_key] = fitting_ec_number_result[metabolite_key]
                        else:
                            for organism_key in fitting_ec_number_result[metabolite_key].keys():
                                if organism_key not in wildcarded_ec_number_dict[metabolite_key].keys():
                                    wildcarded_ec_number_dict[metabolite_key][organism_key] =\
                                        copy.deepcopy(fitting_ec_number_result[metabolite_key][organism_key])
                                else:
                                    wildcarded_ec_number_dict[metabolite_key][organism_key] +=\
                                        copy.deepcopy(fitting_ec_number_result[metabolite_key][organism_key])
                                wildcarded_ec_number_dict[metabolite_key][organism_key] =\
                                    list(set(wildcarded_ec_number_dict[metabolite_key][organism_key]))
                # Create cache files for the searched wildcarded EC numbers \o/
                if wildcarded_ec_number_dict != {}:
                    json_write(cache_basepath + cache_filename, wildcarded_ec_number_dict)
                    wildcard_output[wildcarded_ec_number] = wildcarded_ec_number_dict
            # If the wildcarded EC number is in the cache, load the cache file :D
            else:
                wildcard_output[wildcarded_ec_number] = json_load(cache_basepath + cache_filename)
                #print(f"Loading {cache_filename}...")

        # Continue with the next searched ID batch :D
        batch_start += batch_size

    # If the wildcard level is greater than 0, set the wildcard output as output
    if len(wildcard_output.keys()) > 0:
        output = wildcard_output

    return output

def get_ec_number_kcats_wildcard_search(ec_numbers: List[str],
                                        bigg_id_name_mapping_path: str,
                                        batch_size: int = 5) -> Dict[str, Any]:
    """Returns EC number-dependent kcats using an incremental wildcard level until kcarts were found for all ECs.

    Arguments
    ----------
    *ec_numbers: List[str] ~ The list of checked EC numbers.
    *bigg_id_name_mapping_path: str ~ The full path to a BIGG ID<->Name mapping
    *batch_size: int = 5 ~ How many EC numbers shall be looked up in parallel in SABIO-RK. Do not set this too high!

    Output
    ----------
    A dictionary with the following content:
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
    #print("Starting EC numbers kcat search in SABIO-RK...")

    # We will look-up all EC numbers in SABIO-RK. If an EC number does not have an entry in
    # SABIO-RK, the wildcard level will become higher, e.g. 1.1.1.1 would become 1.1.1.*, and
    # if no entry is found in SABIO-EK, it would become 1.1.*.*, and so on
    # We start with no wildcards :D
    wildcard_level = 0
    # A list containing all found EC numbers, which won't be looked up using a higher wildcard level
    all_found_ec_numbers: List[str] = []
    # The central returned dictionary which will contain all combined kcat entries, divided
    # for organisms and metabolites
    ec_number_kcat_mapping: Dict[str, Any] = {}
    # Since an EC number has a maximum of 5 numbers and there is at least one EC number entry for
    # all EC number major category (e.g. 1.*.*.*), wildcards levels from 0 to 4 are reasonable,
    # and we loop throigh them :D
    for wildcard_level in range(5):
        # Get the list of all EC numbers which we want to search, i.e. all EC numbers which
        # are not already found.
        ec_numbers_to_analyze = list(set(ec_numbers) ^ set(all_found_ec_numbers))  # Difference
        # Add the current wildcard level to the searched EC numbers
        searched_ec_numbers = [_add_wildcard_to_ec_number(x, wildcard_level) for x in ec_numbers_to_analyze]
        # If no searched EC numbers are left with the current wildcard level, quit the for loop since we
        # are done :D
        if searched_ec_numbers == []:
            break
        # If there are searched EC numbers left, get the EC-number associated kcat entries
        #print(f"Wildcard level {wildcard_level}...")
        #print(searched_ec_numbers)
        if wildcard_level < 3:
            # With a low wildcard level, the default batch size of 5 is acceptable and
            # helps to search quicker
            resulting_ec_number_kcat_mapping = get_id_associated_kcats(searched_ec_numbers, "EC",
                                                                       bigg_id_name_mapping_path)
        else:
            # With a high wildcard level, no batching of searched EC numbers should be done since
            # the high number of results would make the search much slower D:
            searched_ec_numbers = list(set(searched_ec_numbers))
            if '1.*.*.*' in searched_ec_numbers:
                searched_ec_numbers.remove('1.*.*.*')
            if '3.*.*.*' in searched_ec_numbers:
                searched_ec_numbers.remove('3.*.*.*')
            #print(searched_ec_numbers)
            resulting_ec_number_kcat_mapping = get_id_associated_kcats(searched_ec_numbers, "EC",
                                                                       bigg_id_name_mapping_path,
                                                                       batch_size=1)

        # In the following loops, the resulting dictionaries for wildcarded EC numbers are mapped to
        # the searched EC numbers, e.g. the if 1.1.1.1233 and 1.1.1.2143 are searched and no results
        # were found with wildcard level 0, the results from wildcard level 1.1.1.* are given
        # to them with the entry that it comes from a wildcard
        resulting_found_ec_numbers = resulting_ec_number_kcat_mapping.keys()
        temp_all_found_ec_numbers: List[str] = []
        for ec_number in ec_numbers_to_analyze:
            if ec_number in all_found_ec_numbers:
                continue
            for found_ec_number in resulting_found_ec_numbers:
                if not is_fitting_ec_numbers(ec_number, found_ec_number, wildcard_level):
                    continue
                kcat = resulting_ec_number_kcat_mapping[found_ec_number]
                if ec_number not in list(ec_number_kcat_mapping.keys()):
                    ec_number_kcat_mapping[ec_number] = kcat
                    temp_all_found_ec_numbers += [ec_number]
                else:
                    ec_number_kcat_mapping[ec_number] = {**ec_number_kcat_mapping[ec_number], **kcat}
                if wildcard_level == 0:
                    ec_number_kcat_mapping[ec_number]["WILDCARD"] = False
                else:
                    ec_number_kcat_mapping[ec_number]["WILDCARD"] = True

        # Continue with the new list of found EC numbers and the next wildcard level :D
        all_found_ec_numbers += temp_all_found_ec_numbers
        wildcard_level += 1

    # Return the found major EC-number<->kcat mapping \o/
    return ec_number_kcat_mapping


"""
Exemplary usage:
kegg_ids = ["R00006", "R00286", "R00086"]
result = get_id_associated_kcats(kegg_ids, "KEGG")
print(result)

ec_numbers = ["1.1.1.213441", "2.7.1.124243", "1.234243.123213.123213", "2.12323.123213.12213"]
result = get_ec_number_kcats_wildcard_search(ec_numbers)
print(result)
"""

def parse_sabio_rk_for_model(model: cobra.Model, json_output_path: str, bigg_id_name_mapping_path: str) -> None:
    """Retrieves kcats from SABIO-RK for the given model and stores it in a JSON for the given model in the given path.

    Algorithm
    ----------
    Using the SABIO-RK REST API (as of 2019/30/04, it is explained under
    http://sabiork.h-its.org/layouts/content/docuRESTfulWeb/RESTWebserviceIntro.gsp),


    Arguments
    ----------
    * model: cobra.Model ~ The model for which kcats shall be retrieved from SABIO-RK.
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
    # GET LIST OF EC NUMBERS
    ec_numbers_list: List[str] = []
    for reaction in model.reactions:
        if "ec-code" not in reaction.annotation.keys():
            continue
        ec_codes = reaction.annotation["ec-code"]
        if type(ec_codes) is str:
            ec_codes = [ec_codes]
        ec_numbers_list += ec_codes
    ec_numbers_list = list(set(ec_numbers_list))

    # GET KCATS FOR EC NUMBERS
    ec_number_kcat_mapping = get_ec_number_kcats_wildcard_search(
        ec_numbers_list, bigg_id_name_mapping_path)

    json_write(json_output_path, ec_number_kcat_mapping)

#Mr.MAO
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
    ec_number_kcat_mapping = get_ec_number_kcats_wildcard_search(
        ec_numbers_list, bigg_id_name_mapping_path)

    json_write(json_output_path, ec_number_kcat_mapping)
    
def parse_sabio_rk_for_model_with_sbml(sbml_path: str, json_output_path: str, bigg_id_name_mapping_path: str) -> None:
    """See this module's parse_sabio_rk_for_model() documentation. This function uses an SBML path.

    Arguments
    ----------
    * sbml_path: str ~ The model's SBML path.
    * json_output_path: str ~ The path of the JSON that shall be created
    """
    # LOAD SBML MODEL
    #model: cobra.Model = cobra.io.read_sbml_model(sbml_path)
    if re.search('\.xml',sbml_path):
        model = cobra.io.read_sbml_model(sbml_path)
    elif re.search('\.json',sbml_path):
        model = cobra.io.json.load_json_model(sbml_path)

    parse_sabio_rk_for_model(model, json_output_path,
                             bigg_id_name_mapping_path)
    
def create_combined_kcat_database(sabio_rk_kcat_database_path: str, brenda_kcat_database_path: str, output_path: str) -> None:
    """Creates a combined JSON of the given SABIO-K and BRENDA kcat databases with non-wildcard entries only.

    Arguments
    ----------
    * sabio_rk_kcat_database_path: str ~ The path to the SABIO-RK kcat database JSON
    * brenda_kcat_database_path: str ~ The path to the BRENDA kcat database JSON
    * output_path: str ~ The outputh path (with filename) of the genreated combined kcat database JSON

    Output:
    A JSON with the following format:
    <pre>
    {
        '$EC_NUMBER': {
            '$BIGG_IDS_OF_SUBSTRATES': {
                '$ORGANISM': {
                    kcat: float
                },
                (...)
            },
            (...),
            'SOURCE': 'SABIO_RK' or 'BRENDA' or 'BRENDA and SABIO-RK',
            'WILDCARD': false
        },
        (...)
    }
    </pre>
    """
    # Load the two given databases as JSONs
    sabio_rk_database = json_load(sabio_rk_kcat_database_path)
    brenda_database = json_load(brenda_kcat_database_path)

    # Get all EC number keys (BRENDA contains all relevant EC numbers)
    ec_number_keys: List[str] = list(brenda_database.keys())
    # Set-up combined kcat database dictionary
    combined_database: Dict[str, Dict[str, Any]] = {}
    # Go through each EC number :D...
    for ec_number_key in ec_number_keys:
        # Get the wildcard status (i.e., found with a * wildcard?) and check if the EC number occurs anywhere...
        # ...for SABIO-RK
        if ec_number_key not in sabio_rk_database.keys():
            #print(f"WARNING: EC number {ec_number_key} could not be found in SABIO-RK, even with wildcards")
            #print("Possible reasons: The EC number format is invalid or there was an SABIO-RK API error")
            is_sabio_rk_from_wildcard = True  # Let the combined database ignore this entry
        else:
            is_sabio_rk_from_wildcard = sabio_rk_database[ec_number_key]["WILDCARD"]
        # ...and for BRENDA
        if ec_number_key not in brenda_database.keys():
            #print(f"WARNING: EC number {ec_number_key} could not be found in SABIO-RK, even with wildcards")
            #print("Possible reason: The EC number format is invalid")
            is_brenda_from_wildcard = True  # Let the combined database ignore this entry
        else:
            is_brenda_from_wildcard = brenda_database[ec_number_key]["WILDCARD"]

        # If both are from wildcards, ignore them :3
        if (is_sabio_rk_from_wildcard) and (is_brenda_from_wildcard):
            continue

        # Set-up dictionary for the EC number since at least one of the two databases
        # is not from a wildcarded search :D
        combined_database[ec_number_key] = {}
        # If both are not from wildcards, combine them :D...
        if (not is_sabio_rk_from_wildcard) and (not is_brenda_from_wildcard):
            # ...by reading their metabolites...
            sabio_rk_metabolite_keys = list(sabio_rk_database[ec_number_key].keys())
            brenda_metabolite_keys = list(brenda_database[ec_number_key].keys())
            metabolite_keys = list(set(sabio_rk_metabolite_keys + brenda_metabolite_keys))
            # ...going through them...
            for metabolite_key in metabolite_keys:
                # ...excluding the WILDCARD key...
                if metabolite_key == "WILDCARD":
                    continue
                # ...and adding the metabolites according to their presence in the databases :D
                is_metabolite_in_brenda: bool = metabolite_key in brenda_metabolite_keys
                is_metabolite_in_sabio_rk: bool = metabolite_key in sabio_rk_metabolite_keys
                if is_metabolite_in_brenda and is_metabolite_in_sabio_rk:
                    sabio_rk_entry = sabio_rk_database[ec_number_key][metabolite_key]
                    brenda_entry = brenda_database[ec_number_key][metabolite_key]
                    combined_database[ec_number_key][metabolite_key] = {**sabio_rk_entry, **brenda_entry}
                elif is_metabolite_in_brenda:
                    brenda_entry = brenda_database[ec_number_key][metabolite_key]
                    combined_database[ec_number_key][metabolite_key] = brenda_entry
                else:
                    sabio_rk_entry = sabio_rk_database[ec_number_key][metabolite_key]
                    combined_database[ec_number_key][metabolite_key] = sabio_rk_entry
            combined_database[ec_number_key]["WILDCARD"] = is_sabio_rk_from_wildcard
            combined_database[ec_number_key]["SOURCE"] = "BRENDA and SABIO-RK"
        # If only the SABIO-RK entry does not come from a wildcard, use it :D
        elif not is_sabio_rk_from_wildcard:
            combined_database[ec_number_key] = sabio_rk_database[ec_number_key]
            combined_database[ec_number_key]["WILDCARD"] = False
            combined_database[ec_number_key]["SOURCE"] = "SABIO-RK"
        # If only the BRENDA entry does not come from a wildcard, use it :-)
        elif not is_brenda_from_wildcard:
            combined_database[ec_number_key] = brenda_database[ec_number_key]
            combined_database[ec_number_key]["WILDCARD"] = False
            combined_database[ec_number_key]["SOURCE"] = "BRENDA"
    json_write(output_path, combined_database)
    
# SCRIPT-WIDE CONSTANTS
WAIT_TIME = .5  # Time to wait for each API call
NCBI_BATCH_SIZE = 20

def get_entrez_id_from_organism_full_name(organism_full_name):
    """Get organism's Entrez numeric identifier.

    This numeric identifier is neccessary for BLAST and NCBI TAXONOMY
    searches.
    This function uses Biopython functions. Returns BLAST-compatible ID as
    txid + NCBI ID + [ORGN].

    Arguments:
    >organism_kegg_id: str ~ The organism's full name, e.g. "Xanthomonas
     campestris pv. campesris B100"
    """
    # An e-mail has to be set, you may change it to yours if you want to
    # be notified if any problems occur.
    Entrez.email = "x@x.x"
    # Set the Entrez search to the NCBI TAXONOMY database.
    handle = Entrez.esearch(db="Taxonomy", term=organism_full_name)
    # Wait in order to not overload the NCBI's server
    time.sleep(WAIT_TIME)
    # Reformat the Entrez search result in order to extract the Entrez ID
    record = Entrez.read(handle)
    organism_ncbi_id = record["IdList"][0]
    # txid+NUMBER+[ORGN] is the form that is used for NCBI BLASTP searches to restrict a search
    # to an organism using the Entrez query constraint input.
    organism_ncbi_id = "txid"+organism_ncbi_id+"[ORGN]"
    # Return the retrieved ID :D
    return organism_ncbi_id

def get_taxonomy_from_organism_ncbi_id(organism_ncbi_id):
    """Get organism's taxonomy from NCBI Taxonomy using Biopython functions.

    The taxonomy is returned as list, starting with the nearest and
    ending with the highest taxonomic level above the organism.

    Arguments:
    >ncbi_organism_id: str ~ The organism's NCBI ID, e.g. retrieved by
     this module's "get_entrez_id_from_organism_full_name" function, in
     the format txid + NCBI ID + [ORGN]
    """
    Entrez.email = "x@x.x"
    handle = Entrez.efetch(db="Taxonomy", id=organism_ncbi_id, retmode="xml")
    records = Entrez.read(handle)
    taxonomy = records[0]["Lineage"].split(";")[::-1]
    taxonomy = [i.lstrip() for i in taxonomy]
    return taxonomy

def get_entrez_id_from_organism_full_name_batch(organism_full_names: List[str]) -> List[str]:
    """Retrieves the Entrez numeric ID of the given organisms.

    This numeric identifier is neccessary for BLAST and NCBI TAXONOMY
    searches.
    This function uses Biopython functions. Returns BLAST-compatible ID as
    txid + NCBI ID + [ORGN].

    Arguments:
    >organism_full_names: List[str] ~ A list of full names of organisms, e.g. "Xanthomonas
     campestris pv. campesris B100"
    """
    batch_start = 0
    organism_ncbi_ids_result: List[str] = []
    # Go through each organism :D
    while batch_start < len(organism_full_names):
        organism_full_names_slice = organism_full_names[batch_start:batch_start+NCBI_BATCH_SIZE]
        query_names = " OR ".join(organism_full_names_slice)
        # An e-mail has to be set, you may change it to yours if you want to
        # be notified if any problems occur.
        Entrez.email = "x@x.x"
        # Set the Entrez search to the NCBI TAXONOMY database.
        handle = Entrez.esearch(db="Taxonomy", term=query_names)
        # Wait in order to not overload the NCBI's server
        time.sleep(WAIT_TIME)
        # Reformat the Entrez search result in order to extract the Entrez ID
        record = Entrez.read(handle)
        organism_ncbi_ids = record["IdList"][::-1]
        # txid+NUMBER+[ORGN] is the form that is used for NCBI BLASTP searches to restrict a search
        # to an organism using the Entrez query constraint input.
        organism_ncbi_ids_result += ["txid"+x +
                                     "[ORGN]" for x in organism_ncbi_ids]

        batch_start += NCBI_BATCH_SIZE
        time.sleep(WAIT_TIME)
    # Return the retrieved IDs :D
    return organism_ncbi_ids_result

def get_taxonomy_from_organism_ncbi_id_batch(organism_ncbi_ids: List[str]) -> Dict[str, List[str]]:
    """Get the taxonomy from NCBI Taxonomy of the given organisms using Biopython functions.

    The taxonomy is returned as Dictionary (Dict[str, List[str]) for each organism,
    where each value is a string list starting with the nearest and
    ending with the highest taxonomic level above the organism.

    Arguments:
    >organism_ncbi_ids: List[str] ~ The list of the NCBI IDs of the organisms,
     e.g. retrieved by this module's "get_entrez_id_from_organism_full_name"
     function, in the format txid + NCBI ID + [ORGN]
    """
    taxonomies: Dict[str, List[str]] = {}
    batch_start = 0
    while batch_start < len(organism_ncbi_ids):
        organism_ncbi_ids_slice = organism_ncbi_ids[batch_start:batch_start+NCBI_BATCH_SIZE]
        query_ids = " OR ".join(organism_ncbi_ids_slice)
        Entrez.email = "x@x.x"
        handle = Entrez.efetch(db="Taxonomy", id=query_ids, retmode="xml")
        records = Entrez.read(handle)
        for record in records:
            taxonomy = record["Lineage"].split(";")[::-1]
            taxonomy = [i.lstrip() for i in taxonomy]
            taxonomies[record["ScientificName"]] = taxonomy
        batch_start += NCBI_BATCH_SIZE
    return taxonomies


def most_taxonomic_similar(base_species: str, taxonomy_dict: Dict[str, List[str]]) -> Dict[str, int]:
    """Returns a dictionary with a score of taxonomic distance from the given organism.

    e.g. if base_species is "Escherichia coli" and taxonomy_dict is
    <pre>
    {
        "Escherichia coli": ["Escherichia", "Bacteria", "Organism"],
        "Pseudomonas": ["Pseudomonas", "Bacteria", "Organism"],
        "Homo sapiens": ["Homo", "Mammalia", "Animalia", "Organism"],
    }
    </pre>
    this function would return
    <pre>
    {
        "Escherichia coli": 0,
        "Pseudomonas": 1,
        "Homo sapiens": 2,
    }
    </pre>

    Arguments
    ----------
    * base_species: str ~ The species to which a relation is made.
    * taxonomy_dict: Dict[str, List[str]] ~ A dictionary with organism names as keys and
      their taxonomic levels (sorted from nearest to farthest) as string list.
    """
    base_taxonomy = taxonomy_dict[base_species]
    level: int = 0
    level_dict: Dict[str, int] = {}
    for taxonomic_level in base_taxonomy:
        level_dict[taxonomic_level] = level
        level += 1

    score_dict: Dict[str, int] = {}
    for species in taxonomy_dict.keys():
        for taxonomic_level in taxonomy_dict[species]:
            if taxonomic_level in list(level_dict.keys()):
                score_dict[species] = level_dict[taxonomic_level]
                break

    return score_dict


"""
# Example:
organism_ncbi_ids = get_entrez_id_from_organism_full_name_batch(["Escherichia coli", "Escherichia fergusonii", "Vibrio natriegens", "Mus musculus"])
taxonomies = get_taxonomy_from_organism_ncbi_id_batch(organism_ncbi_ids)
print(taxonomies)
print(most_taxonomic_similar("Escherichia coli", taxonomies))
"""
    
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

def _get_kcat_from_protein_kcat_database(searched_direction: str, reaction: cobra.Reaction, protein_kcat_database):
    """Returns the kcat from the given protein<->kcat database for the given reaction, if there is one.

    The kcat is minimum of the maximal kcats for each protein of the reaction which can be found in the database.

    Arguments:
    *searched_direction: str ~ The direction in which a kcat is searched
    *reaction: cobra.Reaction ~ The reaction for which the kcat is searched
    *protein_kcat_database ~ The protein<->kcat database

    Output:
    Either the protein database kcat of the reaction in the given direction, or math.nan if there is no
    kcat for this reaction.
    """
    # Get the reaction's gene names
    gene_reaction_rule = reaction.gene_reaction_rule
    gene_reaction_rule = gene_reaction_rule.replace(" or ", "\t")
    gene_reaction_rule = gene_reaction_rule.replace(" and ", "\t")
    gene_names = gene_reaction_rule.split("\t")

    # Get the maximal kcats for each gene name in the given reaction direction
    max_kcats = []
    for gene_name in gene_names:
        if gene_name not in protein_kcat_database.keys():
            continue
        try:
            kcat_direction = protein_kcat_database[gene_name]["direction"][reaction.id]
        except:
            pass
        else:
            max_kcat = max(protein_kcat_database[gene_name]["kcats"])

            if kcat_direction == searched_direction == "forward":
                max_kcats.append(max_kcat)
            else:
                max_kcats.append(max_kcat)

    # Get the minimal maximal kcat and return it :D
    if len(max_kcats) > 0:
        min_max_kcat = min(max_kcats)
    else:
        min_max_kcat = math.nan

    return min_max_kcat


def _get_kcat_list(searched_metabolites: List[str], complete_entry: Dict[str, Any], complete_entry_sa: Dict[str, Any], MW, organism: str,
                   searched_direction: str, reaction: cobra.Reaction,
                   protein_kcat_database) -> List[float]:
    """Returns a list of kcats for the given reaction, created with a taxonomic search and usingg the protein kcat database.

    Algorithm
    ----------
    The given kcat entries are changed so that the kcats are now ordered and associated for each organism and
    the searched metabolies. Using NCBI Taxonomy, the taxonomic distance of the kcat tnry organisms to the given
    organisms is calculated, and - until the highest taxonomic distance or the desired minimal kcat entry length is
    reached - the kcats from the taxonomically nearest organisms are added. If there is a protein database kcat,
    it is added to the list too.

    Arguments
    ----------
    * searched_metabolites: List[str] ~ The list of fitting metabolites for the reaction
    * complete_entry: Dict[str, Any] ~ The complete entry of kcats for this entry brenda和sabio数据
    * organism: str ~ The organism of the model which shall become protein-constraint-enhanced
    * searched_direction: str ~ The affected reaction's direction, e.g. 'forward'
    * reaction: cobra.Reaction ~ The affected reaction for which a kcat shall be calculated
    * protein_kcat_database ~ The protein<->kcat database

    Output
    ----------
    The list of kcats which results from the taxonomic search
    """
    # Create a dictionary with kcat entries for each organism, determined from the given original kcat entry
    species_kcat_mapping: Dict[str, List[float]] = {}
    for searched_metabolite in searched_metabolites:
        species_entries = complete_entry[searched_metabolite]
        for species in list(species_entries.keys()):
            if species not in species_kcat_mapping.keys():
                species_kcat_mapping[species] = []
            species_kcat_mapping[species] += species_entries[species]

    # Create the list of all species
    all_species = list(species_kcat_mapping.keys())
    if organism not in all_species:
        all_species.append(organism)
        organism_added = True
    else:
        organism_added = False

    # Get the taxonomy of the kcat entry's organisms, either by calling
    # NCBI Taxonomy or by reading the cache of previously searched organisms
    cache_basepath = "./_cache/ncbi_taxonomy/"
    ensure_folder_existence("./_cache/")
    ensure_folder_existence(cache_basepath)
    cache_files = get_files(cache_basepath)
    species_to_search: List[str] = []
    taxonomy_dict_cache: Dict[str, List[str]] = {}
    for species in all_species:
        species=species.replace('/','_')
        cache_filename = species + "_taxonomy"
        #print(cache_filename)
        #cache_filename=cache_filename.replace('/','_')
        #print(cache_filename)
        if cache_filename in cache_files:
            cache_filepath = cache_basepath + cache_filename
            taxonomy_dict_cache[species] = pickle_load(cache_filepath)
        elif cache_filename+"_NA" in cache_files:
            cache_filepath = cache_basepath + cache_filename+"_NA"
            taxonomy_dict_cache[species] = pickle_load(cache_filepath)
        else:
            species_to_search.append(species)

    # If there are species which could be searched in NCBI Taxonomy, create a full
    # taxonomy dict which includes 'NOT FOUND' for all non-found organisms and the
    # taxonomies for each found organism.
    if len(species_to_search) > 0:
        #print(species_to_search)
        ncbi_ids = get_entrez_id_from_organism_full_name_batch(species_to_search)
        taxonomy_dict_search = get_taxonomy_from_organism_ncbi_id_batch(ncbi_ids)
        #print(taxonomy_dict_search)
        for searched_species in species_to_search:
            if searched_species not in taxonomy_dict_search.keys():
                taxonomy_dict_search[searched_species] = ["NOT FOUND"]

        for species in list(taxonomy_dict_search.keys()):
            cache_filename = species + "_taxonomy"
            if taxonomy_dict_search[species] == ["NOT FOUND"]:
                cache_filename += "_NA"
            cache_filepath = cache_basepath + cache_filename
            pickle_write(cache_filepath, taxonomy_dict_search[species])
        full_taxonomy_dict = {**taxonomy_dict_search, **taxonomy_dict_cache}
    else:
        full_taxonomy_dict = taxonomy_dict_cache

    # Process the taxonomies in order to find the taxonomic distances of the given organism to the organisms
    # which are present in the kcat entries.
    score_dict = most_taxonomic_similar(organism, full_taxonomy_dict)
    for species in full_taxonomy_dict.keys():
        if species not in score_dict:
            score_dict[species] = max(score_dict.values())+1
    # If we added the organism without kcat entries for it, we delete its distance
    # since there is no kcat which can be retrieved from it
    if organism_added:
        del(score_dict[organism])
    
    #按照物种近源性进行循环（循环终止条件kcat列表长度和物种距离）
    # Loop through the given organisms taxonomically and start with the lowest distance
    # Keep looping to higher taxonomic distances until the highest distance is reached
    # or the desired minimal number of kcat entries is reached
    minimal_distance = min(score_dict.values())
    maximal_distance = max(score_dict.values())
    current_distance = minimal_distance
    num_min_kcat_entries = 10#为什么是保留10条记录？
    kcat_list: List[float] = []
    species_list: List[str] = []  #Mr.Mao 
    kcat_extend='False'#Mr.Mao 
    ori_kcat: List[float] = []#Mr.Mao 
    #while (len(kcat_list) < num_min_kcat_entries) and (current_distance <= maximal_distance):
    #    for species in score_dict.keys():
    #        if species not in species_kcat_mapping.keys():  # e.g. soil bacterium -> bacterium
    #            continue
    #        if score_dict[species] == current_distance:
    #            kcat_list += species_kcat_mapping[species]
    #            species_list.append(species)
    #            kcat_extend='True'
     #       if species==organism:
    #            #print(species)
    #            ori_kcat=species_kcat_mapping[species]#Mr.Mao 
    #    current_distance += 1
    ###change code by Mr.Mao###
    #First,kcat
    if organism in species_kcat_mapping.keys():
        kcat_list=species_kcat_mapping[organism]
        ori_kcat=species_kcat_mapping[organism]
        kcat_extend='Database'
        #print('kcat',organism,kcat_list,ori_kcat)
    #Second,sa
    elif organism in complete_entry_sa.keys() and MW !='none':
        #print(complete_entry_sa[species],MW)
        kcat_list=[sa*MW/60 for sa in complete_entry_sa[organism]]####kg / mol  g/mmol kcat(s-1)
        ori_kcat=complete_entry_sa[organism]#Mr.Mao 
        kcat_extend='Brenda_SA'
        #print('sa',organism,kcat_list,ori_kcat)
    #Last,other species
    else:
        while (len(kcat_list) < num_min_kcat_entries) and (current_distance <= maximal_distance):
            for species in score_dict.keys():
                if species not in species_kcat_mapping.keys():  # e.g. soil bacterium -> bacterium
                    continue
                if score_dict[species] == current_distance and species!=organism:#Mr.Mao
                    kcat_list += species_kcat_mapping[species]
                    species_list.append(species)
                    kcat_extend='Other_species'
                #if species==organism:
                    #print(species)
                #    ori_kcat=species_kcat_mapping[species]#Mr.Mao 
            current_distance += 1
    ##########################
    # Get the protein database kcat for this reaction (if there is one, otherwise it returns math.nan)
    if protein_kcat_database != {}:
        protein_database_kcat = _get_kcat_from_protein_kcat_database(searched_direction, reaction, protein_kcat_database)

        # Add the protein database kcat if there is one, it will influence the resulting kcat since a mean is used
        if protein_database_kcat is not math.nan:
            kcat_list.append(protein_database_kcat)

    return [kcat_list,ori_kcat,species_list,kcat_extend]


def _get_kcat(searched_metabolites, complete_entry, complete_entry_sa, MW, organism: str, searched_direction: str, reaction: cobra.Reaction,
              protein_kcat_database, type_of_kcat_selection: str = "mean") -> float:
    """Returns the kcat for the reaction with the searched metabolites, the given organism and the given complete entry.

    Arguments
    ----------
    * searched_metabolites ~ The metabolites which were selected as substrates.
    * complete_entry ~ The complete entry of the eligible EC numbers for this reaction.
    * organism: str ~ The analyzed organism's name
    * searched_direction: str ~ The direction in which the
    * reaction: cobra.Reaction ~ The viewed reaction.
    * protein_kcat_database ~ The protein-dependent kcat database that may be added if one exists.
    * type_of_kcat_selection: str ~ The type of kcat selection. Can be "mean", "median" or "random". Is "mean" by default.

    Algorithm
    ----------
    _get_kcat_list is called with the given arguments, which returns a list of eligible kcats. If the length of the eligible
    kcat's list is shorter than 10, a new kcat search without metabolite constriants is performed in order to potentially get
    more kcats.
    With the final kcat list, the mean of the kcats is taken and returned as output.
    """
    # Get the list of all eligible kcats
    [kcat_list,ori_kcat,species_list,kcat_extend]= _get_kcat_list(searched_metabolites, complete_entry,complete_entry_sa,MW, organism, searched_direction, reaction, protein_kcat_database)#Mr.Mao
    
    ###change code by Mr.Mao###
    # If the list is shorter than 10, search without any metabolite constraint in order to potentially get more kcats
    #if len(kcat_list) < 10:
    #    searched_metabolites = ["ALL"]
    #    [kcat_list,ori_kcat,species_list,kcat_extend] = _get_kcat_list(searched_metabolites, complete_entry, organism, searched_direction, reaction, protein_kcat_database)#Mr.Mao
    if len(kcat_list) < 1:
        searched_metabolites = ["ALL"]
        [kcat_list,ori_kcat,species_list,kcat_extend] = _get_kcat_list(searched_metabolites, complete_entry, complete_entry_sa, MW, organism, searched_direction, reaction, protein_kcat_database)#Mr.Mao
    ##########################
    
    # Take the eligible kcats using the given selection method
    if type_of_kcat_selection == "mean":
        kcat = statistics.mean(kcat_list)
    elif type_of_kcat_selection == "median":
        kcat = statistics.median(kcat_list)
    elif type_of_kcat_selection == "random":
        kcat = random.choice(kcat_list)
    else:
        print("Wrong type_of_kcat_selection set! Must be 'mean', 'median' or 'random'.")
        sys.exit(-1)

    # Return the determined kcat :D 返回的是10个kcat  'mean', 'median' or 'random'
    return [kcat,kcat_list,ori_kcat,species_list,kcat_extend]


def _get_searched_metabolites(complete_entry, reaction_part_bigg_ids: List[str]) -> List[str]:
    """Returns which metabolites have a valid BIGG ID that ca be found in the complete entry.

    Arguments
    ----------
    * complete_entry ~ The complete kcat entry of all eligible EC numbers of the analyzed reaction
    * reaction_part_bigg_ids: List[str] ~ A list of BIGG IDs of either the products or educts of a reaction

    Output
    ----------
    A list of BIGG IDs which can be used as identifiers, or "ALL" if no fitting metabolite BIGG ID could be found
    """
    # Go through every metabolite in the complete entry
    eligible_metabolites = []
    for complete_entry_key in complete_entry.keys():
        # Check for identical names
        for reaction_part_bigg_id in reaction_part_bigg_ids:
            if reaction_part_bigg_id == complete_entry_key:
                eligible_metabolites.append(complete_entry_key)

        # Check for names in SABIO-RK style, e.g. "etoh_c;atp;h2o"
        num_found_metabolites = 0
        for reaction_part_bigg_id in reaction_part_bigg_ids:
            if (f";{reaction_part_bigg_id};") in (f";{complete_entry_key};"):
                num_found_metabolites += 1
        if num_found_metabolites == len(complete_entry_key.split(";")):
            eligible_metabolites.append(complete_entry_key)

    # Get the unique set of eligible metabolite IDs
    eligible_metabolites = list(set(eligible_metabolites))

    # If no eligible metabolite ID could be found, set "ALL" as pseudo-metabolite
    # indicating that all kcats are usable
    if len(eligible_metabolites) == 0:
        eligible_metabolites = ["ALL"]

    return eligible_metabolites

# PUBLIC FUNCTIONS
def get_reactions_kcat_mapping(sbml_path: str, project_folder: str, project_name: str,organism: str,kcat_database_path: str,
                                sa_database_path: str, reaction_mw_path: str,protein_kcat_database_path: str,
                               type_of_kcat_selection: str) -> None:#Mr.Mao
    """Returns a reaction<->kcat mapping for the given model :D

    The selection of kcats is depending on the affected metabolites of the reaction direction (one
    kcat is given for each the forward and reverse direction), and on the organism (the kcats
    from the taxonomically nearest organism is prefered).

    Arguments
    ----------
    *sbml_path: str ~ Te SBML path to the model
    *project_folder: str ~ The folder in which the model data files are sored
    *project_name: str ~ The name of the used project
    *organism: str ~ The organism's name
    *kcat_database_path: str ~ A path to an already created EC number<->kcats database
    *sa_database_path: str ~ A path to an already created EC number<->sas database#Mr.Mao
    *protein_kcat_database_path: str ~ A path to the custom protein<->kcat database
    *type_of_kcat_selection: str ~ Can be "mean", "median" or "random". Refers to the selection of found kcats of a reaction.
                                   Is "mean" by default.

    Output
    ----------
    A JSON in the given project folder with the name $project_name+'_reactions_kcat_mapping_combined.json' and
    the following structure:
    <pre>
    {
        "$REACTION_NAME": {
            "forward": $forward_kcat,
            "reverse": $reverse_kcat
        },
        (...)
    }
    </pre>
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)
    # Set the path for the output JSON
    basepath = project_folder + project_name
    # Load the combined, EC-number-dependent kcat database :D
    kcat_database = json_load(kcat_database_path)
    # Load the combined, EC-number-dependent sa database :D#Mr.Mao
    sa_database = json_load(sa_database_path)#Mr.Mao
    # If given, load the protein-dependent kcat database :D
    reaction_mw_mapping: Dict[str, float] = json_load(reaction_mw_path)#Mr.Mao
    if protein_kcat_database_path != "none":
        protein_kcat_database = json_load(protein_kcat_database_path)
    else:
        protein_kcat_database = {}

    # Load the given stoichiometric model
    #model = cobra.io.read_sbml_model(sbml_path)
    if re.search('\.xml',sbml_path):
        model = cobra.io.read_sbml_model(sbml_path)
    elif re.search('\.json',sbml_path):
        model = cobra.io.json.load_json_model(sbml_path)

    # Set-up dictionary which will be the content of the output JSON
    reactions_kcat_mapping: Dict[str, Dict[str, float]] = {}
    # Go through each reaction in order to assign kcats for it :D
    for reaction in model.reactions:
        # If no EC number is given in the reaction's annotations,
        # the protein-dependent database is read out in order to
        # find a kcat. This only works if at least one of the assigned
        # enzymes of the reaction's gene rule has a kcat in the
        # protein-dependent database.
        if "ec-code" not in reaction.annotation.keys():
            # 0 means that no kcat can be assigned
            forward_kcat: Any = 0
            reverse_kcat: Any = 0
            
            if protein_kcat_database != {}:
                # Retrieve the kcats from the protein-dependent database :D
                forward_kcat = _get_kcat_from_protein_kcat_database("forward", reaction, protein_kcat_database)
                reverse_kcat = _get_kcat_from_protein_kcat_database("reverse", reaction, protein_kcat_database)

            # If no kcat could be assigned, set the kcat to math.nan
            # which indicates this case
            if forward_kcat == 0.0:
                forward_kcat = math.nan
            if reverse_kcat == 0.0:
                reverse_kcat = math.nan

            # Add the retrieved forward and reverse kcats to the reaction<->kcat mapping dictionary :D
            reactions_kcat_mapping[reaction.id] = {}
            reactions_kcat_mapping[reaction.id]["forward"] = forward_kcat
            reactions_kcat_mapping[reaction.id]["reverse"] = reverse_kcat

            # Print the assigned kcats
            #_print_assigned_kcats(reaction.id, forward_kcat, reverse_kcat)
            continue

        # Retrieve the reaction's associated EC numbers
        reaction_ids = reaction.annotation["ec-code"]
        # If only one EC number is given, set the EC number string to
        # a list in order to make it work with the following code lines
        if type(reaction_ids) is str:
            reaction_ids = [reaction_ids]
        # Get all EC numbers which do not contain a - wildcard, such as
        # in 2.1.1.-
        # These wildcarded EC numbers are in general too permissive in order
        # to get useful kcats
        eligible_reaction_ids = [x for x in reaction_ids if "-" not in x]
        if len(eligible_reaction_ids) == 0:
            eligible_reaction_ids = [x for x in reaction_ids]

        # Create a 'complete entry' from all eligible (i.e., non-wildcarded)
        # EC numbers. This complete entry contains - for every organism
        # and substrate given in the EC number kcat entries - all kcats
        # of all eligible EC numbers. In addition, the pseudo-substrate
        # "ALL" is added which contains all organisms. "ALL" is used
        # later if no fitting substrate can be found.
        complete_entry: Dict[str, Any] = {}
        complete_entry["ALL"] = {}
        complete_entry_sa: Dict[str, Any] = {}#Mr.Mao
        # Go through each reaction ID :D
        
        ###change code by Mr.Mao###
        for reaction_id in eligible_reaction_ids:
            if reaction_id not in sa_database.keys():
                   continue
            reaction_id_entry2 = sa_database[reaction_id]
            for species_key in reaction_id_entry2.keys():
                # Add the metabolite to the complete entry if it does not already occur
                if species_key not in complete_entry_sa:
                    complete_entry_sa[species_key] = []
                complete_entry_sa[species_key] += reaction_id_entry2[species_key]
        ##########################   
        
        for reaction_id in eligible_reaction_ids:
            # If the EC number could not be found in the given EC number<->kcat
            # database, print it and proceed with the next eligible EC number            
            if reaction_id not in kcat_database.keys():
                #print(f"INFO: No entry for EC number {reaction_id}")
                #print("")
                continue
            # Otherwise, get the reaction ID entry from the given database :D
            reaction_id_entry = kcat_database[reaction_id]
            # Exclude all kcat entries which come from a wildcard search
            # with *
            if reaction_id_entry["WILDCARD"]:
                continue
            # Go trough each metabolite in the EC number<->kcat database entries
            for metabolite_key in reaction_id_entry.keys():
                # Ignore the keys which show additional information
                # about the nature of the kcat data
                if metabolite_key in ("WILDCARD", "SOURCE", "TRANSFER"):
                    continue
                # Add the metabolite to the complete entry if it does not already occur
                if metabolite_key not in complete_entry:
                    complete_entry[metabolite_key] = {}
                # Go throudh each species in the currently analyzed EC number
                for species_key in reaction_id_entry[metabolite_key]:
                    # Add the species to the metabolite entry if it does not already occur
                    if species_key not in complete_entry[metabolite_key]:
                        complete_entry[metabolite_key][species_key] = []
                    # ...and do the same for the pseudo-metabolite "ALL"
                    if species_key not in complete_entry["ALL"].keys():
                        complete_entry["ALL"][species_key] = []
                    # Add the list of kcats of the currently analyzed EC number to the current species
                    # and the current metabolite, and for "ALL"
                    complete_entry[metabolite_key][species_key] += reaction_id_entry[metabolite_key][species_key]
                    complete_entry["ALL"][species_key] += reaction_id_entry[metabolite_key][species_key]
        
        #complete_entry实际brenda和sabio两个数据，后续时对这个分析
        
        # If no entries with kcats could be found for any of the eligible EC numbers, continue with the next reaction.
        if complete_entry["ALL"] == {}:
            continue

        # Get the BIGG IDs of the educts and products uusing the SBML's BIGG ID annotation
        educt_bigg_ids: List[str] = []
        for reactant in reaction.reactants:
            if "bigg.metabolite" in reactant.annotation.keys():
                educt_bigg_ids.append(reactant.annotation["bigg.metabolite"])
        product_bigg_ids: List[str] = []
        for product in reaction.products:
            if "bigg.metabolite" in product.annotation.keys():
                product_bigg_ids.append(product.annotation["bigg.metabolite"])
        # If no bigg IDs could be found in the SBML, add the pseudo-metabolite "X"
        # which indicated that "ALL" should be used later.
        if len(educt_bigg_ids) == 0:
            educt_bigg_ids = ["X"]
        if len(product_bigg_ids) == 0:
            product_bigg_ids = ["X"]
        
        #Mr.Mao
        if reaction.id in reaction_mw_mapping.keys():
            MW=reaction_mw_mapping[reaction.id]
        else:
            MW='none'
        ###
        # Get the metabolites which are used in the subsequent forward kcat search
        searched_educts = _get_searched_metabolites(complete_entry, educt_bigg_ids)
        # Get the forward kcat depending on the educts and the organism
        [forward_kcat,forward_kcat_list,ori_forward_kcat,forward_species_list,kcat_extend] = _get_kcat(searched_educts, complete_entry,complete_entry_sa,MW, organism, "forward", reaction, protein_kcat_database, type_of_kcat_selection)#Mr.Mao

        # Get the metabolites which are used in the subsequent forward kcat search
        searched_products = _get_searched_metabolites(complete_entry, product_bigg_ids)
        # Get the reverse kcat depending on the products and the organism
        [reverse_kcat,reverse_kcat_list,ori_reverse_kcat,reverse_species_list,kcat_extend] = _get_kcat(searched_products, complete_entry, complete_entry_sa, MW, organism, "reverse", reaction, protein_kcat_database, type_of_kcat_selection)#Mr.Mao

        # Set the found out kcats in the reactions<->kcat mapping :D
        #print(reaction.id,forward_kcat,forward_kcat_list,ori_forward_kcat,reverse_kcat,reverse_kcat_list,ori_reverse_kcat)
        reactions_kcat_mapping[reaction.id] = {}
        reactions_kcat_mapping[reaction.id]["forward"] = forward_kcat
        reactions_kcat_mapping[reaction.id]["forward_kcat_list"] = forward_kcat_list#Mr.Mao
        reactions_kcat_mapping[reaction.id]["forward_ori_kcat"] = ori_forward_kcat#Mr.Mao
        reactions_kcat_mapping[reaction.id]["forward_species_list"] = forward_species_list#Mr.Mao
        reactions_kcat_mapping[reaction.id]["reverse"] = reverse_kcat
        reactions_kcat_mapping[reaction.id]["reverse_kcat_list"] = reverse_kcat_list#Mr.Mao
        reactions_kcat_mapping[reaction.id]["reverse_ori_kcat"] = ori_reverse_kcat#Mr.Mao
        reactions_kcat_mapping[reaction.id]["reverse_species_list"] = reverse_species_list#Mr.Mao

        reactions_kcat_mapping[reaction.id]["data_type"] =kcat_extend
        # display the found out kcats for this reaction \o/
        #_print_assigned_kcats(reaction.id, forward_kcat, reverse_kcat)

    # Export the kcat mapping results as JSON :D
    json_write(basepath+"_reactions_kcat_mapping_combined.json", reactions_kcat_mapping)

def get_protein_mass_mapping(model: cobra.Model, project_folder: str, project_name: str) -> None:
    """Returns a JSON with a mapping of protein IDs as keys, and as values the protein mass in kDa.

    The protein masses are calculated using the amino acid sequence from UniProt (retrieved using
    UniProt's REST API).

    Arguments
    ----------
    * model: cobra.Model ~ The model in the cobrapy format
    * project_folder: str ~ The folder in which the JSON shall be created
    * project_name: str ~ The beginning of the JSON's file name

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
    uniprot_id_protein_mass_mapping: Dict[str, float] = {}
    # The cache stored UniProt masses for already searched
    # UniProt IDs (each file in the cache folder has the name
    # of the corresponding UniProt ID). This prevents searching
    # UniProt for already found protein masses. :-)
    cache_basepath = "./_cache/uniprot/"
    ensure_folder_existence("./_cache/")
    ensure_folder_existence(cache_basepath)
    cache_files = get_files(cache_basepath)
    # Go through each batch of UniProt IDs (multiple UniProt IDs
    # are searched at once in order to save an amount of UniProt API calls)
    # and retrieve the amino acid sequences and using these sequences, their
    # masses.
    #print("Starting UniProt ID<->Protein mass search using UniProt API...")
    uniprot_ids = list(uniprot_id_protein_id_mapping.keys())
    batch_size = 5
    batch_start = 0
    while batch_start < len(uniprot_ids):
        # Create the batch with all UniProt IDs
        prebatch = uniprot_ids[batch_start:batch_start+batch_size]
        batch = []
        # Remove all IDs which are present in the cache (i.e.,
        # which were searched for already).
        # The cache consists of pickled protein mass floats, each
        # onein a file with the name of the associated protein.
        for uniprot_id in prebatch:
            if uniprot_id not in cache_files:
                batch.append(uniprot_id)
            else:
                cache_filepath = cache_basepath + uniprot_id
                uniprot_id_protein_mass_mapping[uniprot_id] = pickle_load(cache_filepath)
                #print(uniprot_id+":", uniprot_id_protein_mass_mapping[uniprot_id])

        # If all IDs could be found in the cache, continue with the next batch.
        if len(batch) == 0:
            batch_start += batch_size
            continue

        # Create the UniProt query for the batch
        # With 'OR', all given IDs are searched, and subsequently in this script,
        # the right associated masses are being picked.
        query = " OR ".join(batch)
        # uniprot_query_url = f"https://www.uniprot.org/uniprot/?query={query}&format=tab&columns=id,sequence"
        uniprot_query_url = f"https://rest.uniprot.org/uniprotkb/search?query=accession:{query}&format=tsv&fields=accession,sequence"
        #print(f"UniProt batch search for: {query}")

        # Call UniProt's API :-)
        uniprot_data = requests.get(uniprot_query_url).text.split("\n")
        # Wait in order to cool down their server :-)
        time.sleep(2.0)

        # Read out the API-returned lines
        for line in uniprot_data[1:]:
            if line == "":
                continue
            uniprot_id = line.split("\t")[0]
            sequence = line.split("\t")[1]
            # Get the protein mass using biopython's associated function for amino acid sequences
            try:
                mass = ProteinAnalysis(sequence, monoisotopic=False).molecular_weight()
            except ValueError:  # e.g. if an "X" is in a sequence
                continue
            uniprot_id_protein_mass_mapping[uniprot_id] = float(mass)

        # Create the pickled cache files for the searched protein masses
        for uniprot_id in batch:
            cache_filepath = cache_basepath + uniprot_id
            pickle_write(cache_filepath, uniprot_id_protein_mass_mapping[uniprot_id])

        # Continue with the next batch :D
        batch_start += batch_size

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


def get_protein_mass_mapping_with_sbml(sbml_path: str, project_folder: str, project_name: str) -> None:
    """This module's get_protein_mass_mapping() with SBML instead of a cobrapy module as argument.

    Arguments
    ----------
    * sbml_path: str ~ The path to the model's SBML
    * project_folder: str ~ The folder in which the JSON shall be created
    * project_name: str ~ The beginning of the JSON's file name
    """
    #model: cobra.Model = cobra.io.read_sbml_model(sbml_path)
    if re.search('\.xml',sbml_path):
        model = cobra.io.read_sbml_model(sbml_path)
    elif re.search('\.json',sbml_path):
        model = cobra.io.json.load_json_model(sbml_path)
    get_protein_mass_mapping(model, project_folder, project_name)