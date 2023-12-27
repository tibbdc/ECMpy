from ECMpy_function import *
import argparse

def model_main():
    parser = argparse.ArgumentParser(description='Get data from model!')
    parser.add_argument('-m','--model_file', type=str, default = None)
    parser.add_argument('-kcat','--user_provided_kcat_file', type=str, default = None)   
    parser.add_argument('-f','--user_provided_f', type=str, default = None) 
    parser.add_argument('-bigg','--bigg_file', type=str, default = "./data/bigg_models_metabolites.txt")
    parser.add_argument('-gene_abundance','--gene_abundance_file', type=str, default = None)    
    parser.add_argument('-tax_id','--taxonom_id', type=str, default = None)
    parser.add_argument('-org','--organism', type=str, default = None)
    parser.add_argument('-sigma','--sigma', type=str, default = None)  
    parser.add_argument('-ptot','--ptot', type=str, default = None)     
    parser.add_argument('-kcat_method','--kcat_method', type=str, default = None)
    parser.add_argument('-work_folder','--work_folder', type=str, default = "./analysis/get_kcat_mw")  
    parser.add_argument('-brenda','--brenda_file', type=str, default = "./data/brenda_2023_1.txt") 
    parser.add_argument('-uniprot','--uniprot_file', type=str, default = "./data/uniprot_data_accession_key.json") 
    parser.add_argument('-kcat_gap_fill','--kcat_gap_fill', type=str, default = None)
    parser.add_argument('-r_gap_fill','--reaction_gap_fill', type=str, default = None)    
    parser.add_argument('-ecGEM','--ecGEMoutfile', type=str, default='./model/ecGEM.json')
    args = parser.parse_args()
    print(args)
    model_file = args.model_file
    user_kcat_file = args.user_provided_kcat_file
    user_f = args.user_provided_f
    bigg_met_file =  args.bigg_file 
    ecModel_output_file = args.ecGEMoutfile

    # determine whether the data is suitable for constructing an enzyme-constrained model
    [sui_or_not,error_list] = Determine_suitable_ecGEM(model_file,bigg_met_file)
    
    if sui_or_not =='Yes':
        # kcat mw 
        if user_kcat_file == 'No':
            kcat_method=args.kcat_method
            work_folder = args.work_folder
            if kcat_method =='AutoPACMEN':
                print('It\'s time to get reaction kcat_MW using AutoPACMEN!')
                print('-----------------------------------')
                brenda_textfile_path = args.brenda_file
                uniprot_data_file = args.uniprot_file
                organism = args.organism
                kcat_gap_fill = args.kcat_gap_fill
                reaction_gap_fill = args.reaction_gap_fill
                autopacmen_folder = work_folder+'_by_AutoPACMEN/'
                create_file(autopacmen_folder)
                project_name = "model_%s"%kcat_gap_fill
                protein_kcat_database_path = "none"
                reaction_kcat_MW_file = get_reaction_kcatmw_onestop_by_AutoPACMEN(autopacmen_folder,model_file,bigg_met_file,brenda_textfile_path,project_name,uniprot_data_file,organism,protein_kcat_database_path,kcat_gap_fill,reaction_gap_fill)
            elif kcat_method == 'DLKcat':
                print('It\'s time to get reaction kcat_MW using DLKcat!')
                print('-----------------------------------')
                dlkcat_folder = work_folder+'_by_DLKcat/'
                create_file(dlkcat_folder)
                reaction_kcat_MW_file = get_reaction_kcatmw_onestop_by_DLKcat(dlkcat_folder,model_file)
            else:
                print('The method you provided is not supported, please check the spelling!')
        else:
            reaction_kcat_MW_file = user_kcat_file
        #f
        if user_f == 'No':
            print('It\'s time to calculate the enzyme mass fraction (f) based on protein abundance!')
            print('-----------------------------------')
            gene_abundance_file =  args.gene_abundance_file
            taxonom_id =  args.taxonom_id 
            #print(gene_abundance_file,taxonom_id)
            f=calculate_f_v2(model_file, gene_abundance_file,'abundance',taxonom_id)
        else:
            f=user_f
        #ecGEM
        sigma =  float(args.sigma) 
        ptot =  float(args.ptot)
        f = float(f)
        #print(type(f), type(sigma),type(ptot))
        lowerbound = 0   # Lowerbound  of enzyme concentration constraint. 
        upperbound = round(ptot * f * sigma, 3)#total enzyme
        print('It\'s time to construct an enzyme-constrained model!')
        print('-----------------------------------')
        trans_model2enz_json_model_split_isoenzyme(model_file, reaction_kcat_MW_file, f, ptot, sigma, lowerbound, upperbound, ecModel_output_file)

    else:
        print(error_list)
          
if __name__ == '__main__':
    model_main()
