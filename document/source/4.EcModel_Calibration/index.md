# 4.EcModel Calibration


```python
import cobra
import re 
import pandas as pd
#from script.ECMpy_function import *
import sys
sys.path.append(r'./script/')
from ECMpy_function import *
```



```python
#Initial parameters
f = 0.406 
ptot = 0.56 
sigma = 1
lowerbound = 0   
upperbound = round(ptot * f * sigma, 3) 
EC_max_file='./data/EC_kcat_max.json'# https://www.brenda-enzymes.org/brenda_download/file_download.php
method='AutoPACMEN'#DLKcat
use_substrate='EX_glc__D_e'
concentration=10
obj='BIOMASS_Ec_iML1515_core_75p37M'# CG_biomass_cgl_ATCC13032 EX_lys_L_e

#Originl ecmodel and result file
ecModel_file="./model/iML1515_irr_enz_constraint_test.json"
fluxes_infile_ori = './analysis/ECMpy_solution_BIOMASS_Ec_iML1515_core_75p37M_pfba.csv'
reaction_kcat_MW_file = "./analysis/get_kcat_mw_by_%s/reaction_kcat_MW.csv"%method
need_change_reaction_list=[]
changed_reaction_list=[]
round_num=1
reaction_kcat_mw = pd.read_csv(reaction_kcat_MW_file, index_col=0)

#ecmodel and result file
json_output_file = './model/iML1515_irr_enz_constraint_adj_test.json'
reaction_kcat_MW_outfile = './analysis/get_kcat_mw_by_%s/reaction_change_by_enzuse.csv'%method
```

## Calibration kcat according Enzyme usage 


```python
#Calibration
enz_model_obj=0
while (enz_model_obj<0.66 and round_num<50):#maxium growth rate
    [enz_model,reaction_kcat_mw,need_change_reaction_list, changed_reaction_list]=change_enz_model_by_enz_usage(ecModel_file,fluxes_infile_ori,EC_max_file,\
                                                    reaction_kcat_mw,need_change_reaction_list,changed_reaction_list,f, \
                                                    ptot, sigma, lowerbound, upperbound, json_output_file)
    enz_model.objective=obj
    #change original substrate in model
    [ori_obj_id,ori_substrate_id_list,ori_sub_concentration,ori_ATPM]=get_model_substrate_obj(enz_model)
    for eachsubid in ori_substrate_id_list:
        if re.search('_reverse',eachsubid):
            r_id_new=eachsubid.split('_reverse')[0]
            enz_model.reactions.get_by_id(eachsubid).bounds = (0, 0) 
            enz_model.reactions.get_by_id(r_id_new).bounds = (0, 0)  
        else:
            r_id_new=eachsubid+'_reverse'
            enz_model.reactions.get_by_id(eachsubid).bounds = (0, 0) 
            enz_model.reactions.get_by_id(r_id_new).bounds = (0, 0) 
            
    enz_model.reactions.get_by_id(use_substrate).bounds = (-concentration, 0)
    enz_model.reactions.get_by_id(use_substrate+'_reverse').bounds = (0, 0)

    enz_model_obj=enz_model.slim_optimize()
    print('Calibration round %s : '%round_num+str(enz_model_obj))
    round_num=round_num+1
    
reaction_kcat_mw.to_csv(reaction_kcat_MW_outfile)
```

    Need changing reaction: 
    PRFGS
    Changed reaction: 
    ['PRFGS']
    Calibration round 1 : 0.23173512444593233
    Need changing reaction: 
    KARA1_reverse
    Changed reaction: 
    ['PRFGS']
    Calibration round 2 : 0.23173512444593233
    Need changing reaction: 
    FBA_num1
    Changed reaction: 
    ['PRFGS', 'FBA_num1']
    Calibration round 3 : 0.24483790513370257
    Need changing reaction: 
    PFL_num3
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3']
    Calibration round 4 : 0.24723409719596082
    Need changing reaction: 
    3OAS121
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121']
    Calibration round 5 : 0.26291740157273025
    Need changing reaction: 
    3OAS161
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161']
    Calibration round 6 : 0.2806003508042417
    Need changing reaction: 
    3OAS141
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141']
    Calibration round 7 : 0.29935260943266906
    Need changing reaction: 
    GLUDy_reverse
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141']
    Calibration round 8 : 0.29935260943266906
    Need changing reaction: 
    PSERT
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT']
    Calibration round 9 : 0.3019931955364434
    Need changing reaction: 
    PAPSR2_num2
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2']
    Calibration round 10 : 0.3077850507857091
    Need changing reaction: 
    3OAR60
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60']
    Calibration round 11 : 0.3137641168450432
    Need changing reaction: 
    3OAR80
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80']
    Calibration round 12 : 0.3199800849722744
    Need changing reaction: 
    3OAR100
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100']
    Calibration round 13 : 0.3264473195104685
    Need changing reaction: 
    3OAR40
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40']
    Calibration round 14 : 0.3331813808664081
    Need changing reaction: 
    PGCD
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD']
    Calibration round 15 : 0.34013803712390456
    Need changing reaction: 
    PFK_num2
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2']
    Calibration round 16 : 0.3427260533654994
    Need changing reaction: 
    PGM_reverse_num2
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2']
    Calibration round 17 : 0.3427260533654994
    Need changing reaction: 
    PSP_L
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L']
    Calibration round 18 : 0.34773323581020205
    Need changing reaction: 
    ACCOAC
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC']
    Calibration round 19 : 0.3510571116510257
    Need changing reaction: 
    PTAr_num2
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2']
    Calibration round 20 : 0.3578952500160452
    Need changing reaction: 
    3OAR120
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120']
    Calibration round 21 : 0.36237199885300264
    Need changing reaction: 
    3OAR140
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140']
    Calibration round 22 : 0.3669621615609515
    Need changing reaction: 
    HSDy_reverse_num1
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140']
    Calibration round 23 : 0.3669621615609515
    Need changing reaction: 
    3OAR161
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161']
    Calibration round 24 : 0.370397422352087
    Need changing reaction: 
    3OAR121
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121']
    Calibration round 25 : 0.3738976082790336
    Need changing reaction: 
    3OAR141
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141']
    Calibration round 26 : 0.379441609584257
    Need changing reaction: 
    ADSS
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS']
    Calibration round 27 : 0.3827948676142595
    Need changing reaction: 
    3OAR160
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160']
    Calibration round 28 : 0.38391758177386637
    Need changing reaction: 
    ACONTa_num1
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1']
    Calibration round 29 : 0.3866375337831847
    Need changing reaction: 
    IPPS
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS']
    Calibration round 30 : 0.3892744540703348
    Need changing reaction: 
    SADT2
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2']
    Calibration round 31 : 0.39191029957143503
    Need changing reaction: 
    GAPD
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2']
    Calibration round 32 : 0.39191029957143503
    Need changing reaction: 
    GLNS
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS']
    Calibration round 33 : 0.39392994176634
    Need changing reaction: 
    ASNS2_num1
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1']
    Calibration round 34 : 0.3953059242452846
    Need changing reaction: 
    IPPMIa_reverse
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1']
    Calibration round 35 : 0.3953059242452846
    Need changing reaction: 
    IPPMIb_reverse
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1']
    Calibration round 36 : 0.3953059242452846
    Need changing reaction: 
    ACOTA_reverse_num2
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1']
    Calibration round 37 : 0.3953059242452846
    Need changing reaction: 
    AICART
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART']
    Calibration round 38 : 0.39574674998410425
    Need changing reaction: 
    TPI
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART', 'TPI']
    Calibration round 39 : 0.39619632518648024
    Need changing reaction: 
    THD2pp
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART', 'TPI', 'THD2pp']
    Calibration round 40 : 0.3965822171406823
    Need changing reaction: 
    IMPC_reverse
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART', 'TPI', 'THD2pp']
    Calibration round 41 : 0.3965822171406823
    Need changing reaction: 
    PAPPT3
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART', 'TPI', 'THD2pp', 'PAPPT3']
    Calibration round 42 : 0.3975058418569194
    Need changing reaction: 
    AIRC2
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART', 'TPI', 'THD2pp', 'PAPPT3', 'AIRC2']
    Calibration round 43 : 0.39820705457384425
    Need changing reaction: 
    ASPK_num3
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART', 'TPI', 'THD2pp', 'PAPPT3', 'AIRC2', 'ASPK_num3']
    Calibration round 44 : 0.39890635782609746
    Need changing reaction: 
    THRS
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART', 'TPI', 'THD2pp', 'PAPPT3', 'AIRC2', 'ASPK_num3', 'THRS']
    Calibration round 45 : 0.399605175392766
    Need changing reaction: 
    AIRC3_reverse
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART', 'TPI', 'THD2pp', 'PAPPT3', 'AIRC2', 'ASPK_num3', 'THRS']
    Calibration round 46 : 0.399605175392766
    Need changing reaction: 
    CS
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART', 'TPI', 'THD2pp', 'PAPPT3', 'AIRC2', 'ASPK_num3', 'THRS', 'CS']
    Calibration round 47 : 0.4001983560070541
    Need changing reaction: 
    IMPD
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART', 'TPI', 'THD2pp', 'PAPPT3', 'AIRC2', 'ASPK_num3', 'THRS', 'CS', 'IMPD']
    Calibration round 48 : 0.4007962534772659
    Need changing reaction: 
    ANS
    Changed reaction: 
    ['PRFGS', 'FBA_num1', 'PFL_num3', '3OAS121', '3OAS161', '3OAS141', 'PSERT', 'PAPSR2_num2', '3OAR60', '3OAR80', '3OAR100', '3OAR40', 'PGCD', 'PFK_num2', 'PSP_L', 'ACCOAC', 'PTAr_num2', '3OAR120', '3OAR140', '3OAR161', '3OAR121', '3OAR141', 'ADSS', '3OAR160', 'ACONTa_num1', 'IPPS', 'SADT2', 'GLNS', 'ASNS2_num1', 'AICART', 'TPI', 'THD2pp', 'PAPPT3', 'AIRC2', 'ASPK_num3', 'THRS', 'CS', 'IMPD', 'ANS']
    Calibration round 49 : 0.40133927002029507
    


```python

```
