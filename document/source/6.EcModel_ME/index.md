# 6.EcModel_ME


```python
import cobra
#from script.ECMpy_function import *
import sys
sys.path.append(r'./script/')
from ECMpy_function import *
```



```python
ecModel_file="./model/eciML1515.json"
obj='EX_trp__L_e' # EX_lys_L_e
substrate = 'EX_glc__D_e'
substrate_con=10
biomass_id = 'BIOMASS_Ec_iML1515_core_75p37M'
biomass_min=0.1

fluxes_outfile = './analysis/ECMpy_solution_%s_%.2g_pfba.csv'%(obj,substrate_con)
enzcost_cost_file = './analysis/enzcost_%s_%.2g.csv'%(obj,substrate_con)
enzcost_diff_file = './analysis/enzcost_diff_%s_%.2g.csv'%(obj,substrate_con)
enzcost_diff_selcet_file = './analysis/enzcost_diff_%s_%.2g_select.csv'%(obj,substrate_con)
FSEOF_file = './analysis/FESEOF_%s_%.2g.csv'%(obj,substrate_con)
```

## Directly determining targets based on enzyme abundance


```python
enzcost = pd.DataFrame()

enz_model=get_enzyme_constraint_model(ecModel_file)
enz_model.reactions.get_by_id(substrate).bounds=(-substrate_con,0)#EX_glc_e EX_fru_e EX_sucr_e EX_inost_e EX_lac_D_e EX_ac_e
enz_model.reactions.get_by_id('%s_reverse'%substrate).bounds=(0,0)
enz_model.reactions.get_by_id(biomass_id).bounds=(biomass_min,biomass_min)
enz_model.objective=obj
enz_model_pfba_solution = cobra.flux_analysis.pfba(enz_model)

enzcost = get_fluxes_detail_in_model(enz_model,enz_model_pfba_solution,fluxes_outfile,ecModel_file)
enzcost['enz_ratio'] = enzcost['E']/np.sum(enzcost['E'])
enzcost = enzcost.sort_values('enz_ratio', ascending=False)
enzcost_select = enzcost[enzcost['enz_ratio']>0.01]
enzcost_select.to_csv(enzcost_cost_file)
enzcost_select.head(10)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fluxes</th>
      <th>kcat_MW</th>
      <th>E</th>
      <th>equ</th>
      <th>ec-code</th>
      <th>enz_ratio</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ANS</th>
      <td>3.594316</td>
      <td>76.627794</td>
      <td>0.046906</td>
      <td>chor_c + gln__L_c --&gt; anth_c + glu__L_c + h_c ...</td>
      <td>4.1.3.27</td>
      <td>0.206635</td>
    </tr>
    <tr>
      <th>TRPAS2_reverse</th>
      <td>3.594316</td>
      <td>85.270150</td>
      <td>0.042152</td>
      <td>indole_c + nh4_c + pyr_c --&gt; h2o_c + trp__L_c</td>
      <td>NaN</td>
      <td>0.185692</td>
    </tr>
    <tr>
      <th>ANPRT</th>
      <td>3.594316</td>
      <td>158.255810</td>
      <td>0.022712</td>
      <td>anth_c + prpp_c --&gt; ppi_c + pran_c</td>
      <td>2.4.2.18</td>
      <td>0.100053</td>
    </tr>
    <tr>
      <th>TRPS3</th>
      <td>3.594316</td>
      <td>249.207417</td>
      <td>0.014423</td>
      <td>3ig3p_c --&gt; g3p_c + indole_c</td>
      <td>4.1.2.8,4.2.1.20</td>
      <td>0.063537</td>
    </tr>
    <tr>
      <th>PRAIi</th>
      <td>3.594316</td>
      <td>363.692938</td>
      <td>0.009883</td>
      <td>pran_c --&gt; 2cpr5p_c</td>
      <td>5.3.1.24</td>
      <td>0.043537</td>
    </tr>
    <tr>
      <th>IGPS</th>
      <td>3.594316</td>
      <td>363.692938</td>
      <td>0.009883</td>
      <td>2cpr5p_c + h_c --&gt; 3ig3p_c + co2_c + h2o_c</td>
      <td>4.1.1.48</td>
      <td>0.043537</td>
    </tr>
    <tr>
      <th>GLCDpp_num2</th>
      <td>5.074624</td>
      <td>1061.316075</td>
      <td>0.004781</td>
      <td>glc__D_p + h2o_p + q8_c --&gt; glcn_p + h_p + q8h2_c</td>
      <td>NaN</td>
      <td>0.021064</td>
    </tr>
    <tr>
      <th>PRPPS</th>
      <td>3.706289</td>
      <td>841.654904</td>
      <td>0.004404</td>
      <td>atp_c + r5p_c --&gt; amp_c + h_c + prpp_c</td>
      <td>2.7.6.1</td>
      <td>0.019399</td>
    </tr>
    <tr>
      <th>GNK_num2</th>
      <td>5.074624</td>
      <td>1273.803729</td>
      <td>0.003984</td>
      <td>atp_c + glcn_c --&gt; 6pgc_c + adp_c + h_c</td>
      <td>2.7.1.12</td>
      <td>0.017550</td>
    </tr>
    <tr>
      <th>PDH</th>
      <td>2.960121</td>
      <td>784.108302</td>
      <td>0.003775</td>
      <td>coa_c + nad_c + pyr_c --&gt; accoa_c + co2_c + na...</td>
      <td>1.2.1.-,1.2.1.51,1.2.4.1,1.8.1.4,2.3.1.12</td>
      <td>0.016631</td>
    </tr>
  </tbody>
</table>
</div>



##  Determine targets based on the fold changes of enzyme cost 


```python
biomass_max=0.6
FC_threshold = 1.5

enzcost_diff = get_enz_foldchange(ecModel_file,obj,substrate,substrate_con,biomass_id,biomass_min,biomass_max,FC_threshold,fluxes_outfile, enzcost_diff_file)
enzcost_diff_select = enzcost_diff[enzcost_diff['E_μ0.1']>0.0001]
enzcost_diff_select.to_csv(enzcost_diff_selcet_file)

enzcost_diff_select.head(10)
```

    /home/maozt/anaconda3/envs/ECMpy2/lib/python3.7/site-packages/pandas/core/arraylike.py:364: RuntimeWarning: invalid value encountered in log2
      result = getattr(ufunc, method)(*inputs, **kwargs)
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>E_μ0.1</th>
      <th>E_μ0.6</th>
      <th>log2_foldchange(max/min)</th>
      <th>log2_foldchange(min/max)</th>
      <th>type</th>
      <th>foldchange(max/min)</th>
      <th>equ</th>
      <th>ec-code</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>3OAS160_num1</th>
      <td>0.000117</td>
      <td>0.000702</td>
      <td>2.584963</td>
      <td>-2.584963</td>
      <td>weaken_target</td>
      <td>6.0000</td>
      <td>h_c + malACP_c + myrsACP_c --&gt; 3opalmACP_c + A...</td>
      <td>2.3.1.86</td>
    </tr>
    <tr>
      <th>3OAR160</th>
      <td>0.000140</td>
      <td>0.000842</td>
      <td>2.584963</td>
      <td>-2.584963</td>
      <td>weaken_target</td>
      <td>6.0000</td>
      <td>3opalmACP_c + h_c + nadph_c --&gt; 3hpalmACP_c + ...</td>
      <td>1.1.1.100,2.3.1.-,2.3.1.85,2.3.1.86</td>
    </tr>
    <tr>
      <th>3OAR140</th>
      <td>0.000449</td>
      <td>0.002693</td>
      <td>2.584963</td>
      <td>-2.584963</td>
      <td>weaken_target</td>
      <td>6.0000</td>
      <td>3omrsACP_c + h_c + nadph_c --&gt; 3hmrsACP_c + na...</td>
      <td>1.1.1.100,2.3.1.-,2.3.1.85,2.3.1.86</td>
    </tr>
    <tr>
      <th>3OAS140_num1</th>
      <td>0.000115</td>
      <td>0.000690</td>
      <td>2.584963</td>
      <td>-2.584963</td>
      <td>weaken_target</td>
      <td>6.0000</td>
      <td>ddcaACP_c + h_c + malACP_c --&gt; 3omrsACP_c + AC...</td>
      <td>2.3.1.86</td>
    </tr>
    <tr>
      <th>ACONTa_num1</th>
      <td>0.000580</td>
      <td>0.002436</td>
      <td>2.070767</td>
      <td>-2.070767</td>
      <td>weaken_target</td>
      <td>4.2011</td>
      <td>cit_c --&gt; acon_C_c + h2o_c</td>
      <td>4.2.1.3</td>
    </tr>
    <tr>
      <th>EAR121x</th>
      <td>0.000119</td>
      <td>0.000714</td>
      <td>2.584963</td>
      <td>-2.584963</td>
      <td>weaken_target</td>
      <td>6.0000</td>
      <td>h_c + nadh_c + t3c5ddeceACP_c --&gt; cddec5eACP_c...</td>
      <td>1.3.1.9</td>
    </tr>
    <tr>
      <th>3OAR161</th>
      <td>0.000165</td>
      <td>0.000992</td>
      <td>2.584963</td>
      <td>-2.584963</td>
      <td>weaken_target</td>
      <td>6.0000</td>
      <td>3ocpalm9eACP_c + h_c + nadph_c --&gt; 3hcpalm9eAC...</td>
      <td>1.1.1.100</td>
    </tr>
    <tr>
      <th>3OAR141</th>
      <td>0.000165</td>
      <td>0.000992</td>
      <td>2.584963</td>
      <td>-2.584963</td>
      <td>weaken_target</td>
      <td>6.0000</td>
      <td>3ocmrs7eACP_c + h_c + nadph_c --&gt; 3hcmrs7eACP_...</td>
      <td>1.1.1.100</td>
    </tr>
    <tr>
      <th>EAR160y</th>
      <td>0.000153</td>
      <td>0.000918</td>
      <td>2.584963</td>
      <td>-2.584963</td>
      <td>weaken_target</td>
      <td>6.0000</td>
      <td>h_c + nadph_c + tpalm2eACP_c --&gt; nadp_c + palm...</td>
      <td>1.3.1.10</td>
    </tr>
    <tr>
      <th>EAR161y</th>
      <td>0.000179</td>
      <td>0.001072</td>
      <td>2.584963</td>
      <td>-2.584963</td>
      <td>weaken_target</td>
      <td>6.0000</td>
      <td>h_c + nadph_c + t3c9palmeACP_c --&gt; hdeACP_c + ...</td>
      <td>1.3.1.10</td>
    </tr>
  </tbody>
</table>
</div>



## Determine targets based on FESOF


```python
model=get_enzyme_constraint_model(ecModel_file)

FSEOFdf_done = run_FSEOF(model,substrate,substrate_con,biomass_id,obj,FSEOF_file)

FSEOFdf_done.head(10)
```

    ./script/ECMpy_function.py:2653: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      always_down['regulation'] = 'down'
    ./script/ECMpy_function.py:2656: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      always_up['regulation'] = 'up'
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>growth = 0.34</th>
      <th>FC growth = 0.34</th>
      <th>growth = 0.37</th>
      <th>FC growth = 0.37</th>
      <th>growth = 0.4</th>
      <th>FC growth = 0.4</th>
      <th>growth = 0.43</th>
      <th>FC growth = 0.43</th>
      <th>growth = 0.46</th>
      <th>FC growth = 0.46</th>
      <th>...</th>
      <th>FC growth = 0.55</th>
      <th>growth = 0.58</th>
      <th>FC growth = 0.58</th>
      <th>growth = 0.61</th>
      <th>FC growth = 0.61</th>
      <th>WT</th>
      <th>GPR</th>
      <th>FC_mean</th>
      <th>regulation</th>
      <th>reactions</th>
    </tr>
    <tr>
      <th>id</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>H2Otpp_reverse_num2</th>
      <td>39.134690</td>
      <td>1000.000000</td>
      <td>38.784008</td>
      <td>1000.000000</td>
      <td>38.433325</td>
      <td>1000.000000</td>
      <td>38.082642</td>
      <td>1000.000000</td>
      <td>37.731959</td>
      <td>1000.000000</td>
      <td>...</td>
      <td>1000.000000</td>
      <td>36.083873</td>
      <td>1000.000000</td>
      <td>35.649234</td>
      <td>1000.000000</td>
      <td>0.000000</td>
      <td>b0875</td>
      <td>562.598358</td>
      <td>up</td>
      <td>h2o_c --&gt; h2o_p</td>
    </tr>
    <tr>
      <th>H2Otex_reverse_num8</th>
      <td>39.134690</td>
      <td>1000.000000</td>
      <td>38.784008</td>
      <td>1000.000000</td>
      <td>38.433325</td>
      <td>1000.000000</td>
      <td>38.082642</td>
      <td>1000.000000</td>
      <td>37.731959</td>
      <td>1000.000000</td>
      <td>...</td>
      <td>1000.000000</td>
      <td>36.083873</td>
      <td>1000.000000</td>
      <td>35.649234</td>
      <td>1000.000000</td>
      <td>0.000000</td>
      <td>b0929</td>
      <td>562.598358</td>
      <td>up</td>
      <td>h2o_p --&gt; h2o_e</td>
    </tr>
    <tr>
      <th>CO2tex_reverse_num4</th>
      <td>13.119677</td>
      <td>1000.000000</td>
      <td>13.437758</td>
      <td>1000.000000</td>
      <td>13.755840</td>
      <td>1000.000000</td>
      <td>14.073921</td>
      <td>1000.000000</td>
      <td>14.392003</td>
      <td>1000.000000</td>
      <td>...</td>
      <td>1000.000000</td>
      <td>15.630991</td>
      <td>1000.000000</td>
      <td>15.937058</td>
      <td>1000.000000</td>
      <td>0.000000</td>
      <td>b2215</td>
      <td>551.995516</td>
      <td>up</td>
      <td>co2_p --&gt; co2_e</td>
    </tr>
    <tr>
      <th>Htex_reverse_num4</th>
      <td>12.132625</td>
      <td>1000.000000</td>
      <td>12.272701</td>
      <td>1000.000000</td>
      <td>12.412777</td>
      <td>1000.000000</td>
      <td>12.552853</td>
      <td>1000.000000</td>
      <td>12.692929</td>
      <td>1000.000000</td>
      <td>...</td>
      <td>1000.000000</td>
      <td>13.362659</td>
      <td>1000.000000</td>
      <td>13.540216</td>
      <td>1000.000000</td>
      <td>0.000000</td>
      <td>b2215</td>
      <td>551.227208</td>
      <td>up</td>
      <td>h_p --&gt; h_e</td>
    </tr>
    <tr>
      <th>ACtex_reverse_num4</th>
      <td>4.749431</td>
      <td>1000.000000</td>
      <td>4.979100</td>
      <td>1000.000000</td>
      <td>5.208770</td>
      <td>1000.000000</td>
      <td>5.438439</td>
      <td>1000.000000</td>
      <td>5.668108</td>
      <td>1000.000000</td>
      <td>...</td>
      <td>1000.000000</td>
      <td>6.749216</td>
      <td>1000.000000</td>
      <td>7.034352</td>
      <td>1000.000000</td>
      <td>0.000000</td>
      <td>b2215</td>
      <td>548.035665</td>
      <td>up</td>
      <td>ac_p --&gt; ac_e</td>
    </tr>
    <tr>
      <th>TRPtex_reverse_num4</th>
      <td>2.129569</td>
      <td>1000.000000</td>
      <td>1.946947</td>
      <td>1000.000000</td>
      <td>1.764324</td>
      <td>1000.000000</td>
      <td>1.581701</td>
      <td>1000.000000</td>
      <td>1.399079</td>
      <td>1000.000000</td>
      <td>...</td>
      <td>1000.000000</td>
      <td>0.642086</td>
      <td>1000.000000</td>
      <td>0.450470</td>
      <td>1000.000000</td>
      <td>0.000000</td>
      <td>b2215</td>
      <td>546.089712</td>
      <td>up</td>
      <td>trp__L_p --&gt; trp__L_e</td>
    </tr>
    <tr>
      <th>TRPt2rpp_reverse_num1</th>
      <td>2.129569</td>
      <td>1000.000000</td>
      <td>1.946947</td>
      <td>1000.000000</td>
      <td>1.764324</td>
      <td>1000.000000</td>
      <td>1.581701</td>
      <td>1000.000000</td>
      <td>1.399079</td>
      <td>1000.000000</td>
      <td>...</td>
      <td>1000.000000</td>
      <td>0.642086</td>
      <td>1000.000000</td>
      <td>0.450470</td>
      <td>1000.000000</td>
      <td>0.000000</td>
      <td>b3709</td>
      <td>546.089712</td>
      <td>up</td>
      <td>h_c + trp__L_c --&gt; h_p + trp__L_p</td>
    </tr>
    <tr>
      <th>ANS</th>
      <td>2.148896</td>
      <td>55.574245</td>
      <td>1.967978</td>
      <td>50.895402</td>
      <td>1.787061</td>
      <td>46.216558</td>
      <td>1.606144</td>
      <td>41.537715</td>
      <td>1.425226</td>
      <td>36.858872</td>
      <td>...</td>
      <td>22.369518</td>
      <td>0.675055</td>
      <td>17.458107</td>
      <td>0.485145</td>
      <td>12.546695</td>
      <td>0.038667</td>
      <td>b1263 and b1264</td>
      <td>22.007804</td>
      <td>up</td>
      <td>chor_c + gln__L_c --&gt; anth_c + glu__L_c + h_c ...</td>
    </tr>
    <tr>
      <th>TRPS3</th>
      <td>2.148896</td>
      <td>55.574245</td>
      <td>1.967978</td>
      <td>50.895402</td>
      <td>1.787061</td>
      <td>46.216558</td>
      <td>1.606144</td>
      <td>41.537715</td>
      <td>1.425226</td>
      <td>36.858872</td>
      <td>...</td>
      <td>22.369518</td>
      <td>0.675055</td>
      <td>17.458107</td>
      <td>0.485145</td>
      <td>12.546695</td>
      <td>0.038667</td>
      <td>b1260 and b1261</td>
      <td>22.007804</td>
      <td>up</td>
      <td>3ig3p_c --&gt; g3p_c + indole_c</td>
    </tr>
    <tr>
      <th>ANPRT</th>
      <td>2.148896</td>
      <td>55.574245</td>
      <td>1.967978</td>
      <td>50.895402</td>
      <td>1.787061</td>
      <td>46.216558</td>
      <td>1.606144</td>
      <td>41.537715</td>
      <td>1.425226</td>
      <td>36.858872</td>
      <td>...</td>
      <td>22.369518</td>
      <td>0.675055</td>
      <td>17.458107</td>
      <td>0.485145</td>
      <td>12.546695</td>
      <td>0.038667</td>
      <td>b1263</td>
      <td>22.007804</td>
      <td>up</td>
      <td>anth_c + prpp_c --&gt; ppi_c + pran_c</td>
    </tr>
  </tbody>
</table>
<p>10 rows × 25 columns</p>
</div>


