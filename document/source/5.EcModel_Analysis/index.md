# 5.EcModel Analysis



```python
import cobra
import re 
import pandas as pd
import numpy as np
#from script.ECMpy_function import *
import sys
sys.path.append(r'./script/')
from ECMpy_function import *
```

##  Cumulative distribution figure




```python
ecModel_file="./model/iML1515_irr_enz_constraint_adj.json"

method='AutoPACMEN'#DLKcat
reaction_kcat_MW_file = "./analysis/get_kcat_mw_by_%s/reaction_change_by_enzuse.csv"%method
reaction_kcat_MW = pd.read_csv(reaction_kcat_MW_file)
reaction_kcat_MW = round(reaction_kcat_MW,3)
reaction_kcat_dis_file='./analysis/reaction_kcat_distrbution.png'
reaction_mw_dis_file='./analysis/reaction_mw_distrbution.png'
```


```python
reaction_kcat_select = reaction_kcat_MW[reaction_kcat_MW['data_type'] != 'fill']
#Sort values
sorted_data = reaction_kcat_select.sort_values('kcat')
sorted_data = sorted_data.reset_index(drop=True)
y_index = sorted_data.index / (sorted_data.shape[0]  - 1)
data_cdf_data = sorted_data['kcat']
x_name="<b>kcat(1/s)<b>"
y_name="<b>Cummulative distribution<b>"
nticks=1000
fig=draw_cdf_fig(data_cdf_data,reaction_kcat_dis_file,x_name,y_name,y_index,nticks)
fig.show()
```


```python
reaction_kcat_select = reaction_kcat_MW[reaction_kcat_MW['data_type'] != 'fill']
#Sort values
sorted_data = reaction_kcat_select.sort_values('MW')
sorted_data = sorted_data.reset_index(drop=True)
y_index = sorted_data.index / (sorted_data.shape[0]  - 1)
data_cdf_data = sorted_data['MW']
y_index = sorted_data.index / (sorted_data.shape[0]  - 1)
data_cdf_data = data_cdf_data/1000# kDa
x_name="<b>mass(kDa)<b>"
y_name="<b>Cummulative distribution<b>"
nticks=10000
fig=draw_cdf_fig(data_cdf_data,reaction_mw_dis_file,x_name,y_name,y_index,nticks)
fig.show()
```

## Phenotype Phase Plane (PhPP) Analysis

### Input and output files


```python
ecModel_file="./model/iML1515_irr_enz_constraint_adj.json"
obj='BIOMASS_Ec_iML1515_core_75p37M'
z_id='EX_o2_e'
x_id='EX_glc__D_e'
substrate_bound=10
o2_bound=40

GEM_output_file='./analysis/iML1515_glc_o2_df.csv'
ecGEM_output_file='./analysis/eciML1515_glc_o2_df.csv'

PhPP_output_fig_file='./analysis/PhPP_combine.png'
```


```python
GEM_glc_o2_df=get_PhPP_data(ecModel_file, 'GEM', obj, substrate_bound,o2_bound,11, GEM_output_file,x_id,z_id)
GEM_glc_o2_df.drop(0,axis=0,inplace=True)
GEM_glc_o2_df.drop(0,axis=1,inplace=True)

ecGEM_glc_o2_df=get_PhPP_data(ecModel_file, 'ecGEM', obj, substrate_bound,o2_bound,11, ecGEM_output_file,x_id,z_id)
ecGEM_glc_o2_df.drop(0,axis=0,inplace=True)
ecGEM_glc_o2_df.drop(0,axis=1,inplace=True)
```


```python
fig=draw_3d_rbas(GEM_glc_o2_df,substrate_bound,o2_bound,1,11,PhPP_output_fig_file)
fig.show()
```


```python
fig=draw_3d_rbas(ecGEM_glc_o2_df,substrate_bound,o2_bound,1,11,PhPP_output_fig_file)
fig.show()
```


```python
fig=drawphpp(GEM_glc_o2_df,ecGEM_glc_o2_df,substrate_bound,o2_bound,1,11,PhPP_output_fig_file)
fig.show()
```

## Overflow simulation

### Input and output files


```python
# inputfiles
json_model_file="./model/iML1515_irr_enz_constraint_adj.json"
enz_model=get_enzyme_constraint_model(json_model_file)
norm_model=cobra.io.json.load_json_model(json_model_file)
method='AutoPACMEN'#DLKcat
reaction_kcat_MW_file = "./analysis/get_kcat_mw_by_%s/reaction_change_by_enzuse.csv"%method
reaction_kcat_MW = pd.read_csv(reaction_kcat_MW_file,index_col=0)
reaction_kcat_MW = round(reaction_kcat_MW,3)

# outputfiles
overflow_result_figfile="./analysis/pfba_overflow_result.png"
```


```python
use_substrate='EX_glc__D_e'
substrate_name='glucose'
glc_concentration_list = np.arange(1, 10, 0.5)
columns = ['biomass', 'acetate', 'O2', 'CO2', 'pyruvate', 'ethanol']
correspond_rxn = ['BIOMASS_Ec_iML1515_core_75p37M', 'EX_ac_e', 'EX_o2_e_reverse', 'EX_co2_e', 'EX_pyr_e', 'EX_etoh_e']
GEMyield_list = pd.DataFrame()
ecGEMyield_list = pd.DataFrame()

# Calculate yield for growth_model
with enz_model as growth_model:
    for glc_concentration in glc_concentration_list:
        ecGEMyield_list=calculate_yield(growth_model, use_substrate, substrate_name, glc_concentration,columns, correspond_rxn, ecGEMyield_list)

# Calculate yield for norm_model
with norm_model as growth_model:
    for glc_concentration in glc_concentration_list:
        GEMyield_list=calculate_yield(growth_model, use_substrate, substrate_name, glc_concentration,columns, correspond_rxn, GEMyield_list)

```


```python
substrate_name='glucose'
substrate_bound=10
obj_bound=1
secrate_bound=18
column_list=['biomass','CO2','acetate']
y_axis_loc_list=['left','right','right','right']
color_list = generate_random_colors(len(column_list)*2)

fig=draw_overfolw_fig(GEMyield_list,ecGEMyield_list,column_list,y_axis_loc_list,color_list,substrate_name, substrate_bound,obj_bound,secrate_bound,overflow_result_figfile)
fig.show()
```


<div>                            <div id="4a4a3ba5-90d2-4f51-9a31-5294ee828740" class="plotly-graph-div" style="height:600px; width:800px;"></div>            <script type="text/javascript">                require(["plotly"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("4a4a3ba5-90d2-4f51-9a31-5294ee828740")) {                    Plotly.newPlot(                        "4a4a3ba5-90d2-4f51-9a31-5294ee828740",                        [{"line":{"color":"#ca4553","width":3},"marker":{"color":"#ca4553","size":10},"mode":"lines+markers","name":"biomass_GEM","x":[1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5],"y":[0.06384218390161035,0.1087122535335538,0.15348051168052224,0.19824876982749023,0.24301702797445823,0.287785286121426,0.3325535442683941,0.3773218024153618,0.4220900605623304,0.46685831870929845,0.5116265768562658,0.5563948350032342,0.6011630931502019,0.6459313512971697,0.6906996094441386,0.7354678675911055,0.7802361257380747,0.8250043838850418],"type":"scatter"},{"line":{"color":"#8dff1c","width":3},"marker":{"color":"#8dff1c","size":10},"mode":"lines+markers","name":"biomass_ecGEM","x":[1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5],"y":[0.06384218390161003,0.10871225353355472,0.1532503285905529,0.19249399988532345,0.22703761518513996,0.25340754953291406,0.27977748388068563,0.3032547093347451,0.3153727891418566,0.3274908689489677,0.3396089487560784,0.3517270285631926,0.36282028387129345,0.37131079911682985,0.37794932792798575,0.38379698033014087,0.38964463273229244,0.39549228513444823],"type":"scatter"},{"line":{"color":"#9c305c","width":3},"marker":{"color":"#9c305c","size":10},"mode":"lines+markers","name":"CO2_GEM","x":[1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5],"xaxis":"x2","y":[3.379571641831732,4.537861416540408,5.700330084985212,6.86279875343011,8.025267421874993,9.18773609031987,10.350204758764734,11.51267342720962,12.6751420956545,13.83761076409942,15.00007943254431,16.162548100989156,17.325016769434008,18.48748543787895,19.649954106323815,20.812422774768727,21.97489144321361,23.137360111658502],"yaxis":"y2","type":"scatter"},{"line":{"color":"#10beaf","width":3},"marker":{"color":"#10beaf","size":10},"mode":"lines+markers","name":"CO2_ecGEM","x":[1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5],"xaxis":"x2","y":[3.3795716418317574,4.537861416540304,5.709778043367395,7.099005622142966,8.681148907454311,10.598784245956436,12.516419584458838,14.050723922670315,13.654794821234857,13.258865719787599,12.862936618361116,12.467007516932425,12.052914438425068,11.587044681038462,10.683892479013704,11.427936503291878,12.256430418100464,13.0849243431091],"yaxis":"y2","type":"scatter"},{"line":{"color":"#8e5b45","width":3},"marker":{"color":"#8e5b45","size":10},"mode":"lines+markers","name":"acetate_GEM","x":[1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5],"xaxis":"x2","y":[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],"yaxis":"y2","type":"scatter"},{"line":{"color":"#133647","width":3},"marker":{"color":"#133647","size":10},"mode":"lines+markers","name":"acetate_ecGEM","x":[1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5],"xaxis":"x2","y":[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.25103170563648847,1.7003005035425784,3.149569301454742,4.598838099356203,6.048106897258854,7.527489852521172,9.08617640851334,10.237881679387232,11.148729873710707,12.059578068034142,12.970426262357474],"yaxis":"y2","type":"scatter"}],                        {"height":600,"legend":{"font":{"color":"black","family":"Times New Roman","size":20},"x":0.05,"y":0.98},"plot_bgcolor":"white","width":800,"xaxis":{"linecolor":"black","range":[1,10],"tickcolor":"black","tickfont":{"color":"black","family":"Times New Roman","size":15},"ticks":"inside","title":{"font":{"family":"Times New Roman","size":20},"text":"\u003cb\u003eSubstrate uptake rate (mmol\u002fgDW\u002fh)\u003c\u002fb\u003e"}},"yaxis":{"linecolor":"black","range":[0,1],"tickcolor":"black","tickfont":{"color":"black","family":"Times New Roman","size":15},"ticks":"inside","title":{"font":{"family":"Times New Roman","size":20},"text":"\u003cb\u003eGrowth rate (h\u003csup\u003e-1\u003c\u002fsup\u003e)\u003c\u002fb\u003e"}},"xaxis2":{"linecolor":"black","overlaying":"x","range":[1,10],"showticklabels":false,"side":"top"},"yaxis2":{"linecolor":"black","overlaying":"y","range":[0,18],"side":"right","tickcolor":"black","tickfont":{"color":"black","family":"Times New Roman","size":15},"ticks":"inside","title":{"font":{"family":"Times New Roman","size":20},"text":"\u003cb\u003eSecrete rate (mmol\u002fgDW\u002fh)\u003c\u002fb\u003e"}},"template":{"data":{"histogram2dcontour":[{"type":"histogram2dcontour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"choropleth":[{"type":"choropleth","colorbar":{"outlinewidth":0,"ticks":""}}],"histogram2d":[{"type":"histogram2d","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmap":[{"type":"heatmap","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmapgl":[{"type":"heatmapgl","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"contourcarpet":[{"type":"contourcarpet","colorbar":{"outlinewidth":0,"ticks":""}}],"contour":[{"type":"contour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"surface":[{"type":"surface","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"mesh3d":[{"type":"mesh3d","colorbar":{"outlinewidth":0,"ticks":""}}],"scatter":[{"fillpattern":{"fillmode":"overlay","size":10,"solidity":0.2},"type":"scatter"}],"parcoords":[{"type":"parcoords","line":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolargl":[{"type":"scatterpolargl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"scattergeo":[{"type":"scattergeo","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolar":[{"type":"scatterpolar","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"scattergl":[{"type":"scattergl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatter3d":[{"type":"scatter3d","line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattermapbox":[{"type":"scattermapbox","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterternary":[{"type":"scatterternary","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattercarpet":[{"type":"scattercarpet","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"type":"carpet"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}],"barpolar":[{"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"pie":[{"automargin":true,"type":"pie"}]},"layout":{"autotypenumbers":"strict","colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"hovermode":"closest","hoverlabel":{"align":"left"},"paper_bgcolor":"white","plot_bgcolor":"#E5ECF6","polar":{"bgcolor":"#E5ECF6","angularaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"radialaxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"ternary":{"bgcolor":"#E5ECF6","aaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"baxis":{"gridcolor":"white","linecolor":"white","ticks":""},"caxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]]},"xaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"yaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"scene":{"xaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"yaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"zaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"geo":{"bgcolor":"white","landcolor":"#E5ECF6","subunitcolor":"white","showland":true,"showlakes":true,"lakecolor":"white"},"title":{"x":0.05},"mapbox":{"style":"light"}}}},                        {"responsive": true}                    ).then(function(){

var gd = document.getElementById('4a4a3ba5-90d2-4f51-9a31-5294ee828740');
var x = new MutationObserver(function (mutations, observer) {{
        var display = window.getComputedStyle(gd).display;
        if (!display || display === 'none') {{
            console.log([gd, 'removed!']);
            Plotly.purge(gd);
            observer.disconnect();
        }}
}});

// Listen for the removal of the full notebook cells
var notebookContainer = gd.closest('#notebook-container');
if (notebookContainer) {{
    x.observe(notebookContainer, {childList: true});
}}

// Listen for the clearing of the current output cell
var outputEl = gd.closest('.output');
if (outputEl) {{
    x.observe(outputEl, {childList: true});
}}

                        })                };                });            </script>        </div>


## Trade-off simulation


```python
# inputfiles
method='AutoPACMEN'#DLKcat
reaction_kcat_MW_file = "./analysis/get_kcat_mw_by_%s/reaction_change_by_enzuse.csv"%method
reaction_kcat_MW = pd.read_csv(reaction_kcat_MW_file,index_col=0)
reaction_kcat_MW = round(reaction_kcat_MW,3)
json_model_file="./model/iML1515_irr_enz_constraint_adj.json"
enz_model=get_enzyme_constraint_model(json_model_file)
glc_concentration_list = np.arange(1, 10, 0.5)
efficiency_file="./analysis/efficiency_pfba.csv"
trade_off_enzyme_efficiency_figfile="./analysis/trade_off_enzyme_efficiency.png"
use_substrate='EX_glc__D_e'
obj='BIOMASS_Ec_iML1515_core_75p37M'

yield_cost_efficiency_df=get_yield_cost_efficiency(enz_model,glc_concentration_list,use_substrate,obj,reaction_kcat_MW,efficiency_file)
yield_cost_efficiency_df.head()

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
      <th>biomass</th>
      <th>glucose_simu</th>
      <th>glucose_set</th>
      <th>biomass yield</th>
      <th>min enzyme cost</th>
      <th>enzyme efficiency</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1.0</th>
      <td>0.089537</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.497425</td>
      <td>0.099772</td>
      <td>0.897408</td>
    </tr>
    <tr>
      <th>1.5</th>
      <td>0.134305</td>
      <td>1.5</td>
      <td>1.5</td>
      <td>0.497425</td>
      <td>0.149658</td>
      <td>0.897408</td>
    </tr>
    <tr>
      <th>2.0</th>
      <td>0.179073</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>0.497425</td>
      <td>0.199545</td>
      <td>0.897408</td>
    </tr>
    <tr>
      <th>2.5</th>
      <td>0.221429</td>
      <td>2.5</td>
      <td>2.5</td>
      <td>0.492065</td>
      <td>0.227002</td>
      <td>0.975453</td>
    </tr>
    <tr>
      <th>3.0</th>
      <td>0.259839</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>0.481183</td>
      <td>0.227002</td>
      <td>1.144654</td>
    </tr>
  </tbody>
</table>
</div>




```python
yield_cost_efficiency_df = pd.read_csv(efficiency_file)
fig=draw_trade_off(yield_cost_efficiency_df,trade_off_enzyme_efficiency_figfile)
fig.show()
```


<div>                            <div id="c1f74121-59ac-46d6-ae99-c53bec99df09" class="plotly-graph-div" style="height:600px; width:800px;"></div>            <script type="text/javascript">                require(["plotly"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("c1f74121-59ac-46d6-ae99-c53bec99df09")) {                    Plotly.newPlot(                        "c1f74121-59ac-46d6-ae99-c53bec99df09",                        [{"line":{"color":"red","dash":"dash","width":3},"marker":{"color":"red","size":10},"mode":"lines","name":"biomass yield","x":[1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5],"y":[0.4974250905218673,0.4974250905218655,0.4974250905218661,0.4920652548506863,0.4811827774021628,0.4579862429068541,0.4375262337187864,0.4202384916354786,0.3917926089897021,0.368415628292455,0.3489348110447498,0.3324510426043874,0.3157439804733918,0.3006248671742367,0.285896682733143,0.2729012258733541,0.2613497086646535,0.2510141327558489],"type":"scatter"},{"line":{"color":"blue","width":3},"marker":{"color":"blue","size":10,"symbol":5},"mode":"lines","name":"enzyme efficiency","x":[1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5],"xaxis":"x2","y":[0.8974084915316362,0.897408491531657,0.8974084915316339,0.9754529534931816,1.1446539788676235,1.271051694033473,1.3877093637245963,1.4994834093784677,1.553313372635864,1.606693697397809,1.66007390364879,1.7134539913892135,1.7525247427127746,1.7877915660956418,1.8135502340394103,1.839308874408199,1.865067487202114,1.890826013065673],"yaxis":"y2","type":"scatter"}],                        {"height":600,"legend":{"font":{"color":"black","family":"Times New Roman","size":20},"x":0.07,"y":0.95},"plot_bgcolor":"white","width":800,"xaxis":{"linecolor":"black","range":[1,9.6],"tickcolor":"black","tickfont":{"color":"black","family":"Times New Roman","size":15},"ticks":"inside","title":{"font":{"family":"Times New Roman","size":20},"text":"\u003cb\u003eSubstrate uptake rate (mmol\u002fgDW\u002fh)\u003cb\u003e"}},"yaxis":{"linecolor":"black","range":[0.15101413275584888,0.5974250905218673],"tickcolor":"black","tickfont":{"color":"black","family":"Times New Roman","size":15},"ticks":"inside","title":{"font":{"family":"Times New Roman","size":20},"text":"\u003cb\u003eBiomass yield (gDW\u002fg glucose)\u003cb\u003e"}},"xaxis2":{"linecolor":"black","overlaying":"x","range":[1,9.6],"showticklabels":false,"side":"top"},"yaxis2":{"linecolor":"black","overlaying":"y","range":[0.797408491531634,1.9908260130656732],"side":"right","tickcolor":"black","tickfont":{"color":"black","family":"Times New Roman","size":15},"ticks":"inside","title":{"font":{"family":"Times New Roman","size":20},"text":"\u003cb\u003eEnzyme efficiency (gDW\u002fg enzyme)\u003cb\u003e"}},"template":{"data":{"histogram2dcontour":[{"type":"histogram2dcontour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"choropleth":[{"type":"choropleth","colorbar":{"outlinewidth":0,"ticks":""}}],"histogram2d":[{"type":"histogram2d","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmap":[{"type":"heatmap","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmapgl":[{"type":"heatmapgl","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"contourcarpet":[{"type":"contourcarpet","colorbar":{"outlinewidth":0,"ticks":""}}],"contour":[{"type":"contour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"surface":[{"type":"surface","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"mesh3d":[{"type":"mesh3d","colorbar":{"outlinewidth":0,"ticks":""}}],"scatter":[{"fillpattern":{"fillmode":"overlay","size":10,"solidity":0.2},"type":"scatter"}],"parcoords":[{"type":"parcoords","line":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolargl":[{"type":"scatterpolargl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"scattergeo":[{"type":"scattergeo","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolar":[{"type":"scatterpolar","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"scattergl":[{"type":"scattergl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatter3d":[{"type":"scatter3d","line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattermapbox":[{"type":"scattermapbox","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterternary":[{"type":"scatterternary","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattercarpet":[{"type":"scattercarpet","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"type":"carpet"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}],"barpolar":[{"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"pie":[{"automargin":true,"type":"pie"}]},"layout":{"autotypenumbers":"strict","colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"hovermode":"closest","hoverlabel":{"align":"left"},"paper_bgcolor":"white","plot_bgcolor":"#E5ECF6","polar":{"bgcolor":"#E5ECF6","angularaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"radialaxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"ternary":{"bgcolor":"#E5ECF6","aaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"baxis":{"gridcolor":"white","linecolor":"white","ticks":""},"caxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]]},"xaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"yaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"scene":{"xaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"yaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"zaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"geo":{"bgcolor":"white","landcolor":"#E5ECF6","subunitcolor":"white","showland":true,"showlakes":true,"lakecolor":"white"},"title":{"x":0.05},"mapbox":{"style":"light"}}}},                        {"responsive": true}                    ).then(function(){

var gd = document.getElementById('c1f74121-59ac-46d6-ae99-c53bec99df09');
var x = new MutationObserver(function (mutations, observer) {{
        var display = window.getComputedStyle(gd).display;
        if (!display || display === 'none') {{
            console.log([gd, 'removed!']);
            Plotly.purge(gd);
            observer.disconnect();
        }}
}});

// Listen for the removal of the full notebook cells
var notebookContainer = gd.closest('#notebook-container');
if (notebookContainer) {{
    x.observe(notebookContainer, {childList: true});
}}

// Listen for the clearing of the current output cell
var outputEl = gd.closest('.output');
if (outputEl) {{
    x.observe(outputEl, {childList: true});
}}

                        })                };                });            </script>        </div>

