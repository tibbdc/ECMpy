# 1.Model Preview
## Import related functions


```python
#from script.ECMpy_function import *
import sys
sys.path.append(r'./script/')
from ECMpy_function import *
```


```python
model_file = "./data/iML1515_new.xml"
bigg_met_file = './data/bigg_models_metabolites.txt'

Determine_suitable_ecGEM(model_file,bigg_met_file)
```




    'Suitable for constructing enzyme-bound models.'


