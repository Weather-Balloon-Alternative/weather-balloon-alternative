import numpy as np
import pandas as pd
# Import pandas stuff here
data = pd.read_excel("data.xlsx", index_col=0)
data = data.to_numpy()
weights = data[:,0]
data = data[:,1:]
print(data)
print(weights)


Scores = np.nan     # Set equal to imported scores from pandas
Weights_init = np.nan # et equal to imported weights from pandas

n_it = 10
weight_min = 1
weight_max = 20

results = np.zeros((n_it,4))

for i in range(n_it):
    # Select which weight to change
    weightselect = int(np.ceil(np.random.uniform(low=-1,high=len(Weights_init)-1)))
    # Select which weight to change it to
    weightval = int(np.ceil(np.random.uniform(low=weight_min-1,high=weight_max-1)))
    
    # Change weight
    Weights = Weights_init
    Weights[weightselect] = weightval
    # Calculate results
    results[i] = np.mat()