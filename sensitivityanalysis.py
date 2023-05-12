# Import pandas stuff here

import numpy as np

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