# Import pandas stuff here

import numpy as np

Scores = np.nan     # Set equal to imported scores from pandas
Weights_init = np.nan # 

n_it = 1E6
weight_min = 1
weight_max = 20

results = np.zeros((4,n_it))

for i in range(n_it):
    # Select which weight to change
    weightselect = np.ceil(np.random.uniform(low=(weight_min-1),high=weight_max))