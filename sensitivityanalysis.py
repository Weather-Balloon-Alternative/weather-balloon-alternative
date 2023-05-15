import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Import pandas stuff here
data = pd.read_excel("data.xlsx", index_col=0)
data = data.to_numpy()
weights = data[:,0]
data = data[:,1:]

Scores = data     # Set equal to imported scores from pandas
Weights_init = weights # et equal to imported weights from pandas

n_it = 10**6
weight_min = 1
weight_max = 20

results = np.zeros((n_it,5))
wins = np.zeros(5)

np.random.seed(seed=11766343)

for i in range(n_it):
    # Select which weight to change
    weightselect = int(np.ceil(np.random.uniform(low=-1,high=len(Weights_init)-1)))
    # Select which weight to change it to
    weightval = int(np.ceil(np.random.uniform(low=weight_min-1,high=weight_max-1)))
    
    # Change weight
    Weights = Weights_init
    Weights[weightselect] = weightval
    Weights = Weights/sum(Weights)
    # Calculate results
    results[i] = np.mat(Weights) * np.mat(Scores)
    wins[np.where(results[i] == np.max(results[i]))[0][0]] += 1


print(f'Option 1 average: {np.average(results[:,0])} win %: {wins[0]*100/n_it}')
print(f'Option 2 average: {np.average(results[:,1])} win %: {wins[1]*100/n_it}')
print(f'Option 3 average: {np.average(results[:,2])} win %: {wins[2]*100/n_it}')
print(f'Option 4 average: {np.average(results[:,3])} win %: {wins[3]*100/n_it}')
print(f'Option 5 average: {np.average(results[:,4])} win %: {wins[4]*100/n_it}')



names = np.array([f'Powered A/C\n(average: {                   np.round(np.average(results[:,0]),2)})',
                  f'Glider & Stowed \nBalloon\n(average: {         np.round(np.average(results[:,1]),2)})',
                  f'Controlled Deflatable \nBalloon\n(average: {    np.round(np.average(results[:,2]),2)})',
                  f'Deflatable Balloon & \nParasail\n(average: {              np.round(np.average(results[:,3]),2)})',
                  f'Glider & Disposable \nBalloon\n(average: {       np.round(np.average(results[:,4]),2)})'])

plt.bar(names,wins*100/n_it)
plt.ylim((0,100))
plt.ylabel('Winning [%]')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()

results_dict = {}
for i in range(5):
    results_dict[names[i]] = results[:,i]

results_df = pd.DataFrame(results_dict)
results_df.plot(kind='box')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()