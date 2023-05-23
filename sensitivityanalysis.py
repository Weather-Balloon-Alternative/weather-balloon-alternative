import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_excel("data2.xlsx")
names = df.columns.to_list()[1:]
# names = ['Powered A/C',
#          f'Glider &\nStowed Balloon',
#          f'Controlled &\nDeflatable\nBalloon',
#          f'Deflatable\nBalloon &\nParasail',
#          f'Glider &\nDisposable\nBalloon'] # Manual override for the first trade off
data = df.to_numpy()



def monte_carlo(n_iter, n_options, data, names, nchanged, seed=11766343):
    '''
    Args:
        n_iter: int, number of runs
        n_options: int, number of design options
        data: n x n_options array, contains rankings for each concepts for all criteria
        seed: int, random seed generator, must stay the same from run to run

    Returns:
        results: n_iter x n_options array, accumulation of resulting scores for all design options
        wins: 1 x n_options array, tally of wins for each design option
    '''

    results = np.zeros((n_iter, n_options))
    wins = np.zeros(n_options)
    np.random.seed(seed=seed)

    Scores = data[:, 1:]  # Set equal to imported scores from pandas
    Weights_init = data[:, 0]  # Set equal to imported weights from pandas

    weight_min = min(Weights_init)
    weight_max = max(Weights_init)

    for i in range(n_iter):
        # Select which weight to change
        weightselect = np.array([int(np.ceil(np.random.uniform(low=-1, high=len(Weights_init)-1))) for i in range(nchanged)])
        # Select which weight to change it to
        weightval = np.array([int(np.ceil(np.random.uniform(low=weight_min-1, high=weight_max-1))) for i in range(nchanged)])

        # Change weight
        Weights = np.array(Weights_init)
        Weights[weightselect] = weightval

        # Calculate results
        results[i] = np.mat(Weights) * np.mat(Scores) / sum(Weights)
        wins[np.where(results[i] == np.max(results[i]))[0][0]] += 1

    return results, wins


def plotting(results, wins, names, n_iter):
    if not all(isinstance(item,str) for item in names):
        raise ValueError('Names must be of type string')

    for i in range(len(names)):
        names[i] += f'\n(average: {                   np.round(np.average(results[:,i]),2)})'


    plt.bar(names,wins*100/n_iter)
    plt.ylim((0,100))
    plt.ylabel('Winning [%]')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    results_dict = {}
    for i in range(len(names)):
        results_dict[names[i]] = results[:,i]

    results_df = pd.DataFrame(results_dict)
    results_df.plot(kind='box')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    n_iter = int(1e6)
    n_changed = 5 # Amount of weights to be changed in a single run
    results, wins = monte_carlo(n_iter, len(names), data, names, n_changed)
    for i in range(len(names)):
        print(f'{names[i]} average: {np.average(results[:, i])} win %: {wins[i] * 100 / n_iter}')
    plotting(results, wins, names, n_iter)
