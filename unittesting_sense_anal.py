import pytest
import numpy as np
from sensitivityanalysis import monte_carlo, data

'''def test_convergence():
    n_iter1, n_iter2 = int(1e2), int(1e5)
    results1, wins1 = monte_carlo(n_iter1, 5, data)
    results2, wins2 = monte_carlo(n_iter2, 5, data)
    for n in range(5):
        # test to check that resultant averages are less than 5% different between low and high number of iterations
        assert (abs(np.average(results1[:, n]) - np.average(results2[:, n])) / np.average(results2[:, n])) < 0.05'''

def test_equality():
    data[:, 1:] = np.ones((np.shape(data)[0], np.shape(data)[1]-1))
    print(data)
    results, wins = monte_carlo(int(1e2), 5, data)
    assert np.array_equal(results[:, 0].T, results[:, 1].T, results[:, 2].T) == True


