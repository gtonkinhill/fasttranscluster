# Finds the paramters of the dirichlet multinomial distribution
# by maximising the leave-one-out (LOO) likelihood as described in
# Minka, 2000. https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf

import numpy as np

def find_dirichlet_priors(data, alpha=np.full(4, 0.25), max_iter=1000, tol=1e-3):
    total_counts = np.sum(data, 1)
    for i in range(max_iter):
        nalpha = alpha * np.sum(data / (data - 1 + alpha), 0) / np.sum(total_counts/ (total_counts - 1 + np.sum(alpha)), 0)
        print(np.max(np.abs(nalpha - alpha)))
        if np.max(np.abs(nalpha - alpha)) < tol: return(nalpha)
        alpha = nalpha
    print("Reached maximum number of iterations...")
    return (alpha)
