# Finds the paramters of the dirichlet multinomial distribution
# by maximising the leave-one-out (LOO) likelihood or fixed point iteration as described in
# Minka, 2000. https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf

import numpy as np
from scipy.special import psi


def find_dirichlet_priors(data, min_cov=10, max_iter=1000, tol=1e-5, method='LOO'):
    K = np.shape(data)[1]
    data = data[np.count_nonzero(data, 1)>1, ]
    data.sort(axis=1)
    total_counts = np.sum(data, 1)

    alpha=np.mean(data, 0) + 0.5
    nalpha = np.zeros(K)

    if method=='LOO':
        for i in range(max_iter):
            a0 = np.sum(alpha)
            for k in range(K):
                nalpha[k] = alpha[k] * np.sum(data[:,k] / (data[:,k] - 1 + alpha[k]), 0) / np.sum(total_counts/ (total_counts - 1 + a0), 0)
            if np.max(np.abs(nalpha - alpha)) < tol:
                alpha = nalpha
                break
            alpha = nalpha
    else:
        for i in range(max_iter):
            a0 = np.sum(alpha)
            for k in range(K):
                nalpha[k] = (alpha[k] * 
                    np.sum(psi(data[:,k] + alpha[k]) - psi(alpha[k])) / 
                    np.sum(psi(total_counts + a0) - psi(a0)))
            if np.max(np.abs(nalpha - alpha)) < tol: 
                alpha = nalpha
                break
            alpha = nalpha

    alpha[::-1].sort()
    print("Calculated alphas: ", alpha)
    
    return (alpha)