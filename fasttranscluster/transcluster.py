import numpy as np
from scipy.special import gammaln
from scipy.stats import poisson
import functools



def calculate_trans_prob(sparse_snp_dist, sample_dates, K, lamb, beta):
    row_ind = []
    col_ind = []
    lprob = []
    for i,j,d in sparse_snp_dist:
        delta = np.abs(sample_dates[i][1] - sample_dates[j][1])
        lp = lprob_transmission(d, K, delta, lamb, beta)
        row_ind.append(i)
        col_ind.append(j)
        lprob.append(lp)
    
    return row_ind, col_ind, lprob

@functools.lru_cache(maxsize=None)
def lprob_transmission(N, K, delta, lamb, beta):
    lprob = -np.inf
    for k in range(K):
        lprob = np.logaddexp(lprob,
            lprob_k_given_N(N, k, delta, lamb, beta)
        )
    
    return lprob

@functools.lru_cache(maxsize=None)
def lprob_k_given_N(N, k, delta, lamb, beta):
    lprob = (N+1)*np.log(lamb)-delta*(lamb+beta)+k*np.log(beta)-gammaln(k+1)-poisson.logcdf(N, lamb*delta)

    integral = -np.inf
    for i in range(N+k+1):
        integral = np.logaddexp(
            gammaln(N+k+1) - 
            gammaln(i+1) - 
            gammaln(N+k-i+1) + 
            (N+k-i)*np.log(delta) + 
            gammaln(i+1) - 
            (i+1)*np.log(lamb+beta),
            integral
        )
    integral += -gammaln(N+1)
    lprob += integral

    return(lprob)