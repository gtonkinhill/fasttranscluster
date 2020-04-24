import numpy as np
from scipy.special import gammaln
from scipy.stats import poisson

def calculate_trans_prob(sparse_snp_dist, sample_dates, K, lamb, beta):
    trans_dist = []
    for i,j,d in sparse_snp_dist:
        delta = np.abs(sample_dates[i] - sample_dates[j])
        lp = lprob_transmission(d, K, delta, lamb, beta)
        trans_dist.append((i,j,lp))
    
    return(lp)


def lprob_transmission(N, K, delta, lamb, beta):
    lprob = 0
    for k in range(K):
        lprob = np.logaddexp(lprob,
            lprob_k_given_N(N, k, delta, lamb, beta)
        )
    
    print(N, lprob)

    return lprob


def lprob_k_given_N(N, k, delta, lamb, beta):
    lprob = (N+1)*lamb-delta*(lamb+beta)+k*beta-gammaln(k+1)-poisson.logcdf(N, lamb*delta)

    print(N, k, delta, lamb, beta)
    print(poisson.logcdf(N, lamb*delta), (N+1)*lamb-delta*(lamb+beta), k*beta-gammaln(k+1))

    integral = 0
    for i in range(N):
        for j in range(k):
            integral = np.logaddexp((N-i+k-j)*np.log(delta) +
                gammaln(k+1) -
                gammaln(i+1) -
                gammaln(N-i+1) -
                gammaln(j+1) -
                gammaln(k-j+1) +
                gammaln(j+i+1) -
                k*np.log(lamb+beta),
                integral)
    
    print(lprob, integral)
    lprob += integral

    return(lprob)