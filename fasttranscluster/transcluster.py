import numpy as np
from scipy.special import gammaln
from scipy.stats import poisson
from numba import jit


def memoize(f):
    memo = {}

    def helper(N, K, delta, lamb, beta, lgamma):
        x = (N, K, delta, lamb, beta)
        if x not in memo:
            memo[x] = f(N, K, delta, lamb, beta, lgamma)
        return memo[x]

    return helper


def calculate_trans_prob(sparse_snp_dist, sample_dates, K, lamb, beta):
    # precalculate lgamma
    max_nk = max([t[2] for t in sparse_snp_dist]) + K
    lgamma = gammaln(np.arange(max_nk + 2))

    row_ind = []
    col_ind = []
    lprob = []
    for i, j, d in sparse_snp_dist:
        delta = np.abs(sample_dates[i][1] - sample_dates[j][1])
        lp = lprob_transmission(d, K, delta, lamb, beta, lgamma)
        row_ind.append(i)
        col_ind.append(j)
        lprob.append(lp)

    return row_ind, col_ind, lprob


@memoize
def lprob_transmission(N, K, delta, lamb, beta, lgamma):
    lprob = -np.inf
    for k in range(K):
        lprob = np.logaddexp(lprob,
                             lprob_k_given_N(N, k, delta, lamb, beta, lgamma))
    return lprob


@memoize
@jit(nopython=True)
def lprob_k_given_N(N, k, delta, lamb, beta, lgamma):
    lprob = (N + 1) * np.log(lamb) - delta * (
        lamb + beta) + k * np.log(beta) - lgamma[k + 1]

    # ugly poisson cdf but allows for use of numba
    if delta > 0:
        pois_cdf = -lamb * delta
        for i in range(N + 1):
            pois_cdf = np.logaddexp(i * np.log(lamb * delta) - lgamma[i + 1],
                                    pois_cdf)
        lprob -= pois_cdf

    integral = -lgamma[N + 1]
    for i in range(N + k + 1):
        if delta > 0:
            integral = np.logaddexp(
                lgamma[N + k + 1] - lgamma[i + 1] - lgamma[N + k - i + 1] +
                (N + k - i) * np.log(delta) + lgamma[i + 1] -
                (i + 1) * np.log(lamb + beta), integral)
        else:
            integral = np.logaddexp(
                lgamma[N + k + 1] - lgamma[i + 1] - lgamma[N + k - i + 1] +
                lgamma[i + 1] - (i + 1) * np.log(lamb + beta), integral)
    lprob += integral
    return (lprob)
