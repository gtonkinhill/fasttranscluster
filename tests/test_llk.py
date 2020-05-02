# test llk gives the same answer as found using Sage math cloud.
from fasttranscluster.transcluster import lprob_transmission
from fasttranscluster.transcluster import lprob_k_given_N
from scipy.special import gammaln

def test_lprob_k_given_N():
    # compare to the following code in sage
    # from sage.symbolic.integration.integral import definite_integral
    # N=7
    # k=4
    # d=0.16963
    # y=3
    # b=52
    # var('h,j,i')
    # assume(j,'integer')
    # assume(j>0)
    # B = definite_integral(e^(-h*(y+b)) * ((h+d)^k )* sum(h^j *(d^(N-j)/(factorial(j)*factorial(N-j))), j, 0, N) , h, 0,  +Infinity)
    # B = B * (e^(-d*(y+b))*y^(N+1)*b^k)/(factorial(k)*sum((e^(-y*d)*(y*d)^i)/factorial(i), i, 0, N))
    # log(B)

    N=7
    k=4
    delta=0.16963
    lamb=3
    beta=52
    lgamma = gammaln(range(20))

    assert lprob_k_given_N(N, k, delta, lamb, beta, lgamma) == -17.9565184209608

    assert lprob_k_given_N(N, k, 0, lamb, beta, lgamma) + 17.6950323846585 < 1e-6

    return