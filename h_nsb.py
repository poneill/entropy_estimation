"""
Implements the NSB estimator.  See Nemenmen, Shafee and Bialek, 2002.
"""

from collections import Counter
from entropy_estimation_utils import h,integrate_improp
from utils import normalize
from scipy.special import polygamma
from math import log

def h_nsb(xs,alphabet_size,int_points=1000):
    """Compute the Nemenman-Shafree-Bialek Estimator"""
    ns = ns_from_xs(xs,alphabet_size)
    laplace_pseudocount = [1]*alphabet_size
    dalpha = dda_expected_entropy(laplace_pseudocount)
    f = lambda a:expected_entropy_from_alphas_ref([a + n for n in ns]) * dalpha(a) #ns - 1?
    return integrate_improp(f,int_points)

def dda_expected_entropy(qs):
    """return d/da[E[H|a*qs]], a function of alpha"""
    # Agrees with test_diff!
    sum_qs = float(sum(qs))
    h_inf = h(normalize(qs),units='nats')
    h_0 = 0 #expected_entropy_from_alphas([0 for q in qs])
    Z = h_inf#*log(2) # in nats
    return lambda alpha: ((sum_qs*polygamma(1,alpha*sum_qs+1) -
                           sum(qj**2/sum_qs*polygamma(1,alpha*qj+1) for qj in qs))/
                          (Z))

