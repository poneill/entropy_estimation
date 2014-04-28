"""
Utility functions
"""
from math import log
from utils import simplex_sample,inverse_cdf_sample,mean,interpolate

def sample_ps(ps,N):
    """Return a sample from the multinomial distribution given by ps"""
    ks = range(N)
    return [inverse_cdf_sample(ks,ps) for i in xrange(N)]

def plogp(x):
    """Compute lim p-> x p*log(p)"""
    if x > 0:
        return x*log(x)
    else:
        return 0

def unit_factor(units):
    """Multiply entropy in nats by this factor to obtain entropy in desired units"""
    if units == 'nats':
        return 1
    elif units == 'bits':
        return log(2)
    elif units == 'dits':
        return log(10)
    else:
        raise Exception("Unknown units of entropy:",units)
    
def h(ps,units='bits'):
    """compute entropy (in specified units) of a probability distribution ps"""
    return -sum([plogp(p) for p in ps])*unit_factor(units)

def integrate_improp(f,int_points=1000):
    """Integrate f from a to infinity using transform u = 1/(1-t)"""
    interp_start = 0
    interp_stop = 1 - 1.0/int_points
    return mean(f(t/(1-t))*1/(1-t)**2
                for t in (interpolate(interp_start,interp_stop,int_points)))

def expected_entropy_from_alphas_ref(alphas):
    """Compute expectation of entropy (in nats), given alphas.  Eq (9) of AMP"""
    if sum(alphas) == 0:
        return 0
    kappa = float(sum(alphas))
    return (polygamma(0,kappa + 1) - sum(a/kappa*polygamma(0,a+1)
                                        for a in alphas))

def ns_from_xs(xs,alphabet_size):
    counts = Counter(xs)
    ns = [counts[b] for b in counts.keys()] + [0]*(alphabet_size-len(counts))
    return ns
