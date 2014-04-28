"""
Utility functions
"""
from math import log
from utils import simplex_sample,inverse_cdf_sample

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
    """compute entropy (in bits) of a probability distribution ps"""
    return -sum([plogp(p) for p in ps])*unit_factor(units)
