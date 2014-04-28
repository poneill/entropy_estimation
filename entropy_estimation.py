from math import log
from collections import Counter
from utils import frequencies

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

def mle_h(xs,units='bits'):
    """Compute MLE estimator of entropy"""
    p_hats = frequencies(xs)
    return h(p_hats,units=units)
    
    
