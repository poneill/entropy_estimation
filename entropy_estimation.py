from entropy_estimation_utils import h
from utils import frequencies

def h_mle(xs,units='bits'):
    """Compute MLE estimator of entropy"""
    p_hats = frequencies(xs)
    return h(p_hats,units=units)

