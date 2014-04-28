from nsb import nsb
from h_mle import h_mle
from entropy_estimation_utils import *

def test_nsb(alphabet_size=4,sample_size=1000):
    ps = simplex_sample(alphabet_size)
    true_h = h(ps,units='nats')
    xs = sample_ps(ps,sample_size)
    print true_h,h_mle(xs,units='nats'),h_nsb(xs,alphabet_size=alphabet_size)
