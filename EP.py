import numpy as np
from numpy import exp
from numpy.linalg import inv

class EP(object):
    # variables
    constellation = None;
    constellation_len = 0;
    es = 1;                                 # constellation average power
    beta = 0.2;                             # the percentage for taking values from the previous iteration 
    epsilon = 2*5e-7;                       # the default minimal variance is 2*5e-7 (5e-7 for real only)
    l = 10;                                 # maximal iteration
    early_stop = False;
    early_stop_min_diff = 2e-4;             # the elemental minimal difference for mean and variance in Bayesian estimation (1e-4 for real only)
    
    '''
    init
    @constellation:             the constellation, a vector
    @beta:                      the percentage for taking values from the previous iteration
    @epsilon:                   the default minimal variance
    @l:                         maximal iteration
    @early_stop:                whether stop early
    @early_stop_min_diff:       the elemental minimal difference in mean and variance
    '''
    def __init__(self, constellation, *, beta=None, epsilon=None, l=None, early_stop=None, early_stop_min_diff=None):
        # inputs
        constellation = np.asarray(constellation).squeeze();
        if constellation.ndim != 1:
            raise Exception("The constellation must be a vector.");
        else:
            self.constellation = constellation;
            self.constellation_len = len(constellation);
        # optional iputs
        if beta is not None:
            self.beta = beta;
        if epsilon is not None:
            self.epsilon = epsilon;
        if l is not None:
            self.l = l;
        if early_stop is not None:
            self.early_stop = early_stop;
        if early_stop_min_diff is not None:
            self.early_stop_min_diff = early_stop_min_diff;
        