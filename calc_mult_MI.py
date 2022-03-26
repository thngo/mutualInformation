import numpy as np
import pandas as pd
import statistics as st
import re
import csv
import scanpy as sc 
import scanpy.external as sce
import phate
import matplotlib
from matplotlib import pyplot as plt
import warnings
from scipy import stats
from scipy.stats import binom
from scipy.stats import multinomial
import seaborn as sns
from scipy.stats import hypergeom
import warnings
warnings.filterwarnings('ignore')
import pickle
from random import sample
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
import scipy.io as scio
import itertools

# from TN's derivation of multivariate MI

def calc_mult_MI(M,bins):
    totLen = len(M)
    ## Getting probabilities ################
    c = list()
    # getting the histogram distributions 
       
    for i in np.arange(totLen):
        # Get combination of probabilities 
        tmp = list(range(0, totLen)) # range 0 to totLen because python start indexing at 0
        combList = list(itertools.combinations(tmp, i+1)) # get number of iterations
        # Calculate probabilities, for each combinations
        ind_prob = list() # individual probabilities at each bin
        termList = list()
        for j in np.arange(len(combList)):
            # Create subset
            M_sub = list() # initialize subset
            for iSub in combList[j]:
                M_sub.append(M[iSub])
            # Frequency observed in each bin of this M_sub sample
            ind_prob.append(probNorm(M_sub,bins))
        # Calculate each log term
            # i.e. (-1)^1*log(P(X1)*P(X2)*P(X3)*...*P(Xn)) for individual probabilities
            # Sign is (-) if number of variables is ODD when we calculate the probabilities
            # i.e. (-1)^2*log(P(X1,X2)*P(X1,X3)*P(X2,X3)) for joint probabilities
            # Sign is (+) if number of variables is EVEN when we calculate the probabilities
        termList.append(pow(-1,j+1)*np.log(ind_prob))
        # termList.append(pow(-1,j+1)*np.log(np.nanprod(ind_prob))) # power to (j+1) because python start indexing at 0
        
    # Calculate MI
    MI = probNorm(M,bins) * np.nanprod(termList)
    return MI

def probNorm(M,bins):
    c = np.histogramdd(M,bins)[0]
    #normalzing the data to get frequency
    c_normalized = c / float(np.sum(c))
    #c_normalized = c_normalized[np.nonzero(c_normalized)]
    return c_normalized
