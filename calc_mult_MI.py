#from TN's derivation of multivariate MI

def calc_mult_MI(M,bins):
    totLen = len(M)
    ## Getting probabilities ################
    c = list()
    #getting the histogram distributions 
    
    # individual probabilities
    
        
    for i in totLen:
        # Get combination of probabilities 
        tmp = list(range(0, totLen)) # range 0 to totLen because python start indexing at 0
        combList = list(itertools.combinations(tmp, i)) # get number of iterations
        # Calculate probabilities, for each combinations
        ind_prob = list()
        termList = list()
        for j in len(combList):
            # Create subset
            M_sub = list() # initialize subset
            for iSub in test[k]:
                M_sub.append(M[iSub])
            ind_prob.append(probNorm(M_sub,bins))
        # Calculate each log term
            # i.e. (-1)^0*log(P(X1)*P(X2)*P(X3)*...*P(Xn)) for individual probabilities
            # Sign is (+) if number of variables is ODD when we calculate the probabilities
            # i.e. (-1)^1*log(P(X1,X2)*P(X1,X3)*P(X2,X3)) for joint probabilities
            # Sign is (-) if number of variables is EVEN when we calculate the probabilities
        termList.append(pow(-1,j)*np.log(ny.nanprod(ind_prob)))
        
    # Calculate MI
    MI = prob_TN(M,bins) * np.nanprod(termList)
    return MI

def probNorm(M,bins):
    c = np.histogramdd(M,bins)[0]
    #normalzing the data to get frequency
    c_normalized = c / float(np.sum(c))
    #c_normalized = c_normalized[np.nonzero(c_normalized)]
    return c_normalized
