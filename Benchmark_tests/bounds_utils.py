from qubovert import QUBO, PUBO

def only_positive_coefficient(qubo):
    """function for estimating the upper and lower bound

    Keyword arguments:
    qubo -- is the qubo model of the main problem cost function


    Return values:
    lower bound -- is the function lower bound
    upper bound -- is the function upper bound
    or return false if not all the involve coefficients are negative

    """
    lowerBound = 0.0
    upperBound = 0.0
    for key in qubo:
        if len(key) > 0:
            if qubo[key] > 0:
                upperBound += qubo[key]
            elif qubo[key] < 0:
                return False
        else:
            # offset
            upperBound += qubo[key]
            lowerBound += qubo[key]
    return lowerBound, upperBound

def bounds_naive(qubo):
    """function for estimating the upper and lower bound

    Keyword arguments:
    qubo -- is the qubo model of the main problem cost function


    Return values:
    lower bound -- is the function lower bound
    upper bound -- is the function upper bound

    """
    upperBound = 0.0
    lowerBound = 0.0
    for key in qubo:
        if len(key) > 0:
            if qubo[key] > 0:
                upperBound += qubo[key]
            else:
                lowerBound += qubo[key]
        else:
            # offset
            upperBound += qubo[key]
            lowerBound += qubo[key]
    return lowerBound, upperBound

def bounds_pos_neg(qubo):
    """function for estimating the upper and lower bound

    Keyword arguments:
    qubo -- is the qubo model of the main problem cost function

    Return values:
    lower bound -- is the function lower bound
    upper bound -- is the function upper bound

    """
    p_sum = {}
    n_sum = {}
    offset = 0
    for var in qubo.variables:
        p_sum[var] = 0
        n_sum[var] = 0
    for key in qubo:
        if len(key) == 1:
            p_sum[key[0]] += qubo[key]
            n_sum[key[0]] += qubo[key]
        elif len(key) > 1:
            if qubo[key] < 0:
                for elem in key:
                    p_sum[elem] += qubo[key]
            elif qubo[key] > 0:
                for elem in key:
                    n_sum[elem] += qubo[key]
        else:
            offset = qubo[key]
    lowerBound = offset
    upperBound = offset
    for key in p_sum:
        if p_sum[key] < 0:
            lowerBound += p_sum[key]
        if n_sum[key] > 0:
            upperBound += n_sum[key]
    return lowerBound, upperBound
