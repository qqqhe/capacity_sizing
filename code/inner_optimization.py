# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 15:20:42 2018

Author: Qie He
"""

import numpy as np
    
def int2bin(k, n):
    """Transform a non-negative integer k into an array of n elements that represents its binary expansion"""
    binary_expansion = np.zeros(n, dtype=int)
    position = n-1    
    while k > 0:
        if k % 2 == 1:            
            binary_expansion[position] = 1
        k = int(k/2)
        position -=1
    return binary_expansion
    
def getx(v, lb, ub, i, B):
    """Return a vector x such that 
    (1) all the components of x sum up to 1
    (2) for each $j != i$, x[j]=lb[j] if v[j]=0 and x[j]=ub[j] if v[j]=1
    (3) lb[i] <= x[i] <= ub[i]
    """
    x = lb + np.multiply((ub - lb), v)
    x[i] = B - (x.sum() - x[i])
    # Test if variable x[i] is within the bounds
    if x[i] <= ub[i] and x[i] >= lb[i]:
        return x
    else:
        return np.array([])


def max_sum_xlogx(n, B, lb, ub):
    """
    Find the optimal solution of max x1*log(x1) + ... + xn*log(xn)
    subject to constraints x1 + ... xn = B, lb_i <= xi <= ub_i.
    
    Parameters
    ----------
    n: int
        Number of variables.
    B: float
        Total number of resources.
    lb: array of floats
        Lower bounds of variables.
    ub: array of floats
        Upper bounds of variables.
        
    Returns
    -------    
    opt_obj: float
        Optimal objective value.
    opt_sol: an array of floats
        Optimal solution.
    """    
    # Initialize the optimal solution and optimal objective
    opt_obj = - n/np.e
    opt_sol = np.array([])
    
    # First select a variable whose value may not be at bound, indexed by idx-var_interior
    for idx_var_interior in range(n):
        for idx in range(pow(2, n-1)):
            idx_binary_expansion = int2bin(idx, n-1)
            # Insert element 0 into position idx_var_interior
            idx_binary_expansion = np.insert(idx_binary_expansion, idx_var_interior, 0)
            # Compute the solution x with all (but one) variables at bounds
            x = getx(idx_binary_expansion, lb, ub, idx_var_interior, B)
            
            if x.size > 0:
                obj = np.multiply(x, np.log(x)).sum()
                if obj >= opt_obj:
                    opt_obj = obj
                    opt_sol = x
    return [opt_obj, opt_sol]
                                    

# Additional test data
# lowerbound = np.array([0.5, 0.47, 0.45, 0.43, 0.41, 0.39, 0.37])
# upperbound = np.array([0.66, 0.70, 0.71, 0.73, 0.84, 0.88, 0.91])
# res_bound = 4.2
#num_var = len(lowerbound)

#num_var = 10
#a1 = np.random.rand(num_var)
#a2 = np.random.rand(num_var)
#lowerbound = np.minimum(a1, a2)
#upperbound = np.maximum(a1, a2)
#res_bound = lowerbound.sum() + np.random.rand()*(upperbound.sum() - lowerbound.sum())  # total resources available
#
#[opt_obj, opt_sol] = max_sum_xlogx(num_var, res_bound, lowerbound, upperbound)
#print("The optimal solution is:\n", opt_sol)
#print("The optimal objective is:\n", opt_obj)