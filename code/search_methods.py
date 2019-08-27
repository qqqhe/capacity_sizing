# -*- coding: utf-8 -*-
"""
Created on Tue May 22 15:42:06 2018

Author: Qie He
"""


def bisection_search(fun, lb, ub, threshold, mode, *fun_args):
    '''
    A bi-section search method to find a root of a monotone function within a given range.
    
    Parameters
    ----------
    fun: function name
    lb : float
        Lower bound of the search range.
    ub : float
        Upper bound of the search range.
    threshold: list of floats
        The first element is the tolerance for checking if the function value is close to zero.
        The second element is to check if the width of the search range is close to zero.
    mode: boolean
        True if monotonically nondecreasing, False if monotonically nonincreasing.    
    *fun_args: list of arguments for the function
        For the zeta functions, the list is [network, totalflow, route_idx, fun_idx]
        For the zeta_homo_noparking function, the list is [network]
        For the diff_totalflow_sumflow function, the list is [network]      

    Returns
    ----------
    x: float
        The root of the function.
        
    
    
    
# lb: lower bound
# ub: upper bound
# threshold: tolerance for stopping the search
# mode = true: monotone increasing function; decreasing otherwise
# fun: the function for which you need to find a root
#
# *fun_args: additional variables defining the function:
# (1) fun_args = [1, total_flow, route]
# (2) fun_args = [2, network]
# (3) fun_args = [3, total_flow, network, function_index]
    '''
    x = (lb + ub) / 2     
   
    network = fun_args[0]
    # For the zeta function 
    if len(fun_args) > 1:               
        total_flow = fun_args[1]
        route_idx = fun_args[2]    
        fun_idx = fun_args[3]    
        fun_val = fun(network, route_idx, x, total_flow, fun_idx)
    else:
        fun_val = fun(network, x)
    
    if mode:
        while (abs(fun_val) > threshold[0]) or ((ub - lb) > threshold[1]):
            if fun_val > 0:
                ub = x
            else:
                lb = x
            x = (lb + ub) / 2
            if len(fun_args) > 1:
                fun_val = fun(network, route_idx, x, total_flow, fun_idx)
            else:
                fun_val = fun(network, x)                
    else:
        while (abs(fun_val) > threshold[0]) or ((ub - lb) > threshold[1]):
            if fun_val > 0:
                lb = x
            else:
                ub = x
            x = (lb + ub) / 2
            if len(fun_args) > 1:
                fun_val = fun(network, route_idx, x, total_flow, fun_idx)
            else:
                fun_val = fun(network, x)                    
    return x
    
'''      
    if fun_args[0] == 1:
        total_flow = fun_args[1]
        route = fun_args[2]
        fun_val = fun(x, total_flow, route)
        if mode:
            while (abs(fun_val) > threshold[0]) or ((ub - lb) > threshold[1]):
                if fun_val > 0:
                    ub = x
                else:
                    lb = x
                x = (lb + ub) / 2         
                fun_val = fun(x, total_flow, route)
        else:
            while (abs(fun_val) > threshold[0]) or ((ub - lb) > threshold[1]):
                if fun_val > 0:
                    lb = x
                else:
                    ub = x
                x = (lb + ub) / 2
                fun_val = fun(x, total_flow, route)
        return x
    elif fun_args[0] == 2:
        network = fun_args[1]
        [fun_val, flows] = fun(x, network)
        if mode:
            while abs(fun_val) > threshold[0] or (ub - lb) > threshold[1]:
                if fun_val > 0:
                    ub = x
                else:
                    lb = x
                x = (lb + ub) / 2
                [fun_val, flows] = fun(x, network)
        else:
            while abs(fun_val) > threshold[0] or (ub - lb) > threshold[1]:
                if fun_val > 0:
                    lb = x
                else:
                    ub = x
                x = (lb + ub) / 2
                [fun_val, flows] = fun(x, network)
        return [x, flows]
    elif fun_args[0] == 3:

'''
        


        
    
