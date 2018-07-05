# -*- coding: utf-8 -*-
"""
Created on Tue May 22 15:42:06 2018

Author: Qie He
"""


# Bi-section search method to find a root of a monotone function, if it exists.
# lb: lower bound
# ub: upper bound
# epsilon: tolerance for stopping the search
# mode = true: monotone increasing function; decreasing otherwise
# fun: the function for which you need to find a root
# '''
# *args: additional variables defining the function:
# (1) args = [1, total_flow, route]
# (2) args = [2, network]
# (3) args = [3, total_flow, network, function_index]

def bisection_search(lb, ub, epsilon, mode, fun, *args): 
    x = (lb + ub) / 2  
    if args[0] == 1:
        total_flow = args[1]
        route = args[2]
        fun_val = fun(x, total_flow, route)
        if mode:
            while (abs(fun_val) > epsilon[0]) or ((ub - lb) > epsilon[1]):
                if fun_val > 0:
                    ub = x
                else:
                    lb = x
                x = (lb + ub) / 2         
                fun_val = fun(x, total_flow, route)
        else:
            while (abs(fun_val) > epsilon[0]) or ((ub - lb) > epsilon[1]):
                if fun_val > 0:
                    lb = x
                else:
                    ub = x
                x = (lb + ub) / 2
                fun_val = fun(x, total_flow, route)
        return x
    elif args[0] == 2:
        network = args[1]
        [fun_val, flows] = fun(x, network)
        if mode:
            while abs(fun_val) > epsilon[0] or (ub - lb) > epsilon[1]:
                if fun_val > 0:
                    ub = x
                else:
                    lb = x
                x = (lb + ub) / 2
                [fun_val, flows] = fun(x, network)
        else:
            while abs(fun_val) > epsilon[0] or (ub - lb) > epsilon[1]:
                if fun_val > 0:
                    lb = x
                else:
                    ub = x
                x = (lb + ub) / 2
                [fun_val, flows] = fun(x, network)
        return [x, flows]
    elif args[0] == 3:
        network = args[1]
        total_flow = args[2]
        fun_idx = args[3]
        route_idx = args[4]
        fun_val = fun(x, total_flow, network, fun_idx, route_idx)
        if mode:
            while (abs(fun_val) > epsilon[0]) or ((ub - lb) > epsilon[1]):
                if fun_val > 0:
                    ub = x
                else:
                    lb = x
                x = (lb + ub) / 2         
                fun_val = fun(x, total_flow, network, fun_idx, route_idx)
        else:
            while (abs(fun_val) > epsilon[0]) or ((ub - lb) > epsilon[1]):
                if fun_val > 0:
                    lb = x
                else:
                    ub = x
                x = (lb + ub) / 2
                fun_val = fun(x, total_flow, network, fun_idx, route_idx)
        return x

        


        
    
