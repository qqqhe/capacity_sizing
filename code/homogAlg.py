# This script file finds the unique equilibrium point defined in \bq = \bF(\bq)
# with any given \bC and under the MNL choice model, which is equivalent to
# \log(q_{j}^{0}) = \log\left(1-\|\bq^{0}\|_{1}\right)+b_{j}-\beta (q_{j}^{0})^{\theta}-\varphi q_{j}/C_{j}
# \ \ \forall \ j \in [J].
#                                            By: Xinchang Wang
#                                1st version: April 3, 2018, 2nd version: April 4, 2018
# import scipy.special as sp -- not used
import numpy as np


# Define product-specific function for which we need to find the unique root
# product_flow: number of commuters picking up the route
# product_idx: the index number of the route
# flow_norm: \|\bq\|_{1}, the sum of all the flows. Its domain is [(1-e^(\beta-b_{product_idx}))^{+},1]
# para_list: the input parameter list
def product_spec_fun(product_flow, product_idx, flow_norm, para_list):
    beta = para_list[2]
    theta = para_list[3]
    varphi = para_list[7]
    intrinsic_u = para_list[8]
    capacity = para_list[9]
    f = beta * product_flow ** theta - intrinsic_u[product_idx] + np.log(product_flow) \
        + varphi*product_flow/capacity[product_idx] - np.log(1 - flow_norm)
    return f


# Define grand function -- the big Omega function
# omegas: the roots for product-specific functions, it is a list
def inverse_fun_total(flow_norm, paras):
    J = paras[1]
    omegas = np.zeros(J)
    for product_idx in range(J):
        omegas[product_idx] = bisection_search(0, 1, True, product_spec_fun, True,
                                               product_idx, flow_norm, para_list=paras)
    f = np.sum(omegas) - flow_norm
    return [f, omegas]


# Bi-section search method to find a root of a monotone function, if it exists.
# lb: lower bound
# ub: upper bound
# mode = true: monotone increasing function; decreasing otherwise
# fun: the function for which you need to find a root
# '''
# *args: additional variables defining the function:
# (1) args = [True, product_idx, flow_norm] ==> fun takes "product_spec_fun"
# (2) args = [False] ==> fun takes "inverse_fun_total"
# '''
def bisection_search(lb, ub, mode, fun, *args, para_list):
    epsilon = [para_list[5], para_list[6]]
    if args[0]:
        product_idx = args[1]
        flow_norm = args[2]
        x = (lb + ub) / 2
        fun_val = fun(x, product_idx, flow_norm, para_list)
        if mode:
            while (abs(fun_val) > epsilon[0]) or ((ub - lb) > epsilon[1]):
                if fun_val > 0:
                    ub = x
                else:
                    lb = x
                x = (lb + ub) / 2
                fun_val = fun(x, product_idx, flow_norm, para_list)
        else:
            while abs(fun_val) > epsilon[0] or (ub - lb) > epsilon[1]:
                if fun_val > 0:
                    lb = x
                else:
                    ub = x
                x = (lb + ub) / 2
                fun_val = fun(x, product_idx, flow_norm, para_list)
        return x
    else:
        x = (lb + ub) / 2
        [fun_val, omegas] = fun(x, para_list)
        if mode:
            while abs(fun_val) > epsilon[0] or (ub - lb) > epsilon[1]:
                if fun_val > 0:
                    ub = x
                else:
                    lb = x
                x = (lb + ub) / 2
                [fun_val, omegas] = fun(x, para_list)
        else:
            while abs(fun_val) > epsilon[0] or (ub - lb) > epsilon[1]:
                if fun_val > 0:
                    lb = x
                else:
                    ub = x
                x = (lb + ub) / 2
                [fun_val, omegas] = fun(x, para_list)
        return [x, omegas]


# Class: equilibrium flows for input and output
class EquilibriumFlows(object):
    def __init__(self, seed, J, beta, theta, scale, epsilon_1, epsilon_2, varphi, intrinsic_u, capacity):
        self.seed = seed
        self.J = J
        self.beta = beta
        self.theta = theta
        self.scale = scale
        self.epsilon_1 = epsilon_1
        self.epsilon_2 = epsilon_2
        self.intrinsic_u = intrinsic_u
        self.varphi = varphi
        self.capacity = capacity

    def print_root(self):
        flow_norm_lb = max(np.maximum(1 - np.exp(self.beta - self.intrinsic_u), 0)) + self.epsilon_1
        paras = [self.seed, self.J, self.beta, self.theta, self.scale, self.epsilon_1, self.epsilon_2,
                 self.varphi, self.intrinsic_u, self.capacity]
        [flow_norm_star, flows] = bisection_search(flow_norm_lb, 1, False, inverse_fun_total, False, para_list=paras)
        all_flows = np.transpose(np.asmatrix(flows))
        return [flow_norm_star, all_flows]


# input para_list: seed, J, beta, theta, scale, epsilon_1, epsilon_2, varphi, intrinsic_u, capacity
# seed = 10
# J = 10
# beta = 0.1
# theta = 0.5
# scale = 1
# epsilon_1 = 10 ** (-6)
# epsilon_2= 10 ** (-6)
# varphi = 0.2
input_paras = [100, 10, 0.1, 0.5, 1, 10 ** (-6), 10 ** (-6), 0.2]
# np.random.seed()
np.random.seed(input_paras[0])
# Inverse transform sampling: b_j's follow UNIFORM[-1,1]
intrinsic = input_paras[4] * (2 * np.random.rand(input_paras[1]) - 1)
# Generate the lower bound for \|\bq\|_{1}
flow_norm_lb = max(np.maximum(1 - np.exp(input_paras[2] - intrinsic), 0)) + input_paras[5]
# capacities are given and randomly generated following Uniform[0,1]
capacities = np.random.rand(input_paras[1])
result = EquilibriumFlows(input_paras[0], input_paras[1], input_paras[2], input_paras[3], input_paras[4],
                          input_paras[5], input_paras[6], input_paras[7], intrinsic, capacities)
# Output:
# flow_norm_1: \|\bq\|_{1}
# vector_flows: \bq
[flow_norm_1, vector_flows] = result.print_root()
print("\|bq\| is: " + str(flow_norm_1))
print("The vector of flows at equilibrium is: \n" + str(vector_flows))
print("The sum of flows is equal to: \n" + str(vector_flows.sum()))