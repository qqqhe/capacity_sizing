# -*- coding: utf-8 -*-
"""
Created on Mon May 21 17:18:21 2018

Author: Qie He
"""

import time

import numpy as np

from search_methods import bisection_search

from settings import Settings

import inner_optimization as ip

def diff_lhs_rhs_equilibrium(flow, totalflow, route):
    """An auxiliary function $f_z(q) = \log{q_j} - \log{F_j(q)}$ computing the difference between the lhs and 
    rhs of the equilibrium equation $\log{q_j} = \log{F_j(q)}$ over a route, with the total flow over the network being z. 
    If $f_z(q^*)=0$, then $q^*$ is the equilibrium flow."""        
    difference = np.log(flow) + route.beta*(flow**route.congestion_power) - route.b +\
                route.parking_attraction*flow/route.parking_capacity - np.log(1-totalflow)    
    return difference
    
def diff_flowsum_totalflow(totalflow, network):
    """An auxiliary function Omega(z) computing the difference between the sum of the
    flows $\sum_{j} q_j(z)$ and the given total flow z.
    If Omega(z*)=0, then z* is the sum of the flows at the equilibrium of the network."""
    flows = np.zeros(network.num_routes)
    for j in range(network.num_routes):
        # Compute the flow over the route given the total flow over the network
        # --version 1---
        #network.routes[j].compute_flow(totalflow)
        #flows[j] = network.routes[j].flow    
        
        # --version 2---
        flows[j] =  bisection_search(0, 1, [1e-6, 1e-6], True, diff_lhs_rhs_equilibrium, 1, totalflow, network.routes[j])            
    difference = np.sum(flows) - totalflow
    return [difference, flows]
        
# Needs to rewrite
def zeta(x, z, nw_setting, fun_idx, route_idx):
    """An auxiliary function $\zeta(x, z)$ used to compute the bound of the flow over a route
    given by route_idx."""
    fun_val = np.log(x) + nw_setting.beta[route_idx] * x**nw_setting.congestion_power[route_idx] -\
                nw_setting.b[route_idx] - np.log(1-z) 
    if fun_idx == 1:
        fun_val += nw_setting.parking_attraction
    elif fun_idx == 2:
        fun_val += nw_setting.parking_attraction*x/nw_setting.lb_capacity[route_idx]
    elif fun_idx == 3:
        fun_val += nw_setting.parking_attraction*x/nw_setting.ub_capacity[route_idx]
    return fun_val 


class Route():
    """Define a class object for a single route"""
    
    def __init__(self, nw_setting, idx=0):
        """Initialize the parameters of the idx-th route over the network"""        
        # Route-specific parameters
        self.congestion_scalar = nw_setting.congestion_scalar[idx] # the scalar in the time-flow function, alpha_prime
        self.congestion_power = nw_setting.congestion_power[idx] # the power in the time-flow function, theta
        self.time_flowfree = nw_setting.time_flowfree[idx] # flow-free travel time, a 
        self.parking_capacity = 100 # the default scaled parking lot capacity, set to a large number
        
        # User-specific parameters for a route
        self.intrinsic_utility = nw_setting.intrinsic_utility[idx] # the intrinsic utility of a route, c
        
        # User-specific parameters for a network
        self.time_sensitivity = nw_setting.time_sensitivity # the unit disutility of travel time, alpha
        self.parking_attraction = nw_setting.parking_attraction # the coefficient for available parking lots in the utility function, phi

        # Auxiliary parameters in the user utility function
        self.beta = nw_setting.beta[idx]
        self.b =  nw_setting.b[idx]

        # Attributes that need to be computed over a network        
        self.flow = 0
        
    def get_travel_time(self):
        """Compute the trave time over the route with the given flow"""
        travel_time = self.time_flowfree +\
            self.congestion_scalar * self.flow ** self.congestion_power
        return travel_time        
        
    def get_route_utility(self):
        """Compute the utility of the route with the given flow"""
        route_utility = self.b - self.beta*(self.flow**self.congestion_power) -\
             self.parking_attraction*self.flow/self.parking_capacity
        return route_utility 
            
    def compute_flow(self, totalflow):
        """Compute the flow q over the route that satisfies the equation $q_j = F_j(q)$,
        with the given total flow over the network."""
        if totalflow < 1 - np.exp(self.beta - self.b + self.parking_attraction/self.parking_capacity):
            print("The value of the total flow is too small. Please select a larger value.")
        else:
            self.flow = bisection_search(0, 1, [1e-6, 1e-6], True, diff_lhs_rhs_equilibrium, 1, totalflow, self)
    
    def update_capacity(self, capacity):
        """Update the parking capacity over the route"""
        self.parking_capacity = capacity
    
    def update_flow(self, flow):
        """Update the equiplibrium flow over the route"""
        self.flow = flow
    
    def compute_capacity(self, totalflow):
        """Given the flow over the route and the total flow over the network at equilibrium, compute the corresponding parking capacity"""
        self.parking_capacity = self.parking_attraction * self.flow / (self.b \
        - self.beta*(self.flow**self.congestion_power) + np.log(1-totalflow) - np.log(self.flow))
        

class Network():
    """Define a class object for a road network with several parallel routes between a pair of origin and destination"""
    
    def __init__(self, nw_setting, homo_routes=False):
        """Initialize the network"""
        self.num_routes = nw_setting.num_routes
        self.routes = []        
        self.social_welfare = 0

        # User-specific parameters for a network
        self.time_sensitivity = nw_setting.time_sensitivity # the unit disutility of travel time, alpha
        self.parking_attraction = nw_setting.parking_attraction # the coefficient for available parking lots in the utility function, phi

        self.create_routes(homo_routes)
        self.compute_social_welfare()        
        
    def create_routes(self, homo_routes):
        """Create a set of routes for the network"""
        if homo_routes:
            route = Route(nw_setting)
            for i in range(self.num_routes):
                self.routes.append(route)
        else:
            for i in range(self.num_routes):            
                route = Route(nw_setting, i)
                self.routes.append(route)            

    def compute_social_welfare(self):
        """Compute the social welfare over the network"""
        self.social_welfare = 0
        for i in range(self.num_routes):
            self.social_welfare += self.routes[i].flow * self.routes[i].get_route_utility()
   
    def compute_lb_totalflow(self):
        """Compute the appropriate lower bound of the total flow over the network"""
        lb_totalflow = 0
        for i in range(self.num_routes):
            temp = 1 - np.exp(self.routes[i].beta - self.routes[i].b +\
                              self.routes[i].parking_attraction/self.routes[i].parking_capacity)
            if lb_totalflow < temp:
                lb_totalflow = temp
        return lb_totalflow
        
    def compute_equilibrium(self):
        """Compute the flow equilibrium over the network"""
        # Set the appropriate lower bound for the total flow over a network
        lb_totalflow = self.compute_lb_totalflow()
        # Use binary search to find the total flow at equilibrium
        [totalflow_at_equilibrium, flow_at_equilibrium] =\
            bisection_search(lb_totalflow, 1, [1e-6, 1e-6], False, diff_flowsum_totalflow, 2, self)
        
        for i in range(self.num_routes):
            self.routes[i].flow = flow_at_equilibrium[i]



def search_optimal_capacities(nw_setting, step_size, epsilon, network):
    """Find parking lot capacities in a network to maximize the social welfare at equilibrium"""
    
    # Compute the lower bound for the total flow over a network
    lb_totalflow = np.amax(1 - np.exp(nw_setting.beta - nw_setting.b +\
           nw_setting.parking_attraction*np.minimum(1, 1/nw_setting.lb_capacity)))
    lb_totalflow = max(lb_totalflow, epsilon)    
    
    # Initialize the value of total flow
    totalflow = lb_totalflow
    
    # An auxiliary threshold of the total flow computed based on the capacity upper bounds
    threshold_ub_capacity = 1 - np.exp(nw_setting.beta - nw_setting.b \
                            + nw_setting.parking_attraction/nw_setting.ub_capacity)
    #print(threshold_ub_capacity)
    
    # Initialize the bounds for flow over each route
    ub_flow = np.zeros(nw_setting.num_routes)
    lb_flow = np.zeros(nw_setting.num_routes)
    
    # Initialize the optimal social welfare over the network
    opt_socialwelfare = - nw_setting.num_routes/np.e 
    opt_totalflow = None
    opt_flows = np.array([])
    opt_capacity = np.zeros(nw_setting.num_routes)

    # For debugging
    lower_bound = np.zeros(nw_setting.num_routes)
    upper_bound = np.array(nw_setting.num_routes)
    count = 0
        
    while totalflow < 1 - epsilon:
        flag_nofeasibleflow = False
        for i in range(nw_setting.num_routes):
            if totalflow >= threshold_ub_capacity[i]:
            # Line 5-7 of Algorithm 3
                x3_star = bisection_search(0, 1, [epsilon, epsilon], True, zeta, 3, nw_setting, totalflow, 3, i)
                if x3_star > nw_setting.ub_capacity[i]:
                    flag_nofeasibleflow = True
                    break                    
                else:
                    ub_flow[i] = x3_star 
            else:                               
                ub_flow[i] = 1
            
            # Computing x1_star and x2_star 
            x1_star = bisection_search(0, 1, [epsilon, epsilon], True, zeta, 3, nw_setting, totalflow, 1, i)
            x2_star = bisection_search(0, 1, [epsilon, epsilon], True, zeta, 3, nw_setting, totalflow, 2, i)
            lb_flow[i] = max(x1_star, x2_star)
        
        if not flag_nofeasibleflow:
            # The implementation of line 11 to 18
            [opt_obj, opt_x] = ip.max_sum_xlogx(nw_setting.num_routes, totalflow, lb_flow, ub_flow)        
            
            # Compute the value of $\pi(z)$ and update the optimal z if possible
            if opt_x.size > 0:
                temp = opt_obj - totalflow * np.log(1-totalflow)
                if temp > opt_socialwelfare:
                    opt_socialwelfare = temp
                    opt_totalflow = totalflow
                    print("Update optimal flow")
                    print(opt_x)
                    print(lb_flow)
                    print(ub_flow)
                    opt_flows = opt_x
                    
                    # For debugging
                    lower_bound = lb_flow
                    upper_bound = ub_flow 
                    count += 1
                
        totalflow += step_size
    
    # Line 20 of ALgorithm 3
    if opt_flows.size > 0:          
        for i in range(nw_setting.num_routes):
            network.routes[i].update_flow(opt_flows[i])
            network.routes[i].compute_capacity(opt_totalflow)
            opt_capacity[i] = network.routes[i].parking_capacity
#            print("The optimal flow over the " + str(i+1) + "-th route is " + str(opt_flows[i]))
#            print("The optimal capacity over the " + str(i+1) + "-th route is " + str(network.routes[i].parking_capacity))
        print("The optimal flow is: ")
        print(opt_flows)
        print("The optimal parking capacity is: ")
        print(opt_capacity)            
        print("--The optimal total flow is " + str(opt_totalflow))
        print("The maximum social welfare is " + str(opt_socialwelfare) +".")  
        
    else:
        print("No optimal solution is found!")
        
    
    # For debugging
    temp1 = np.zeros(nw_setting.num_routes)
    temp2 = np.zeros(nw_setting.num_routes)
    temp3 = np.zeros(nw_setting.num_routes)
    for i in range(nw_setting.num_routes):
        temp1[i] = zeta(opt_flows[i], opt_totalflow, nw_setting, 1, i)
        temp2[i] = zeta(opt_flows[i], opt_totalflow, nw_setting, 2, i)
        temp3[i] = zeta(opt_flows[i], opt_totalflow, nw_setting, 3, i)
    print("The function value of zeta at the optimal flow: ")
    print(temp1)
    print(temp2)
    print(temp3)
        
   # For debugging
    print("The capacity upper bound when optimal flow is found: ")
    print(lower_bound)
    print("The capacity lower bound when optimal flow is found: ")
    print(upper_bound)
    print(str(count))
#--------------------Testing----------------------------------------------------

"""       
route = Route()
print("The travel time over the route is " + str(route.get_travel_time()))
print("The commuter's utility over the route is " + str(route.get_route_utility())) 
"""   
       
                

np.random.seed(10)    
nw_setting = Settings()
network = Network(nw_setting)
#for i in range(nw_setting.num_routes):
#    print("The travel time over the route is " + str(network.routes[i].get_travel_time()))
#    print("The commuter's utility over the route is " + str(network.routes[i].get_route_utility())) 
#    #print("The equilibrium flow over the route is " + str(network.routes[i].compute_flow(0.5)))
#    network.routes[i].compute_flow(0.5)
#    print("The equilibrium flow over the route is " + str(network.routes[i].flow))
#    
#print("The social welfare over the network is " + str(network.social_welfare))
#
#[diff, flows] = diff_flowsum_totalflow(0.8, network)
#for i in flows:
#    print("The flow is " + str(i) + ".")
#    
#network.compute_equilibrium()
#
#for i in range(nw_setting.num_routes):
#    print("The travel time over the route is " + str(network.routes[i].get_travel_time()))
#    print("The commuter's utility over the route is " + str(network.routes[i].get_route_utility())) 
#    #print("The equilibrium flow over the route is " + str(network.routes[i].compute_flow(0.5)))    
#    print("The equilibrium flow over the route is " + str(network.routes[i].flow))
#
#network.compute_social_welfare()
#print("The social welfare over the network is " + str(network.social_welfare))

print("The upper bound on the capacity is: ")
print(nw_setting.ub_capacity)
print("The lower bound on the capacity is: ")
print(nw_setting.lb_capacity)


print("----The result of optimal capacity sizing-------")
start_time = time.time()
search_optimal_capacities(nw_setting, 1e-4, 1e-6, network)
print("--- %s seconds ---" % (time.time() - start_time))
#network.compute_equilibrium()
#flow = np.zeros(nw_setting.num_routes)
#for i in range(nw_setting.num_routes):
#    flow[i] = network.routes[i].flow
#print("The flow recalculated: ")
#print(flow)