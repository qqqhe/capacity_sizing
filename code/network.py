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

import matplotlib.pyplot as plt



        
     
def diff_totalflow_sumflow(network, z):
    """
    An auxiliary function $\Omega(z)$ computing the difference between the given total flow z over a network and the sum of the flow $\sum_{j} q_j(z)$, where $q_j(z)$ is computed from the equilibrium function for the j-th route. If $\Omega(z^*)=0$, then $z^*$ is the sum of the flows at the equilibrium of the network.
    The function $\Omega(z)$ is a monotone non-increasing function, $\Omega(1)=1$ and $\Omega(0)<0$.
    Version: Modified on Jan 28, 2019
    
    Parameters
    ----------
    network: class object
        Store all info of the given network.
    z: float
        The value of total flow over the network            
    
    Returns
    -------
    diff: float or []
        The difference between $z$ and the sum of $q_j(z)$
    """
    beta = network.beta
    b = network.b
    C = network.capacity
    phi = network.phi
    theta = network.theta
    
    ## Testing
    print("The current value of z is: " + str(z) + "\n")
    ##
    
    # First make sure with the given z value, we can compute the equilibrium q over each route
    if z >= np.amax(1 - np.multiply(C, np.exp(beta * np.power(C, theta) + phi - b))) and z < 1:
        # Compute the equilibrium $q_j(z)$
        q = np.zeros(network.num_routes)
        for j in range(network.num_routes):
            q[j] = bisection_search(zeta, 0, C[j], [1e-10, 1e-10], True, network, z, j, 4)
        
        diff = z - np.sum(q)
        
        ## Testing
        print("The difference is: " + str(diff) + "\n")
        ##
        return diff
    else:
        print("The given z value is not a valid total flow.")



def compute_social_welfare(flow):
    """
    Compute the social welfare over the flow over the network under the MNL model
    
    Parameters
    ----------
    flow: an array of floats
        The i-th component is the equilibrium flow over the i-th route        
    
    Returns
    -------
    social_welfare: float
    """
    total_flow = np.sum(flow)
    if 0 < total_flow < 1:
        social_welfare = np.sum(np.multiply(flow, np.log(flow))) - total_flow * np.log(1 - total_flow)
        return social_welfare                
    else:
        print("The given flow vector is not feasible.")
        return 0
        
def diff_equilibrium_equation(flow, capacity, network, route_idx):
    """
    An auxiliary function computing the difference between the lhs and rhs of the equilibrium equation $\log{q_j} = \log{F_j(q)}$ for the i-th route of a network.
    
    Parameters
    ----------
    flow: an array of floats
        A flow vector.
    capacity: an array of floats
        A vector storing parking capacity over all routes.
    network: class Network object
    route_idx: int
        The index of the route.
    
    Returns       
    -------
    diff: float
        The difference between the lhs and the rhs.    
    """  
    total_flow = np.sum(flow)
    if 0 < total_flow < 1:
        diff = np.log(flow[route_idx]) + network.beta * (flow[route_idx]**network.theta) - network.b[route_idx] + \
                network.phi * flow[route_idx]/capacity[route_idx] - np.log(1-total_flow)  
        return diff                
    else:
        print("The given flow vector is not feasible.")
        return 0
        
def check_feasibility(network, flow, capacity):
    """
    Check if a solution (flow, capacity) is feasible for the given instance; three sets of constraints are checked: (1) equilibrium equations; (2) bound constraints for capacities; (3) bound constraints for the flow.    
    
    Parameters
    ----------
    network: class Network object
        An object contains all information about a capacity allocation instance.
    flow: an array of floats
        Flow over the network
    capacity: an array of floats
        Capacity over the network
    
    Returns
        Boolean.
        True if feasible, False otherwise.
    -------
    """
    diff_eq = np.zeros(network.num_routes)
    for i in range(network.num_routes):
        diff_eq[i] = abs(diff_equilibrium_equation(flow, capacity, network, i))        
        
    diff_cap = max(np.amax(capacity - network.u), np.amax(network.l - capacity))
    diff_flow = np.amax(flow - capacity)        
    
    if np.amax(diff_eq) < 1e-5 and diff_cap < 1e-5 and diff_flow < 1e-5:
        return True
    else:
        return False
    
 
    

'''
class Route():
    """Define a class object for a single route"""
    
    def __init__(self, network, capacity_init_value=100, idx=0):
        """Initialize the parameters of the idx-th route over the network"""        
        #---Route-specific parameters---
        self.congestion_multiplier = network.congestion_multiplier[idx] 
        self.congestion_power = network.congestion_power[idx]
        self.time_flowfree = network.time_flowfree[idx]  
        
        # Attributes that need to be computed over a network        
        self.flow = 0        

        # The default scaled parking lot capacity, set to a large number
        self.parking_capacity = capacity_init_value
        
        # User-specific parameters
        # Parameter specific for a route
        self.intrinsic_utility = network.intrinsic_utility[idx] 
        # Parameters for a network
        self.time_sensitivity = network.time_sensitivity 
        self.parking_attraction = network.parking_attraction 

        # Auxiliary parameters in the user utility function
        self.beta = network.beta[idx]
        self.b =  network.b[idx]       
        
    def get_travel_time(self):
        """Compute the trave time over the route with the given flow"""
        travel_time = self.time_flowfree +\
            self.congestion_multiplier * self.flow ** self.congestion_power
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
'''
        

class Network():
    """Define a class object for a road network with several parallel routes 
    between a pair of origin and destination"""
    
    def __init__(self, num_routes, beta, phi, theta, b, l, u):
        """Initialize the network parameters"""        
        self.num_routes = num_routes
        self.beta = beta
        self.phi = phi
        self.theta = theta
        self.b = b
        self.l = l
        self.u = u

        self.lb_totalflow = self.compute_lb_totalflow()
        # Initialize the flow over each route to be 0
        self.flow = np.zeros(num_routes) 
        # Initialize the capacity over a route to be 0
        self.capacity = np.zeros(num_routes)                

    def compute_social_welfare(self):
        """Compute the social welfare over the network"""
        #self.social_welfare = 0
        #for i in range(self.num_routes):
        #    self.social_welfare += self.routes[i].flow * self.routes[i].get_route_utility()
   
    def compute_lb_totalflow(self):
        """Compute the appropriate lower bound of the total flow over the network"""
        lb_totalflow = np.amax(1 - np.exp(self.beta - self.b + self.phi * np.minimum(1, 1/self.l)))        
        return max(0, lb_totalflow)
        
    def update_flow(self, flow):
        """Update the equiplibrium flow over the network"""
        self.flow = flow
        
    def compute_capacity(self, totalflow, i):
        """Given the flow over the i-th route and the total flow over
        the network at equilibrium, compute the corresponding parking capacity"""
        self.capacity[i] = self.phi * self.flow[i] / (self.b[i] \
        - self.beta*(self.flow[i]**self.theta) + np.log(1-totalflow) - np.log(self.flow[i]))
        
    
    def compute_equilibrium(self):
        """
        Compute the flow equilibrium over the network with the given capacities. 
        Version: Jan 28, 2019.
        """        
        # First compute a valid lower bound for the total flow
        totalflow_lb = max(0, np.amax(1 - np.multiply(self.capacity, np.exp(self.beta * np.power(self.capacity, self.theta) + self.phi - self.b))))
        if totalflow_lb > 1:
            print("This network does not have equilibrium!")
        else:
            # Compute the total flow at equilibrium z_star
            z_star = bisection_search(diff_totalflow_sumflow, totalflow_lb, 1, [1e-10, 1e-10], True, self)   
            # Compute the flow over each route at equilibrium 
            for i in range(self.num_routes):
                self.flow[i] = bisection_search(zeta, 0, self.capacity[i], [1e-10, 1e-10], True, self, z_star, i, 4)
    


def zeta_homo_noparking(network, x):
    """
    An auxiliary function to help evaluate the difference between the lhs and the rhs of the flow equilibrium function of the homogeneous instance when parking information is not given and the flow over each route is x.
    This functioin is a non-decreasing function of x.
    
    Parameters
    ----------
    network: contains all parameters related to the network
       
    Returns
    -------
    fun_val: float
        The value of the difference between the lhs and rhs of the flow equilibrium function.        
    """        
    fun_val = network.beta * x**network.theta + np.log(x) - np.log(1 - network.num_routes*x) - network.b[0]
    return fun_val 


def zeta(network, i, x, z, fun_idx):
    """
    An auxiliary function used to compute the bound of the flow over the i-th route.
    
    Parameters
    ----------
    network: contains all parameters related to the network
    i: the index of the route
    x: the flow over the i-th route
    z: the total flow over the network
    fun_idx: indicate which function to use to compute the bound
    
       
    Returns
    -------
    fun_val: float or []
        The value of the zeta function        
    """
        
    fun_val = network.beta * x**network.theta + np.log(x) - np.log(1-z) - network.b[i]

    if fun_idx == 1:
        fun_val += network.phi
    elif fun_idx == 2:
        fun_val += network.phi * x / network.l[i]
    elif fun_idx == 3:
        fun_val += network.phi * x / network.u[i]
    elif fun_idx == 4:
        fun_val += network.phi * x / network.capacity[i]
    else:
        print("Wrong index for the zeta function")
        fun_val = []
    return fun_val 


def search_optimal_capacities(network, step_size, tolerance, filename):
    """
    Find optimal parking-lot capacities in a network to maximize the social welfare at equilibrium
    
    Parameters
    ----------
    network: Network object
        Store all information needed for an instance.
    step_size: float
        The incremental step size used to search the optimal total flow over the network.
    tolerance: float
        The tolerance for terminating the bisection search
    filename: string
        Auxiliary argument to store the plot h(z) into a file 
        
    
    Returns
    opt_flow: array of floats
        Store the optimal flow over each route
    opt_cap: array of floats
        The optimal capacity over each route
    opt_socialwelfare: float
        The optimal social-welfare    
    ----------    
    """
    ## Initialization
    # Initialize the value of total flow over the network
    totalflow = max(network.lb_totalflow, step_size)
    
    # An auxiliary threshold of the total flow computed based on the capacity upper bounds, used in Line 4 of Algorithm 3.
    aux_bound = 1 - np.exp(network.beta - network.b + network.phi/network.u)
    
       
    # Initialize the bounds for flow over each route
    ub_flow = np.zeros(network.num_routes)
    lb_flow = np.zeros(network.num_routes)
    
    # Initialize the optimal solution over the network
    opt_socialwelfare = np.array([])
    opt_totalflow = 0
    opt_flows = np.array([])
    opt_capacity = np.zeros(network.num_routes)
    

#    # For debugging only
#    lower_bound = np.zeros(network.num_routes)
#    upper_bound = np.zeros(network.num_routes)
#    count = 0
    
    # Try to plot out the (totalflow, social_welfare) scatter plot
    z = []
    hz = []
#    # End of debugging

    ## Start the search
    while totalflow < 1 - tolerance:
        flag_nofeasibleflow = False
       
        # Compute the bounds for the flow.
        for i in range(network.num_routes):
             # Line 3-8 of Algorithm 3. Compute the upper bounds for the flow.
            if totalflow >= aux_bound[i]:            
                x3_star = bisection_search(zeta, 0, 1, [tolerance, tolerance], True, network, totalflow, i, 3)              
                if x3_star > network.u[i]:
                    flag_nofeasibleflow = True
                    break                    
                else:
                    ub_flow[i] = x3_star 
            else:                               
                ub_flow[i] = 1            
            # Line 9-10 of Algorithm 3. Compute the lower bounds of the flow.
            x1_star = bisection_search(zeta, 0, 1, [tolerance, tolerance], True, network, totalflow, i, 1)
            x2_star = bisection_search(zeta, 0, 1, [tolerance, tolerance], True, network, totalflow, i, 2)
            lb_flow[i] = max(x1_star, x2_star)
        
            
        if not flag_nofeasibleflow:
            # Check feasibility of the flow based on the current total flow, lower and upper bounds of the flow
            if totalflow < np.sum(lb_flow) or totalflow > np.sum(ub_flow):                
                totalflow += step_size                

#                # For debugging only
#                print("\nThe current total flow is: " + str(totalflow))
#                print("\nThe capacity upper bound when optimal flow is found: ")
#                print(upper_bound)
#                print("\nThe capacity lower bound when optimal flow is found: ")
#                print(lower_bound)
#                print(str(count))
#                # Eng of debugging
#                
                continue
                        
            # The implementation of line 11 to 18. Find the optimal flow given the current value of z.
            [opt_obj, opt_x] = ip.max_sum_xlogx(network.num_routes, totalflow, lb_flow, ub_flow)        
                      
            
            # Line 18 of Algorithm 3. Compute the social welfare given the current z and optimal q(z).
            temp = opt_obj - totalflow * np.log(1-totalflow)

            ##### Testing: to plot out the function of h(z)
            z.append(totalflow)
            hz.append(temp)
            ##### End of Testing: to plot out the function of h(z)
               
            if opt_socialwelfare.size == 0 or temp > opt_socialwelfare:
                opt_socialwelfare = temp
                opt_flows = opt_x
                opt_totalflow = totalflow                
                
                # For debugging only
#                print("\nUpdate optimal flow")
#                print(opt_x)
#                print(lb_flow)
#                print(ub_flow)
#                print("Total flow is " + str(opt_totalflow))                
                
                # For debugging
#                np.copyto(lower_bound, lb_flow)                
#                np.copyto(upper_bound, ub_flow) 
#                count += 1
#                    print("The lower and upper bounds are: ")
#                    print(lb_flow)
#                    print(lower_bound)
#                    print("\n")
#                    print(ub_flow)
#                    print(upper_bound)
#                    print("\n")
                
        totalflow += step_size    

    
        
#    # For debugging only
#    print("\n----------------\n Exiting the while loop.")
#    print("\nThe capacity upper bound when optimal flow is found: ")
#    print(upper_bound)
#    print("\nThe capacity lower bound when optimal flow is found: ")
#    print(lower_bound)
#    print(str(count))   
#    # Eng of debugging
    
    # Line 20 of ALgorithm 3
    if opt_flows.size > 0:
        network.update_flow(opt_flows)          
        for i in range(network.num_routes):            
            network.compute_capacity(opt_totalflow, i)
            opt_capacity[i] = network.capacity[i]
        print("\n--------------\nThe optimal flow is: ")
        print(opt_flows)
        print("\n--------------\nThe optimal parking capacity is: ")
        print(opt_capacity)            
        print("\n--------------\nThe optimal total flow is " + str(opt_totalflow))
        print("\n--------------\nThe maximum social welfare is " + str(opt_socialwelfare) +".")
        
        
        ##### Testing: to plot out the function of h(z)
        #plt.scatter(z, hz, c='r', marker='r')
        plt.plot(z, hz, '-', linewidth=0.5)
        #plt.xlim(0.5, 1)
        plt.savefig(filename + '.png', bbox_inches='tight')
        ##### End of Testing: to plot out the function of h(z)
        
       
        
#        # For debugging
#        temp1 = np.zeros(network.num_routes)
#        temp2 = np.zeros(network.num_routes)
#        temp3 = np.zeros(network.num_routes)
#        for i in range(network.num_routes):        
#            temp1[i] = zeta(network, i, opt_flows[i], opt_totalflow, 1)
#            temp2[i] = zeta(network, i, opt_flows[i], opt_totalflow, 2)
#            temp3[i] = zeta(network, i, opt_flows[i], opt_totalflow, 3)
#        print("The function value of zeta at the optimal flow: ")
#        print(temp1)
#        print(temp2)
#        print(temp3)
#            
#       # For debugging
#        print("\nThe capacity upper bound when optimal flow is found: ")
#        print(upper_bound)
#        print("\nThe capacity lower bound when optimal flow is found: ")
#        print(lower_bound)
#        print(str(count))
#       # End of debugging
               
        return opt_flows, opt_capacity, opt_socialwelfare       
    else:
        print("\nNo optimal solution is found!")
        return np.array([]), opt_capacity, opt_socialwelfare
    
    
    


#--------------------Testing----------------------------------------------------

"""       
route = Route()
print("The travel time over the route is " + str(route.get_travel_time()))
print("The commuter's utility over the route is " + str(route.get_route_utility())) 
"""   
       
                
'''
np.random.seed(10)    
network = Settings()
network = Network(network)
#for i in range(network.num_routes):
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
#for i in range(network.num_routes):
#    print("The travel time over the route is " + str(network.routes[i].get_travel_time()))
#    print("The commuter's utility over the route is " + str(network.routes[i].get_route_utility())) 
#    #print("The equilibrium flow over the route is " + str(network.routes[i].compute_flow(0.5)))    
#    print("The equilibrium flow over the route is " + str(network.routes[i].flow))
#
#network.compute_social_welfare()
#print("The social welfare over the network is " + str(network.social_welfare))

print("The upper bound on the capacity is: ")
print(network.capacity_ub)
print("The lower bound on the capacity is: ")
print(network.capacity_lb)


print("----The result of optimal capacity sizing-------")
start_time = time.time()
search_optimal_capacities(network, 1e-4, 1e-6, network)
print("--- %s seconds ---" % (time.time() - start_time))
#network.compute_equilibrium()
#flow = np.zeros(network.num_routes)
#for i in range(network.num_routes):
#    flow[i] = network.routes[i].flow
#print("The flow recalculated: ")
#print(flow)
'''