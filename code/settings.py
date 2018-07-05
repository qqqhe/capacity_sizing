# -*- coding: utf-8 -*-
"""
Created on Mon May 21 17:06:45 2018

Author: Qie He
"""

import numpy as np

class Settings():
    """A class to store all the settings for the optimal capacity problem"""
    
    def __init__(self, num_routes=5, a_lb=1, a_ub=10, c_lb=1, c_ub=10, alpha_lb=1, alpha_ub=2, phi_lb=0, phi_ub=2, C_lb=0.5, C_ub=2):
        """Initialize all the parameters"""
        #np.random.seed(10)
        
        # Parameters for a network
        self.num_routes = num_routes
        
        # User-specific parameters for a network
        self.time_sensitivity = alpha_lb + np.random.rand()*(alpha_ub - alpha_lb) # the unit disutility of travel time, alpha
        self.parking_attraction = phi_lb + np.random.rand()*(phi_ub - phi_lb) # the coefficient for available parking lots in the utility function, phi
        
        # User-specific parameters for a route
        self.intrinsic_utility = c_lb + np.random.rand(num_routes)*(c_ub - c_lb) # the intrinsic utility of a route, c
        
        # Route-specific parameters
        self.congestion_scalar = np.ones(num_routes) # the scalar in the time-flow function, alpha_prime
        self.congestion_power = 2*np.ones(num_routes) # the power in the time-flow function, theta
        self.time_flowfree = a_lb + np.random.rand(num_routes)*(a_ub - a_lb)  # flow-free travel time, a 
        self.lb_capacity = C_lb + np.random.rand(num_routes)*(C_ub - C_lb)
        self.ub_capacity = self.lb_capacity +\
            np.random.rand(num_routes)*(C_ub - self.lb_capacity)
        
        # Auxiliary parameters in the user utility function
        self.beta = self.congestion_scalar*self.time_sensitivity
        #self.b =  self.intrinsic_utility - self.time_sensitivity*self.time_flowfree + self.parking_attraction
        self.b = 2*np.random.rand(num_routes)
        
#setting = Settings()