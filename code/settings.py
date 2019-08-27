"""
A class to store all the route and user related parameters of a  road network

Created on Mon May 21 17:06:45 2018

Author: Qie He
"""

import numpy as np

class Settings():
    """A class to store all the parameters for the optimal capacity problem"""
    
    def __init__(self, num_routes, beta, phi, theta, b, l, u):
        """Initialize all the parameters"""
        self.num_routes = num_routes
        self.beta = beta
        self.phi = phi
        self.theta = theta
        self.b = b
        self.l = l
        self.u = u
        
        '''
        #---Route-specific parameters---        
        # Define the scalar in the time-flow function, $\alpha'$ in the paper,
        # and set its value to 1.
        self.congestion_multiplier = np.ones(num_routes) 
        
        # Define the power in the time-flow function, $\theta$ in the paper,
        # and set its value to 2.
        self.congestion_power = theta
        
        # Define flow-free travel time, $a$ in the paper
        self.time_flowfree = time_flowfree_lb + np.random.rand(num_routes)*(time_flowfree_ub - time_flowfree_lb)  
        
        #---User-specific parameters for a network---
        # Define the unit disutility of travel time, $\alpha$ in the paper
        self.time_sensitivity = time_sensitivity_lb + np.random.rand()*(time_sensitivity_ub - time_sensitivity_lb) 
       
        # Define a user's intrinsic utility of each route
        self.intrinsic_utility = intrinsic_utility_lb + np.random.rand(num_routes)*(intrinsic_utility_ub - intrinsic_utility_lb) # the intrinsic utility of a route, c
            
        #---Auxiliary parameters in the user utility function---
        self.beta = self.congestion_multiplier*self.time_sensitivity
        self.b =  self.intrinsic_utility - self.time_sensitivity*self.time_flowfree + self.parking_attraction
        #self.b = 2*np.random.rand(num_routes)
        '''
        
    def display_settings(self):
        '''Display the parameters in the setting'''
        print('The parameters of the newly generated road network.\n')
        print('The number of routes is {0}.'.format(str(self.num_routes)))
        print('\nThe free flow time over the network:')
        for i in range(self.num_routes):
            print('Route {0}: {1}'.format(str(i+1), str(self.time_flowfree[i])))
        
        print('\nThe intrinsic utility of a user:')
        for i in range(self.num_routes):
            print('Route {0}: {1}'.format(str(i+1), str(self.intrinsic_utility[i])))
    
        print('\nThe range of parking capacity to search over each route:')
        for i in range(self.num_routes):
            print('Route {0}: [{1}, {2}]'.format(str(i+1), str(self.capacity_lb[i]), str(self.capacity_ub[i])))
        
       
if __name__ == '__main__':
    s = Settings()
    s.display_settings()