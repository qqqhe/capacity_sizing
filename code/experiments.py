# -*- coding: utf-8 -*-
"""
The main code for running the optimal capacity design problem

Created on Wed Oct 24 13:56:02 2018

Author: Qie He
"""
import os
import numpy as np
from network import Network
from network import search_optimal_capacities
import matplotlib.pyplot as plt


def clear():
    '''
    Clear the screen
    '''
    if os.name in ('nt','dos'):
        os.system('cls')
    elif os.name in ('linux','osx','posix'):
        os.system('clear')
    else:
        print("\n") * 120


def generate_instance(instance_name, num_routes, param_digits, homo):
    '''
    Generate a test instance for the general optimal capacity allocation problem,
    and store it in a file.
        
    Parameters
    ----------
    instance_name: string
        The full name of the file include the path of the file
    num_routes: integer
        Numer of routes in the network
    param_digits: integer
        The number of decimal points kept for each parameter
    homo: binary
        Whether the instance is homogeneous or not
    
    Returns
    -------
    None.
    '''
    
    
    # Step 1: Generate an instance
    beta = np.round(np.random.uniform(0, 2), param_digits)
    phi = np.round(np.random.uniform(0, 2), param_digits)
    theta = 2
    if not homo:
        b = np.round(np.random.uniform(-4, 4, num_routes), param_digits)
        l = np.round(np.random.uniform(0, 1, num_routes), param_digits)
        u = np.round(np.random.uniform(l, 2, num_routes), param_digits)
    else:
        # Each route has the same value of b, l, and u.
        b = np.ones(num_routes) * np.round(np.random.uniform(-4, 4), param_digits)
        temp = np.round(np.random.uniform(0, 1), param_digits)
        l = np.ones(num_routes) * temp
        u = np.ones(num_routes) * np.round(np.random.uniform(temp, 2), param_digits)
    
    # Step 2: Determine if the instance is feasible   
    
    
    # Step 3: If the instance is feasible, then create a file for the instance and store the parameters
    with open(instance_name, 'w') as f1:
        if not homo:
            f1.write("Nonhomogeneous instance\n")
        else:
            f1.write("Homogeneous instance\n")
        f1.write("num_routes " + str(num_routes) + "\n")            
        f1.write("beta " + str(beta) + "\n")            
        f1.write("phi " + str(phi) + "\n")
        f1.write("theta " + str(theta) + "\n")
        f1.write("b\n")
        # Write a 1-d array of floats into one line of text file,
        # with the entries of the array separated by a comma
        f1.write(",".join([str(item) for item in b.tolist()]) + "\n")
        f1.write("l\n")
        f1.write(",".join([str(item) for item in l.tolist()]) + "\n")
        f1.write("u\n")
        f1.write(",".join([str(item) for item in u.tolist()]) + "\n")
        
    # Step 4: Generate an AMPL data file
    instance_name += '.dat'
    with open(instance_name, 'w') as f1:
        f1.write("param J := " + str(num_routes) + ";\n")
        f1.write("param beta :=" + str(beta) + ";\n")
        f1.write("param phi :=" + str(phi) + ";\n")
        f1.write("param theta :=" + str(theta) + ";\n\n")
        f1.write("param: b l u:=\n")
        for i in range(num_routes):
            f1.write(str(i+1) + " " + str(b[i]) + " " + str(l[i]) + " " + str(u[i]) + "\n")
        f1.write(";\n\n")
        f1.write("var q :=\n")
        for i in range(num_routes):
            f1.write(str(i+1) + " 0.001\n")
        f1.write(";")

        
        
        
       

def read_instance(instance_name):
    '''
    Read an instance file of the optimal capacity allocation problem, and create a class object to store all the parameters.
    
    Parameters
    ----------
    instance_name: string
        The full name of the file include the path of the file
    
    Returns
    -------
    None.
    '''
    ## Step 1: read the parameter of an instance    
    with open(instance_name) as f1:
        lines = f1.readlines()
    
    beta = float(lines[2].split()[1])
    phi = float(lines[3].split()[1])
    theta = int(lines[4].split()[1])
    
    # read an array from a text into a numpy array
    b = np.fromstring(lines[6], dtype=float, sep=',')
    l = np.fromstring(lines[8], dtype=float, sep=',')
    u = np.fromstring(lines[10], dtype=float, sep=',')

    ## Step 2: Create the instance
    instance = Network(num_routes, beta, phi, theta, b, l, u)
    return instance


def store_all_instance_name(num_routes, num_instances, homo):
    """
    Generate a string containing all instance names, which will be used in an AMPL batch command
    
    Parameters
    ----------
    num_routes: integer
        Number of the routes contained in the generated instances
    num_instances: integer
        Total number of instances generated
    homo: boolean
        True if the generated instance has homogeneous parameters
    
    Returns
    -------
    all_names: string
        A string containing all instance names
    """
    all_names = "{"
    if not homo:
        for instance in range(1, num_instances+1):
            all_names += "'general_r" + str(num_routes) + "_i" + str(instance) + ".dat', "            
    else:
        for instance in range(1, num_instances+1):
            all_names += "' homo_r'" + str(num_routes) + "_i" + str(instance) + ".dat', "
    all_names += "}"
    return all_names
    
    
  
##---------------------------Start the main program-------------------------
clear()

        

##-------------------------- Instance generation----------------------------

## Generate instances
# Initial setup
gen_homo_instance = False
param_digits = 3 # Precision of the randomly generated parameters
num_instances = 1
min_num_routes = 3
max_num_routes = 3


if not gen_homo_instance:
    instance_path = "./instances/general/"
else:
    instance_path = "./instances/homogeneous/"
    

for num_routes in range(min_num_routes, max_num_routes+1):
    for instance in range(1, num_instances+1):
        if not gen_homo_instance:
            instance_name = instance_path + 'general_r' + str(num_routes) + '_i' + str(instance)            
        else:
            instance_name = instance_path + 'homo_r' + str(num_routes) + '_i' + str(instance)            
        generate_instance(instance_name, num_routes, param_digits, gen_homo_instance)
    # Store the names of the instances in a string
    all_names = store_all_instance_name(num_routes, num_instances, gen_homo_instance)




##------------------------Running the base case S0--------------------------
## Solve the capacity design problem by considering both the effects of parking
## information and congestion effect


### Testing: check whether the optimal capacity is in the interior of the bounds
test_interior = np.empty((max_num_routes-min_num_routes+1, num_instances))
### End of testing

for num_routes in range(min_num_routes, max_num_routes+1):
    for instance in range(1, num_instances+1):
        ## Step 1: Read the instance from the file
        if not gen_homo_instance:
            filename = instance_path + 'general_r' + str(num_routes) + '_i' + str(instance)
        else:
            filename = instance_path + 'homo_r' + str(num_routes) + '_i' + str(instance)            
  
        
        ### Testing special instances
     
        network = read_instance(instance_path + 'general_r3_interior1')
        # Compute the value of the equilibrium flow for the homogeneous case when no parking information is available
        #q0 = bisection_search(zeta_homo_noparking, 1e-6, 1-1e-6, [1e-10, 1e-10], True, network)
        #print("The value of q0 is: " + str(q0))
        ###  End of testing special instances
        
        #network = read_instance(filename)
        
        ## Step 2: Optimize the capacity 
        
        opt_flow, opt_cap, opt_socialwelfare = search_optimal_capacities(network, 1e-4, 1e-10, filename)
        
         ##### Testing: Check whether there exists optimal capacity at the interior of the capacity bounds
        if np.amax(np.minimum(np.absolute(opt_cap - network.l), np.absolute(opt_cap - network.u))) > 0.1:
            test_interior[num_routes - min_num_routes, instance-1] = 1
        else:
            test_interior[num_routes - min_num_routes, instance-1] = 0
            
        ##### End of Testing
        
        # Close existing plots
        plt.close("all")
        
  
        ## Step 3: Record the results
        



##------------------------Running experiment S1-----------------------------





##------------------------Running experiment S2-----------------------------


##------------------------Running experiment S3-----------------------------



##------------------------Running experiment S4-----------------------------
