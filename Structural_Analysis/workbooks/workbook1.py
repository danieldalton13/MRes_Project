import os
import sys
sys.path.append(os.getcwd())  

import numpy as np
import TrussSolver as ts  # Import TrussSolver module

# Load definitions with fixed nodes annotated
define_loads = [
    (0, [0, 25]),   # Support 1 (Fixed Node)
    (4, [0, 25]),   # Support 2 (Fixed Node)
    (10, [0, -10]), # Fixed Node
    (11, [0, -10]), # Fixed Node
    (12, [0, -10]), # Fixed Node
    (13, [0, -10]), # Fixed Node
    (14, [0, -10])  # Fixed Node
]

define_nodes = np.array([    
    [0,0],      #0
    [4,0],      #1
    [8,0],      #2
    [12,0],     #3
    [16,0],     #4
    [0,4],      #5
    [4,4],      #6
    [8,4],      #7
    [12,4],     #8
    [16,4],     #9
    [0,8],      #10
    [4,8],      #11
    [8,8],      #12
    [12,8],     #13
    [16,8],     #14
    ])

# Extract indices of fixed nodes from the load definitions
fixed_nodes_indices = [load[0] for load in define_loads]

# All node indices
all_node_indices = range(len(define_nodes))

# Identify variable nodes as those not in fixed_nodes_indices
variable_nodes_indices = [idx for idx in all_node_indices if idx not in fixed_nodes_indices]

# Separate fixed and variable nodes based on indices
fixed_nodes = define_nodes[fixed_nodes_indices]
variable_nodes = define_nodes[variable_nodes_indices]

# Node connections (topology of the truss)
define_connections = [
    [1,5,6],        # Node 0 connects to Nodes 1,5,6
    [0,2,6],                    # Node 1 
    [1,3,6,7,8],                # Node 2 
    [2,4,8],                    # Node 3  
    [3,8,9],                    # Node 4 
    [0,6,10],                   # Node 5  
    [0,1,2,5,7,10,11,12],       # Node 6  
    [2,6,8,12],                 # Node 7 
    [2,3,4,7,9,12,12,14],       # Node 8 
    [4,8,14],                   # Node 9  
    [5,6,11],                   # Node 10 
    [6,10,12],                  # Node 11  
    [6,7,8,11,13],              # Node 12  
    [8,12,14],                  # Node 13 
    [8,9,13]                    # Node 14 
]

# Setup Truss using TrussSolver
t1 = ts.Truss(define_nodes, define_connections)

# Adding loads to fixed nodes
for position, load in define_loads:
    t1.addLoad(position, load)

# Solve the truss system
t1.solve()

# Output the results
t1.print_summary()

# Drawing the truss configuration
trussdraw = ts.PgTruss(t1, 1600)
trussdraw.drawNodes()
