import os
import sys
sys.path.append(os.getcwd())

import numpy as np
import TrussSolver as ts  # Import TrussSolver module

# Load and node definitions
define_loads = [
    (0, [0, 25]),
    (4, [0, 25]),
    (10, [0, -10]),
    (11, [0, -10]),
    (12, [0, -10]),
    (13, [0, -10]),
    (14, [0, -10])
]

define_nodes = np.array([
    [0, 0],
    [4, 0],
    [8, 0],
    [12, 0],
    [16, 0],
    [0, 4],
    [4, 4],
    [8, 4],
    [12, 4],
    [16, 4],
    [0, 8],
    [4, 8],
    [8, 8],
    [12, 8],
    [16, 8],
])

# Initial connection setup
define_connections = [
    [1, 5, 6],
    [0, 2, 6],
    [1, 3, 6, 7, 8],
    [2, 4, 8],
    [3, 8, 9],
    [0, 6, 10],
    [0, 1, 2, 5, 7, 10, 11, 12],
    [2, 6, 8, 12],
    [2, 3, 4, 7, 9, 12, 13, 14],
    [4, 8, 14],
    [5, 6, 11],
    [6, 10, 12],
    [6, 7, 8, 11, 13],
    [8, 12, 14],
    [8, 9, 13]
]

# Extract indices of fixed nodes from the load definitions
fixed_nodes_indices = [load[0] for load in define_loads]

# All node indices
all_node_indices = range(len(define_nodes))

# Identify variable nodes as those not in fixed_nodes_indices
variable_nodes_indices = [idx for idx in all_node_indices if idx not in fixed_nodes_indices]

# Separate fixed and variable nodes based on indices
fixed_nodes = define_nodes[fixed_nodes_indices]
variable_nodes = define_nodes[variable_nodes_indices]

best_weight = float('inf')
best_configuration = None

# Iterate over each node's connection list
for i in range(len(define_connections)):
    for j in range(len(define_connections[i])):
        # Create a modified connection list for this test
        test_connections = [list(conn) if k != i else conn[:j] + conn[j+1:] for k, conn in enumerate(define_connections)]
        
        # Setup Truss using TrussSolver
        t1 = ts.Truss(define_nodes, test_connections)
        
        # Adding loads to fixed nodes
        for position, load in define_loads:
            t1.addLoad(position, load)
        
        # Run the truss solver and capture output
        try:
            t1.solve()
            t1.truss_weight()  # This will calculate and update the weight attribute
            weight = t1.weight
            
            if weight < best_weight:
                best_weight = weight
                best_configuration = test_connections
                
        except Exception as e:
            print(f"Failed configuration with node {i} connection {define_connections[i][j]} removed: {e}")


# Output the results
t1.print_summary()

# Drawing the truss configuration
trussdraw = ts.PgTruss(t1, 1600)
trussdraw.drawNodes()
