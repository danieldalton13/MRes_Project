import os
import sys
sys.path.append(os.getcwd())

import numpy as np
import TrussSolver as ts  # Import TrussSolver module

# Load and node definitions
define_loads = [
    (0, [0, 25]),
    (4, [0, 25]),
    (12, [0, -50])
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

# Indices of fixed nodes from the load definitions
fixed_nodes_indices = [load[0] for load in define_loads]

best_weight = float('inf')
best_configuration = None

# Function to check if every fixed node has at least one connection
def valid_configuration(connections, fixed_nodes_indices):
    node_connection_counts = {index: 0 for index in fixed_nodes_indices}
    for node_index in fixed_nodes_indices:
        for connection in connections:
            if node_index in connection:
                node_connection_counts[node_index] += 1
    return all(count > 0 for count in node_connection_counts.values())

# Test every possible combination of connections
from itertools import combinations

for L in range(1, len(define_connections)+1):
    for subset in combinations(define_connections, L):
        test_connections = [list(conn) for conn in subset]
        
        if not valid_configuration(test_connections, fixed_nodes_indices):
            continue

        # Setup Truss using TrussSolver
        t1 = ts.Truss(define_nodes, test_connections)
        
        # Adding loads to fixed nodes
        for position, load in define_loads:
            t1.addLoad(position, load)
        
        # Run the truss solver and capture output
        try:
            t1.solve()
            t1.truss_weight()
            weight = t1.weight
            
            if weight < best_weight:
                best_weight = weight
                best_configuration = test_connections
                
        except Exception as e:
            print(f"Failed configuration: {e}")

# Reinitialize with the best configuration to summarize and draw
t1 = ts.Truss(define_nodes, best_configuration)
for position, load in define_loads:
    t1.addLoad(position, load)
t1.solve()

# Output the results
t1.print_summary()

# Drawing the truss configuration
trussdraw = ts.PgTruss(t1, 1600)
trussdraw.drawNodes()
