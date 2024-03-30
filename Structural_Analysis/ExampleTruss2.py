import os
import sys
sys.path.append(os.getcwd())  # Ensure the current working directory is in the sys.path

import numpy as np
import TrussSolver as ts  # Import your custom TrussSolver module

# Directly defining nodes using numpy array
nodes = np.array([
    [0, 0],    # Node 0
    [80, 0],   # Node 1, assuming direct input based on previous L[0]
    [130, 0],  # Node 2, assuming direct input based on previous L[1]
    [250, 0],  # Node 3, directly using 'c'
    [370, 0],  # Node 4, example direct input, adjust as necessary
    [420, 0],  # Node 5, example direct input, adjust as necessary
    [500, 0],  # Node 6, 2*c
    [80, 30],  # Node 7, direct input based on previous H[0] and L[0]
    [130, 50], # Node 8, direct input based on previous H[1] and L[1]
    [250, 50], # Node 9, using 'c' and height
    [370, 50], # Node 10, example, adjust as necessary
    [420, 30]  # Node 11, example, adjust as necessary
])

# Keep connections as a list of lists, not converting to numpy array
connections = [
    [1, 7],             # 0
    [0, 7, 2],          # 1
    [1, 7, 8, 3],       # 2
    [2, 8, 9, 10, 4],   # 3
    [3, 10, 11, 5],     # 4
    [4, 11, 6],         # 5
    [5, 11],            # 6
    [0, 1, 2, 8],       # 7
    [7, 2, 3, 9],       # 8
    [8, 3, 10],         # 9
    [9, 3, 4, 11],      # 10
    [10, 4, 5, 6]       # 11
]

t1 = ts.Truss(nodes, connections)

t1.addLoad(0, [0, 0.5])
t1.addLoad(6, [0, 0.5])
t1.addLoad(9, [0, -1])

print(t1.bars)

t1.solve()

t1.print_summary()

trussdraw = ts.PgTruss(t1, 1600)
trussdraw.drawNodes()
