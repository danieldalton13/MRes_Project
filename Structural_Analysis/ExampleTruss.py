import os
import sys
sys.path.append(os.getcwd())  # Ensure the current working directory is in the sys.path

import numpy as np
import TrussSolver as ts  # Import your custom TrussSolver module

# Truss dimensions and specifications
c = 250
side = c/3
height = 50
h1 = 30
h2 = 20
l1 = 80
l2 = 50

L = [l1, l1+l2, c]
H = [h1, h1 + h2, height]

# Defining nodes using numpy array
nodes = np.array([    
    [0, 0],             # 0
    [L[0], 0],          # 1
    [L[1], 0],          # 2
    [L[2], 0],          # 3
    [2*c - L[1], 0],    # 4
    [2*c - L[0], 0],    # 5
    [2*c, 0],           # 6
    [L[0], H[0]],       # 7
    [L[1], H[1]],       # 8
    [L[2], H[2]],       # 9
    [2*c - L[1], H[1]], # 10
    [2*c - L[0], H[0]]  # 11
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
