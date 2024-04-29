import os
import sys
sys.path.append(os.getcwd())  

import numpy as np
import TrussSolver as ts  # Import TrussSolver module

# Define Node Coordinates
nodes = np.array([    
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

# Define what is connecting to what
connections = [
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


t1 = ts.Truss(nodes, connections)

#Add Loads 
t1.addLoad(0, [0, 25])  #Support 1
t1.addLoad(4, [0, 25])  #Support 2
t1.addLoad(10, [0, -10])
t1.addLoad(11, [0, -10])
t1.addLoad(12, [0, -10])
t1.addLoad(13, [0, -10])
t1.addLoad(14, [0, -10])

print(t1.bars)

t1.solve()

t1.print_summary()

trussdraw = ts.PgTruss(t1, 1600)
trussdraw.drawNodes()
