from math import gcd, ceil                                                          # Imports gcd for checking connectivity, ceil for rounding up in 'member adding' 
import itertools                                                                    # Used for combinations in creating the 'ground structure'
from scipy import sparse                                                            # Sparse matrix utilities for efficient memory usage in matrix operations
import numpy as np                                                                  # Fundamental package for scientific computing with Python
import cvxpy as cvx                                                                 # For solving optimization problems
import matplotlib.pyplot as plt                                                     # For plotting the truss structures
from shapely.geometry import Point, LineString, Polygon                             # Geometric objects for checking intersections and domain shapes

#Calculate matrix B (the equilibrium matrix)
def calcB(Nd, Cn, dof):
    m, n1, n2 = len(Cn), Cn[:,0].astype(int), Cn[:,1].astype(int)                   # Extract node indices and the number of connections
    l, X, Y = Cn[:,2], Nd[n2,0]-Nd[n1,0], Nd[n2,1]-Nd[n1,1]                         # Lengths and direction components of each member
    d0, d1, d2, d3 = dof[n1*2], dof[n1*2+1], dof[n2*2], dof[n2*2+1]                 # DOF (Degrees of Freedom) flags for the ends of each member
    s = np.concatenate((-X/l * d0, -Y/l * d1, X/l * d2, Y/l * d3))                  # Stress contributions from each DOF
    r = np.concatenate((n1*2, n1*2+1, n2*2, n2*2+1))                                # Row indices for the sparse matrix (equilibrium constraints)
    c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m)))    # Column indices for the sparse matrix
    return sparse.coo_matrix((s, (r, c)), shape = (len(Nd)*2, m))                   # Return the sparse matrix in COO (Coordinate list) format

#Solve linear programming problem
def solveLP(Nd, Cn, f, dof, st, sc, jc):
    l = [col[2] + jc for col in Cn]                                                 # Adjusted member lengths incorporating joint costs
    B = calcB(Nd, Cn, dof)                                                          # Calculate the equilibrium matrix B
    a = cvx.Variable(len(Cn))                                                       # Decision variable for cross-sectional areas
    obj = cvx.Minimize(np.transpose(l) * a)                                         # Objective function to minimize volume
    q, eqn, cons= [],  [], [a>=0]                                                   # Force variables, equilibrium equations, and non-negativity constraints
    for k, fk in enumerate(f):
        q.append(cvx.Variable(len(Cn)))                                             # Force variable for each load case
        eqn.append(B * q[k] == fk * dof)                                            # Equilibrium constraint for each load case
        cons.extend([eqn[k], q[k] >= -sc * a, q[k] <= st * a])                      # Adding constraints for stress limits
    prob = cvx.Problem(obj, cons)                                                   # Define the CVXPY (convex optimisation) problem
    vol = prob.solve()                                                              # Solve the optimization problem
    q = [np.array(qi.value).flatten() for qi in q]                                  # Retrieve force distribution for each load case
    a = np.array(a.value).flatten()                                                 # Retrieve optimal cross-sectional areas
    u = [-np.array(eqnk.dual_value).flatten() for eqnk in eqn]                      # Retrieve dual values (used for checking violations)
    return vol, a, q, u

# Check dual violation (ensuring the solution respects the dual constraint)
def stopViolation(Nd, PML, dof, st, sc, u, jc):
    lst = np.where(PML[:,3]==False)[0]                                              # Indices of inactive members
    Cn = PML[lst]                                                                   # Extract inactive members
    l = Cn[:,2] + jc                                                                # Apply joint costs
    B = calcB(Nd, Cn, dof).tocsc()                                                  # Calculate equilibrium matrix for inactive members
    y = np.zeros(len(Cn))                                                           # Initialize violation measure
    for uk in u:
        yk = np.multiply(B.transpose().dot(uk) / l, np.array([[st], [-sc]]))        # Calculate stress ratios
        y += np.amax(yk, axis=0)                                                    # Max stress ratio across all load cases
    vioCn = np.where(y>1.0001)[0]                                                   # Find members where the stress ratio exceeds the limit
    vioSort = np.flipud(np.argsort(y[vioCn]))                                       # Sort these members by violation magnitude
    num = ceil(min(len(vioSort), 0.05*max( [len(Cn)*0.05, len(vioSort)])))          # Determine number of members to activate
    for i in range(num):
        PML[lst[vioCn[vioSort[i]]]][3] = True                                       # Activate the most violated members
    return num == 0                                                                 # Return True if no members were activated, indicating no violations

#Visualize truss
def plotTruss(Nd, Cn, a, q, threshold, str, update = True):
    plt.ion() if update else plt.ioff()                                             # Interactive mode on for updates, off for final plot
    plt.clf(); plt.axis('off'); plt.axis('equal');  plt.draw()                      # Clear the plot and set up the axes
    plt.title(str)                                                                  # Title of the plot
    tk = 5 / max(a)                                                                 # Thickness scaling based on maximum area
    for i in [i for i in range(len(a)) if a[i] >= threshold]:                       # Loop over all members above the area threshold
        if all([q[lc][i]>=0 for lc in range(len(q))]): c = 'r'                      # Color red if all internal forces are tensile
        elif all([q[lc][i]<=0 for lc in range(len(q))]): c = 'b'                    # Color blue if all internal forces are compressive
        else: c = 'tab:gray'                                                        # Color gray if forces vary across load cases
        pos = Nd[Cn[i, [0, 1]].astype(int), :]                                      # Position of the member's nodes
        plt.plot(pos[:, 0], pos[:, 1], c, linewidth = a[i] * tk)                    # Draw the member
    plt.pause(0.01) if update else plt.show()                                       # Update or show the plot

#Main function
def trussopt(width, height, st, sc, jc):
    poly = Polygon([(0, 0), (width, 0), (width, height), (0, height)])              # Define the domain as a simple rectangle
    convex = True if poly.convex_hull.area == poly.area else False                  # Check if the domain is convex
    xv, yv = np.meshgrid(range(width+1), range(height+1))                           # Create a grid of nodes within the domain
    pts = [Point(xv.flat[i], yv.flat[i]) for i in range(xv.size)]                   # Convert grid points to Shapely points
    Nd = np.array([[pt.x, pt.y] for pt in pts if poly.intersects(pt)])              # Keep only points within the polygon
    dof, f, PML = np.ones((len(Nd),2)), [], []                                      # Initialize degree of freedom and force arrays, and the potential member list (PML)

    #Load and support conditions
    for i, nd in enumerate(Nd):
        if nd[0] == 0: dof[i,:] = [0, 0]                                            # Set zero degrees of freedom for nodes on the left boundary (supports)
        f += [0, -1] if (nd == [width, height/2]).all() else [0, 0]                 # Apply a downward force at the midpoint of the right boundary

    #Create the 'ground structure'
    for i, j in itertools.combinations(range(len(Nd)), 2):
        dx, dy = abs(Nd[i][0] - Nd[j][0]), abs(Nd[i][1] - Nd[j][1])                 # Calculate horizontal and vertical distances between nodes
        if gcd(int(dx), int(dy)) == 1 or jc != 0:                                   # Check for connectivity using GCD unless joint costs modify connectivity
            seg = [] if convex else LineString([Nd[i], Nd[j]])                      # Create a line segment between nodes if domain is not convex
            if convex or poly.contains(seg) or poly.boundary.contains(seg):         # Check if the line is within the domain
                PML.append( [i, j, np.sqrt(dx**2 + dy**2), False] )                 # Add the member to the potential list, initially inactive
    PML, dof = np.array(PML), np.array(dof).flatten()                               # Convert lists to numpy arrays for better performance
    f = [f[i:i+len(Nd)*2] for i in range(0, len(f), len(Nd)*2)]                     # Format the force array properly for optimization
    print('Nodes: %d Members: %d' % (len(Nd), len(PML)))                            # Print the current number of nodes and members in the ground structure
    for pm in [p for p in PML if p[2] <= 1.42]:                                     # Activate members whose length is less than or equal to 1.42 units
        pm[3] = True

    #Start the 'member adding' loop
    for itr in range(1, 100):                                                       # Limit iterations to 100 for practicality
        Cn = PML[PML[:,3] == True]                                                  # Filter active members
        vol, a, q, u = solveLP(Nd, Cn, f, dof, st, sc, jc)                          # Solve the linear programming problem
        print("Itr: %d, vol: %f, mems: %d" % (itr, vol, len(Cn)))                   # Print current iteration results
        plotTruss(Nd, Cn, a, q, max(a) * 1e-3, "Itr:" + str(itr))                   # Plot the current truss configuration
        if stopViolation(Nd, PML, dof, st, sc, u, jc): break                        # Check for and handle dual violations, break if none
    print("Volume: %f" % (vol))                                                     # Print the final volume after optimization
    plotTruss(Nd, Cn, a, q, max(a) * 1e-3, "Finished", False)                       # Display the final optimized truss structure

#Execution function when called directly by Python
if __name__ =='__main__':
    trussopt(width = 20, height = 10, st = 1, sc =1, jc = 0)                        # Call the main function with specified parameters
