import numpy as np
import numpy.linalg as la 
import math as mt 
import pygame as pg 

## 'checkInMembers' checks a given pair of nodes (def by i and j) ##
## form a member that is not alraedy listed in the 'members' list ##

def checkInMembers(members, relation):  # Function named 'checkInMembers' takes two parameters 'members' and relation'
    i = relation[0] # assigns the first elemnt of relation and assigns it to 'i'
    j = relation[1] # assigns the second elemnt of relation and assigns it to 'j'
    output = False # initialises the output function to False, thus if a new member then result will be 'False'
    for n in members: # loops through members to see if any matches
        if (i in n) and (j in n):
            output = True # if the new member is already found in the list then no "new" member found thus 'True'
    return output



class Truss:
    def __init__(self,nodes_, connections_):  #constructor takes 'nodes_' and 'connections_' as inputs and assigns them to 'self.nodes' and 'self.connections' respectively
        self.nodes = nodes_
        self.connections = connections_

        N = len(nodes_) #determines the number of nodes

        self.loads = np.zeros(N*2) #initialises a zero array for loads applied to nodes. n=2 as loads have x and y component
        self.bars = [] #initialises an empty list to store the bars (members)
        j = 0
        for connection in connections_: #the nested loops iterate over each connection, checking with 'checkInMembers' if a bar between nodes 'i' and 'j' is already listed in 'self.bars'. If not, the bar is added. 
            for i in connection:
                if not checkInMembers(self.bars,[i,j]):
                    self.bars.append([i,j])
            j += 1
        nForces = len(self.bars)
        self.G = np.zeros((2*N,nForces)) #global stiffness matrix initialisation - matrix with zeros, with 2*N rows (for x and y components at each node) and nForces columns (each representing a bar)
        self.F = np.zeros(len(self.bars)) #force vector initialisation - zero array for internal forces in bars, sized by the number of bars

        for i in range(2*N): #the loop constructs the global stiffness matrix by iterating over each potential force (represented by bars) and determining its contribution to each nodes equilibrium.
            if (i%2 == 0): #checks if the current index 'i' is even, which corresponds to an x-compoent of force at a node (because of the arrangement in the 'loads' array where x and y components alternate)
                iNode = int(i/2) #calculates the actual node index from 'i'. eg, if 'i' is 0, or 2 'iNode' will be 0, indicating its dealing with the x component of the first nodes force
                for iForce in range(nForces): #iterates over each each potential force in the truss, which corresponds to each bar or membeer. 
                    if iNode in self.bars[iForce]:  #checks if the current node is part of the bar being considered
                        bar = list(self.bars[iForce])   #used to identify the other node 'j' connected to 'iNode' by the current bar. 
                        bar.remove(iNode)
                        j = bar[0]
                        self.G[i,iForce] = self.c(iNode,j) #calculates the cosine of the angle between the bat and the x-axis
                    else:
                        self.G[i,iForce] = 0
            else : #corresponds to the odd indicies of 'i' which deal with the y-component at each force
                iNode = int((i-1)/2) 
                for iForce in range(nForces):
                    if iNode in self.bars[iForce]:
                        bar = list(self.bars[iForce])
                        bar.remove(iNode)
                        j = bar[0]
                        self.G[i,iForce] = self.s(iNode,j)
                    else:
                        self.G[i,iForce] = 0
                                           

    def addLoad(self,node, load): #apply external loads to the nodes of the truss structure
        iX = node*2 #calculates the index of the 'self.loads' array where th x-component of the force for the specified node should be stored
        iY = node*2 + 1
        self.loads[iX] = -load[0]   #applies the load
        self.loads[iY] = -load[1]
    
    def c(self, i, j):
        x0 = self.nodes[i][0]
        x1 = self.nodes[j][0]
        y0 = self.nodes[i][1]
        y1 = self.nodes[j][1]

        d = mt.sqrt( (x0 - x1)**2 + (y0 - y1)**2 )

        if d == 0:
            return 0
        else:
            return (x1 - x0)/d

    def s(self, i, j):
        x0 = self.nodes[i][0]
        x1 = self.nodes[j][0]
        y0 = self.nodes[i][1]
        y1 = self.nodes[j][1]

        d = mt.sqrt( (x0 - x1)**2 + (y0 - y1)**2 )
        if d == 0:
            return 0
        else:
            return (y1 - y0)/d
    
    def k(self, i, j):
        if j in self.connections[i]:
            return 1
        else:
            return 0

    def solve(self):

        lst = la.lstsq(self.G, self.loads , rcond=None)
        fmax = 0
        fmax_tension = 0
        for force in lst[0]:
            if abs(force) > fmax:
                fmax = abs(force)
        for force in lst[0]:
            if force > fmax_tension:
                fmax_tension = force

        self.result = lst[0]
        self.fmax = fmax
        self.fmax_tension = fmax_tension


    def truss_weight(self):
        density = 1
        total_weight = 0
        for bar in self.bars:
            total_weight += density * mt.sqrt( (  self.nodes[bar[0]][0] - self.nodes[bar[1]][0] )**2 + (  self.nodes[bar[0]][1] - self.nodes[bar[1]][1]  )**2 )
        self.weight = total_weight
        self.ltw_tension = 1000 * 1/(self.fmax_tension * self.weight)
        self.ltw = 1000 * 1/(self.fmax * self.weight)
        return total_weight

    def print_summary(self):
        print("Truss Summary:")
        self.truss_weight()  # Make sure this calculates and stores weight in self.weight
        print(f"Total Weight: {self.weight:.2f} units")  # Assuming weight is calculated in some units

        # Assuming self.solve() has been called and it sets self.fmax and self.fmax_tension
        print(f"Maximum Force in any member: {self.fmax:.2f} units")
        print(f"Maximum Tension Force in any member: {self.fmax_tension:.2f} units")

        # If you have other features to print, add them here
    

class PgTruss:
    def __init__(self,truss_,screenSize_):
        self.truss = truss_
        self.xmaxS = screenSize_
        self.margin = 100
        xmaxR = 0
        ymaxR = 0
        for nodei in truss_.nodes:
            for nodej in truss_.nodes:
                if abs(nodei[0] - nodej[0]) > xmaxR:
                    xmaxR = abs(nodei[0] - nodej[0])
                if abs(nodei[1] - nodej[1]) > ymaxR:
                    ymaxR = abs(nodei[1] - nodej[1])

        self.xmaxR = xmaxR
        self.ymaxR = ymaxR 
        self.scale =(self.xmaxS - 2*self.margin)/xmaxR
        self.ymaxS = int((self.scale * self.ymaxR) + self.margin*2)
        pg.init()
        self.screen = pg.display.set_mode((self.xmaxS,self.ymaxS))

    def drawNodes(self):
        
        black = (60,60,60)
        white = (255,255,255)
        red = (255, 147, 117)
        blue = (149, 206, 252)
        green = (87, 255, 126)
        pg.font.init()
        myfont = pg.font.SysFont('Lucida Sans', 20)
        self.screen.fill(black)
        screenNodes = []
        nodeIndx = 0
        for node in self.truss.nodes:
            nodex = int((node[0] * self.scale) + self.margin)
            nodey = int(self.ymaxS -((node[1] * self.scale) + self.margin))
            screenNodes.append([nodex,nodey])
            pg.draw.circle(self.screen,white,(nodex,nodey),4, 1)
            
            xloadI = 2*nodeIndx
            yloadI = xloadI + 1
            if self.truss.loads[xloadI] != 0 or self.truss.loads[yloadI] != 0 :
                loadx = 0 - self.truss.loads[xloadI]
                loady = 0 - self.truss.loads[yloadI]
                text = " P" + str(nodeIndx) + " = " + str(round(loady,3)) 
                loadText = myfont.render(text, True, green)
                self.screen.blit(loadText,(nodex - 10, nodey - 25))
            
            textsurface = myfont.render(str(nodeIndx), True, white)
            self.screen.blit(textsurface,(nodex + 10, nodey + 5))

            nodeIndx += 1
            
            
        self.screenPoints = screenNodes

        barIndx = 0
        for bar in self.truss.bars:
            pointA = self.screenPoints[bar[0]]
            pointB = self.screenPoints[bar[1]]
            midPointx = int((pointA[0] + pointB[0])/2)
            midPointy = int((pointA[1] + pointB[1])/2)

            forceBar = round(self.truss.result[barIndx],2)
            if forceBar>0:
                color = red
            else:
                color = blue

            pg.draw.line(self.screen,white,pointA, pointB, 1)

            textsurface = myfont.render(str(forceBar), True, color)
            self.screen.blit(textsurface,(midPointx, midPointy - 20))
            barIndx += 1
        
        
        
        
        pg.display.flip()
        run = True
        while run:
            event = pg.event.wait()
            if event.type == pg.QUIT:
                run = False  
            elif event.type == pg.KEYDOWN and event.key == pg.K_ESCAPE:
                run = False
