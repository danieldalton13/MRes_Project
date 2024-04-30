import numpy as np
import numpy.linalg as la
import math as mt
import pygame as pg

def checkInMembers(members, relation):
    i = relation[0]
    j = relation[1]
    output = False
    for n in members:
        if (i in n) and (j in n):
            output = True
    return output

class Truss:
    def __init__(self, nodes_, connections_):
        self.nodes = nodes_
        self.connections = connections_
        N = len(nodes_)
        self.loads = np.zeros(N*2)
        self.bars = []
        j = 0
        for connection in connections_:
            for i in connection:
                if not checkInMembers(self.bars, [i, j]):
                    self.bars.append([i, j])
            j += 1

        nForces = len(self.bars)
        self.G = np.zeros((2*N, nForces))
        self.F = np.zeros(len(self.bars))

        for i in range(2*N):
            if (i % 2 == 0):  # x-component
                iNode = int(i/2)
                for iForce in range(nForces):
                    if iNode in self.bars[iForce]:
                        bar = list(self.bars[iForce])
                        bar.remove(iNode)
                        j = bar[0]
                        self.G[i, iForce] = self.c(iNode, j)
            else:  # y-component
                iNode = int((i-1)/2)
                for iForce in range(nForces):
                    if iNode in self.bars[iForce]:
                        bar = list(self.bars[iForce])
                        bar.remove(iNode)
                        j = bar[0]
                        self.G[i, iForce] = self.s(iNode, j)

    def addLoad(self, node, load):
        iX = node*2
        iY = node*2 + 1
        self.loads[iX] = -load[0]
        self.loads[iY] = -load[1]

    def c(self, i, j):
        x0, y0 = self.nodes[i]
        x1, y1 = self.nodes[j]
        d = mt.sqrt((x0 - x1)**2 + (y0 - y1)**2)
        return (x1 - x0) / d if d != 0 else 0

    def s(self, i, j):
        x0, y0 = self.nodes[i]
        x1, y1 = self.nodes[j]
        d = mt.sqrt((x0 - x1)**2 + (y0 - y1)**2)
        return (y1 - y0) / d if d != 0 else 0

    def solve(self):
        lst = la.lstsq(self.G, self.loads, rcond=None)
        fmax = max(abs(force) for force in lst[0])
        fmax_tension = max(force for force in lst[0] if force > 0)
        self.result = lst[0]
        self.fmax = fmax
        self.fmax_tension = fmax_tension

    def truss_weight(self):
        density = 1
        total_weight = sum(density * mt.sqrt((self.nodes[bar[0]][0] - self.nodes[bar[1]][0])**2 + (self.nodes[bar[0]][1] - self.nodes[bar[1]][1])**2) for bar in self.bars)
        self.weight = total_weight
        self.ltw_tension = 1000 / (self.fmax_tension * self.weight)
        self.ltw = 1000 / (self.fmax * self.weight)
        return total_weight

    def print_summary(self):
        print("Truss Summary:")
        self.truss_weight()  # Ensure weight is calculated before printing
        print(f"Total Weight: {self.weight:.2f} units")
        print(f"Maximum Force in any member: {self.fmax:.2f} units")
        print(f"Maximum Tension Force in any member: {self.fmax_tension:.2f} units")

class PgTruss:
    def __init__(self, truss_, screenSize_):
        self.truss = truss_
        self.xmaxS = screenSize_
        self.margin = 100
        xmaxR = 0
        ymaxR = 0

        # Get maximum x and y range across all nodes
        x_coords = [node[0] for node in truss_.nodes]
        y_coords = [node[1] for node in truss_.nodes]
        if x_coords:
            xmaxR = max(x_coords) - min(x_coords)
        if y_coords:
            ymaxR = max(y_coords) - min(y_coords)

        self.xmaxR = xmaxR
        self.ymaxR = ymaxR

        # Avoid division by zero
        if xmaxR == 0:
            self.scale = 0  # or set a default scale if needed
        else:
            self.scale = (self.xmaxS - 2 * self.margin) / xmaxR
        
        self.ymaxS = int(self.scale * self.ymaxR + 2 * self.margin)
        pg.init()
        self.screen = pg.display.set_mode((self.xmaxS, self.ymaxS))

    def drawNodes(self):
        black, white, red, blue, green = (60, 60, 60), (255, 255, 255), (255, 147, 117), (149, 206, 252), (87, 255, 126)
        pg.font.init()
        myfont = pg.font.SysFont('Lucida Sans', 20)
        self.screen.fill(black)
        screenNodes = []
        nodeIndx = 0
        for node in self.truss.nodes:
            nodex = int(node[0] * self.scale + self.margin)
            nodey = int(self.ymaxS - (node[1] * self.scale + self.margin))
            screenNodes.append([nodex, nodey])
            pg.draw.circle(self.screen, white, (nodex, nodey), 4, 1)
            xloadI, yloadI = 2 * nodeIndx, 2 * nodeIndx + 1
            if self.truss.loads[xloadI] != 0 or self.truss.loads[yloadI] != 0:
                loadx = -self.truss.loads[xloadI]
                loady = -self.truss.loads[yloadI]
                text = f"P{nodeIndx} = {round(loady, 3)}"
                loadText = myfont.render(text, True, green)
                self.screen.blit(loadText, (nodex - 10, nodey - 25))
            textsurface = myfont.render(str(nodeIndx), True, white)
            self.screen.blit(textsurface, (nodex + 10, nodey + 5))
            nodeIndx += 1
        self.screenPoints = screenNodes
        barIndx = 0
        for bar in self.truss.bars:
            pointA, pointB = self.screenPoints[bar[0]], self.screenPoints[bar[1]]
            midPointx, midPointy = (pointA[0] + pointB[0]) // 2, (pointA[1] + pointB[1]) // 2
            forceBar = round(self.truss.result[barIndx], 2)
            color = red if forceBar > 0 else blue
            pg.draw.line(self.screen, white, pointA, pointB, 1)
            textsurface = myfont.render(str(forceBar), True, color)
            self.screen.blit(textsurface, (midPointx, midPointy - 20))
            barIndx += 1
        pg.display.flip()
        run = True
        while run:
            event = pg.event.wait()
            if event.type == pg.QUIT:
                run = False
            elif event.type == pg.KEYDOWN and event.key == pg.K_ESCAPE:
                run = False

# Use the Truss and PgTruss classes with node and connection definitions as follows
nodes = np.array([
    [0, 0], [4, 0], [8, 0], [12, 0], [16, 0],
    [0, 4], [4, 4], [8, 4], [12, 4], [16, 4],
    [0, 8], [4, 8], [8, 8], [12, 8], [16, 8],
])

connections = [
    [1, 5, 6], [0, 2, 6], [1, 3, 6, 7, 8], [2, 4, 8],
    [3, 8, 9], [0, 6, 10], [0, 1, 2, 5, 7, 10, 11, 12],
    [2, 6, 8, 12], [2, 3, 4, 7, 9, 12, 12, 14], [4, 8, 14],
    [5, 6, 11], [6, 10, 12], [6, 7, 8, 11, 13], [8, 12, 14],
    [8, 9, 13]
]

t1 = Truss(nodes, connections)
t1.addLoad(0, [0, 25])
t1.addLoad(4, [0, 25])
t1.addLoad(10, [0, -10])
t1.addLoad(11, [0, -10])
t1.addLoad(12, [0, -10])
t1.addLoad(13, [0, -10])
t1.addLoad(14, [0, -10])

t1.solve()
t1.truss_weight()
t1.print_summary()

trussdraw = PgTruss(t1, 1600)
trussdraw.drawNodes()
