import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from pylab import *
from scipy.optimize import fsolve
from scipy.optimize import minimize

# ---------------------------------------------


# make local files accessible to import statements
import sys, os
sys.path.append(os.path.join(os.getcwd(), '.'))
from source.mesh_factory import Regular2D
from source.constants import *
from source.w import *


class Matrices:

    def __init__(self, mesh, mi, la):
        self.m = mesh
        self.W11 = w11(mesh)
        self.W12 = w12(mesh)
        self.W22 = w22(mesh)
        self.W21 = w12(mesh).T

        # print("self.W11")
        # print(self.W11)
        # print("self.W12")
        #  print(self.W12)
        # print("self.W22")
        # print(self.W22)

        self.B11 = (2 * mi + la) * self.W11 + mi * self.W22
        self.B12 = mi * self.W21 + la * self.W12
        self.B21 = la * self.W21 + mi * self.W12
        self.B22 = mi * self.W11 + (2 * mi + la) * self.W22

        for pointId in range(len(self.m.point)):
            if (mesh.point.dirichlet[pointId]):
                self.B11[pointId][pointId] = 1
                self.B12[pointId][pointId] = 1
                self.B21[pointId][pointId] = 1
                self.B22[pointId][pointId] = 1


# ---------------------------------------------
class Drawer:

    def __init__(self, solv):
        self.solv = solv
        self.mesh = solv.mesh

    def draw(self):
        txt = 'C:\\Users\\user\\PycharmProjects\\ContactMechanics\\source\\output\\CROSS OPTIMIZATION NEW'

        plt.axes().set_aspect('equal', 'box')
        plt.scatter(self.mesh.point[:, 0], mesh.point[:, 1])

        for edge in self.mesh.edge.points:
            x1, y1 = self.mesh.point[int(edge[0])][0], self.mesh.point[int(edge[0])][1]
            x2, y2 = self.mesh.point[int(edge[1])][0], self.mesh.point[int(edge[1])][1]
            plt.plot([x1, x2], [y1, y2], 'k-', lw=0.5)

        txtpng = txt + '.png'

        plt.scatter(self.solv.DisplacedPoints[:, 0], self.solv.DisplacedPoints[:, 1], marker='o')

        plt.savefig(txtpng, transparent=True, bbox_inches='tight', pad_inches=0, dpi=300)  # DPI 2000
        plt.close()

    # plt.scatter(self.solv.DisplacedPoints[:,0],self.solv.DisplacedPoints[:,1], marker='o')

# ---------------------------------------------
class Solver:

    def __init__(self, mesh, F0, FN, mi, la):

        self.mesh = mesh
        self.numberOfPoints = len(mesh.point)
        self.mi = mi
        self.la = la

        self.M = Matrices(mesh, mi, la)
        self.F = F(mesh, F0, FN)

        self.u = np.zeros([self.numberOfPoints, 2])
        self.DisplacedPoints = np.zeros([self.numberOfPoints, 2])

        for i in range(0, self.numberOfPoints):
            self.DisplacedPoints[i] = self.mesh.point.coordinates[i]

    def norm(self, v):
        return sqrt(v[0] * v[0] + v[1] * v[1])

    def length(self, p1, p2):
        return float(sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1])))

    '''
    def nDown(self, e):
        #[0,-1]
        e1 = int(self.s.Edges[e][0])
        e2 = int(self.s.Edges[e][1])
        dx = self.s.Points[e2][0]-self.s.Points[e1][0]
        dy = self.s.Points[e2][1]-self.s.Points[e1][1]
        norm = sqrt(dx*dx + dy*dy)
        n = np.array([float(dy)/norm, float(-dx)/norm])
        if(n[1]> 0):
            n = -n
        return n

    def nDownAtContactBoundary(self):
        N = np.zeros([self.s.indNumber(),2])  

        for i in range (0,self.s.indNumber()):
            for e in range (-self.s.BorderEdgesD-self.s.BorderEdgesN -self.s.BorderEdgesC, -self.s.BorderEdgesD-self.s.BorderEdgesN):
                e1 = int(self.s.Edges[e][0])
                e2 = int(self.s.Edges[e][1])
                if(i == e1 or i == e2):
                    n = self.nDown(e)
                    if(N[i][0] == 0 and N[i][1] == 0):
                        N[i] = n
                    else:
                        N[i] = 0.5*N[i] + 0.5*n
        return N
    '''

    def Jwu(self):
        J = float(0.)
        return J
        '''
        for i in range (0,self.s.indNumber()):
            for e in range (-self.s.BorderEdgesD-self.s.BorderEdgesN -self.s.BorderEdgesC, -self.s.BorderEdgesD-self.s.BorderEdgesN):
                e1 = int(self.s.Edges[e][0])
                e2 = int(self.s.Edges[e][1])
                if(i == e1 or i == e2):
                    umL =  np.zeros(2) # u at mL
                    wmL =  np.zeros(2) # w at mL
                    if(e1 < self.s.indNumber()):
                        umL += self.u[e1]*0.5
                        wmL += self.w[e1]*0.5
                    if(e2 < self.s.indNumber()):
                        umL += self.u[e2]*0.5
                        wmL += self.w[e2]*0.5

                    p1 = self.s.Points[int(e1)][0:2]
                    p2 = self.s.Points[int(e2)][0:2]
                    mL = (p1 + p2) * 0.5
                    L = self.length(p1,p2)
                    nmL = self.nDown(e)  #n at mL

                    uNmL = umL[0] * nmL[0] + umL[1] * nmL[1]
                    wNmL = wmL[0] * nmL[0] + wmL[1] * nmL[1]
                    uTmL = umL - uNmL * nmL

                    J += L * 0.5 * (self.jN(uNmL) + self.h(wNmL) * self.jT(uTmL)) ###########

        return J
        '''

    '''
    def setWVector(self, wVector):
        self.w[:,0] = wVector[0:self.s.indNumber()]
        self.w[:,1] = wVector[self.s.indNumber():2*self.s.indNumber()]  
    '''

    def iterate(self, uVector):
        print("self.u")
        print(self.u)
        self.u[:, 0] = uVector[0:self.numberOfPoints]
        self.u[:, 1] = uVector[self.numberOfPoints:2 * self.numberOfPoints]

        for i in range(0, self.numberOfPoints):
            self.DisplacedPoints[i][0] = self.mesh.point.coordinates[i][0] + self.u[i][0]
            self.DisplacedPoints[i][1] = self.mesh.point.coordinates[i][1] + self.u[i][1]

    def fu(self):
        Fu = float(0.)
        for i in range(0, self.numberOfPoints):
            Fu += np.dot(self.F.F[i], self.u[i, :])
        return Fu

    def Buu(self):
        return np.dot((np.dot(self.M.B11, self.u[:, 0]) + np.dot(self.M.B12, self.u[:, 1])), self.u[:, 0]) \
               + np.dot((np.dot(self.M.B21, self.u[:, 0]) + np.dot(self.M.B22, self.u[:, 1])), self.u[:, 1])

    def L2(self, uVector):
        self.u[:, 0] = uVector[0:self.numberOfPoints]
        self.u[:, 1] = uVector[self.numberOfPoints:2 * self.numberOfPoints]

        l = (0.5 * float(self.Buu())) + float(self.Jwu()) - float(self.fu())

        return l * 100000

    ########################################################

    knu = 10.
    delta = 0.1

    def jN(self, uN):  # uN - scalar
        if (uN <= 0):
            return 0
        if (uN <= self.delta):
            return 0.5 * self.knu * uN * uN
        # return self.knu * self.delta * (uN - 0.5 * self.delta)
        return 0.5 * self.knu * self.delta * self.delta

    def h(self, uN):
        if (uN <= 0):
            return 0
        return 8. * uN
        # return 16.*uN

    def jT(self, uT):  # uT - vector
        # return self.norm(uT)
        return log(self.norm(uT) + 1)

########################################################


# ---------------------------------------------
class F:

    def __init__(self, mesh, F0, FN):
        self.F0 = F0
        self.FN = FN
        self.mesh = mesh
        self.numberOfPoints = len(mesh.point)
        self.F = np.zeros([self.numberOfPoints, 2])
        self.Zero = np.zeros([self.numberOfPoints])
        self.One = np.zeros([self.numberOfPoints])

    ########################################################

    def f0(self, x):
        return self.F0

    def fN(self, x):
        return self.FN

    ########################################################

    def setF(self):
        # halfLongTriangleSide = self.mesh.halfLongTriangleSide
        # quaterLongTriangleSide = self.mesh.halfLongTriangleSide/ 2

        self.F = np.zeros([self.numberOfPoints, 2])

        for i, element in enumerate(mesh.element):
            a = mesh.edge[element[EDGE_0]][0]
            b = mesh.edge[element[EDGE_0]][1]
            c = mesh.edge[element[EDGE_1]][0]
            a_ok = not mesh.point.dirichlet[a]
            b_ok = not mesh.point.dirichlet[b]
            c_ok = not mesh.point.dirichlet[c]
            xa, ya = mesh.point[a][0], mesh.point[a][1]
            xb, yb = mesh.point[b][0], mesh.point[b][1]
            xc, yc = mesh.point[c][0], mesh.point[c][1]

            efield = 1 / 4 * mesh.edge.length[element[EDGE_0]] ** 2

            if (a_ok):
                self.F[a] += (efield / 6.) * (
                            self.f0([(xa + xb) * 0.5, (ya + yb) * 0.5]) + self.f0([(xa + xc) * 0.5, (ya + yc) * 0.5]))
            if (b_ok):
                self.F[b] += (efield / 6.) * (
                            self.f0([(xb + xa) * 0.5, (yb + ya) * 0.5]) + self.f0([(xb + xc) * 0.5, (yb + yc) * 0.5]))
            if (c_ok):
                self.F[c] += (efield / 6.) * (
                            self.f0([(xc + xa) * 0.5, (yc + ya) * 0.5]) + self.f0([(xc + xb) * 0.5, (yc + yb) * 0.5]))

        print("self.F Z 0")
        print(self.F)

        for border in mesh.subarea['border'][NEUMANN]:
            for edgeId in border:
                p1Id = mesh.edge.points[edgeId][0]
                p2Id = mesh.edge.points[edgeId][1]

                x1, y1 = mesh.point[p1Id][0], mesh.point[p1Id][1]
                x2, y2 = mesh.point[p2Id][0], mesh.point[p2Id][1]

                edgeLength = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                if (not mesh.point.dirichlet[p1Id]):
                    self.F[p1Id] += ((edgeLength * 0.5) * self.fN([(x1 + x2) * 0.5, (y1 + y2) * 0.5]))
                if (not mesh.point.dirichlet[p2Id]):
                    self.F[p2Id] += ((edgeLength * 0.5) * self.fN([(x1 + x2) * 0.5, (y1 + y2) * 0.5]))

        self.Zero = self.F[:, 0]
        self.One = self.F[:, 1]

# ---------------------------------------------
########################################################
F0 = np.array([0.5, 0.4])  # np.array([-1.2,-1.]) -0.9
FN = np.array([-0.3, 0.])
mi = 4
la = 4
########################################################

mesh = Regular2D.construct(1, 0.25, left=DIRICHLET, top=NEUMANN, right=NEUMANN, bottom=CONTACT)
# print(mesh.point.coordinates)
# print(mesh.subarea[""]["dirichlet"])
matrices = Matrices(mesh, 4, 4)

solver = Solver(mesh, F0, FN, mi, la)
solver.F.setF()
drawer = Drawer(solver)

uVector = np.zeros([2 * solver.numberOfPoints])

while True:
    uVector = minimize(solver.L2, uVector, method='Powell' \
                       , options={'disp': True, 'xtol': 0.000000001,
                                  'ftol': 0.00000001}).x  # 'xtol': 0.000000001, 'ftol': 0.00000001 #'Powell' 'BFGS' 'neldermead'
    print(uVector)
    # print(".", end=" ")
    # print("FIXED DIFF NORM: " + str(np.linalg.norm(np.subtract(uVector,wVector))))
    # if (np.linalg.norm(np.subtract(uVector,wVector)) < 0.00001): #0.00001
    break

solver.iterate(uVector)
drawer.draw()
print("DONE")


