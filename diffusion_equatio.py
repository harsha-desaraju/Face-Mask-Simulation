# This is a program to solve the diffusion equation using finite difference method

from math import *
import matplotlib.pyplot as plt

class FiniteDifference:

    def generateMesh(self,T,Nt,L,Nx):
        ''' Generates a mesh from 0 to t with Nt uniform steps
           and from 0 to x with Nx uniform steps. Returns a list
           of lists with x horizontal and t vertical with all
           values as zeroes. '''

        self.T = T
        self.Nt = Nt
        self.L = L
        self.Nx = Nx

        mesh = []
        for i in range(Nt):
            mesh.append([])
            for j in range(Nx):
                mesh[i].append(0)

        self.mesh = mesh

    def boundaryConditions(self,f1,f2,f3):
        ''' f1 is a function which gives initial conditions for p(x,t=0)
            f2 is a function which gives boundary conditons for p(x=0,t)
            f3 is a function which gives boundary condtions for p(x=l,t)    '''

        delx = self.L/self.Nx
        delt = self.T/self.Nt

        for i in range(self.Nt):
            self.mesh[i][0] = f2(i*delt)


        for i in range(self.Nt):
            self.mesh[i][-1] = f3(i*delt)

        for i in range(self.Nx):
            self.mesh[0][i] = f1(i*delx,self.L,self.Nx)

    def finiteDifference(self,D):

        delx = self.L/self.Nx
        delt = self.T/self.Nt

        s = (D*delt)/(delx**2)

# The main algorithm for solving the PDE.
#These are the equations of finite difference method using forward propagation.

        for i in range(0,self.Nt-1):
            for j in range(1,self.Nx-1):
                 self.mesh[i+1][j] = self.mesh[i][j]*(1-2*s) + s*(self.mesh[i][j+1]+self.mesh[i][j-1])

        return self.mesh

def f1(x,L,Nx):
    if x == 0:
        return 1
    else:
        return 0

def f2(t):
    if t!=0:
        temp = 1/sqrt(4*pi*0.8*t)
        #print(t,'\t',temp)

        return temp

def f3(t):
    return 0

#Time
T = 10
#Length
L = 2
#Number of time steps
Nt = 100
#Number of steps in x
Nx = 50
#Diffusion constant
D = 0.005


def variance(lst,L,Nx):
    '''This is a function to calculate the variance in x at a given time.
        The function return the variances in x at different times as a list.'''
    var = []
    delx = L/Nx
    for i in range(len(lst)):
        temp = 0
        for j in range(len(lst[0])):
            temp = temp+ (j*delx)**2*lst[i][j]
        var.append(temp)

    return var

m = FiniteDifference()

m.generateMesh(T,Nt,L,Nx)
m.boundaryConditions(f1,f2,f3)
prob = m.finiteDifference(D)


var = variance(prob,L,Nx)

t = []
for i in range(Nt):
    delt = T/Nt
    t.append(i*delt)

line = []
for i in t:
    line.append(2*D*i)


plt.plot(t,var,label= 'calculated',marker= "*")
plt.plot(t,line, label='Theoretical')
plt.xlabel("time")
plt.ylabel("Variance")
plt.legend()
plt.show()
