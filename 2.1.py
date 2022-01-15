# This code incorporates the electric effects in the filtration process of a mask.

from math import *
import random
from scipy.integrate import quad


class Mask:

    def __init__(self,f):
        PropList = f.readlines()
        self.length  = float(str(PropList[0]).split('=')[1])*1e4
        self.bredth  = float(str(PropList[1]).split('=')[1])*1e4
        self.threadWidth = float(str(PropList[2]).split('=')[1])
        self.threadGap = float(str(PropList[3]).split('=')[1])
        self.chargeDensity = float(str(PropList[4]).split('=')[1])*1e-6    ######### Take care of units
        self.particleRadius = float(str(PropList[5]).split('=')[1])
        self.particleNum = int(str(PropList[6]).split('=')[1])
        self.z = float(str(PropList[7]).split('=')[1])*1e4
        self.ParticleMass = float(str(PropList[8]).split('=')[1])
        self.particleCharge = float(str(PropList[9]).split('=')[1])
        self.particleSpeed = float(str(PropList[10]).split('=')[1])*1e6

        def generateGrid(dist):
            temp = []
            temp.append(self.threadWidth)
            while(temp[-1]<dist):
                temp.append(temp[-1]+self.threadGap+self.threadWidth)
            return temp

        self.xCoordinates = generateGrid(self.length)
        self.yCoordinates = generateGrid(self.bredth)
        self.length = self.xCoordinates[-1]
        self.bredth = self.yCoordinates[-1]

    def throwParticles(self):
        self.passed = 0   # Stores the number of particles that passed the mask.

#The below function  generates a particle in xy plane.
#In this code z is fixed and is taken from props file.
        def particle():
            x = random.uniform(0,self.length)
            y = random.uniform(0,self.bredth)
            return [x,y]

#The below function is for finding the fibers above whom the particle appears.
#The function uses bisection method.
        def coor(lst,a):
            low,high = 0,len(lst)-1
            while(high-low>0):
                mid = int((low+high)/2)
                if lst[mid]< a and lst[mid+1]>a:
                    return [lst[mid],lst[mid+1]]
                else:
                    if lst[mid] < a:
                        low = mid
                    else:
                        high = mid


#The below function checks if the particle is blocked by process of interception.
#If the particle touches the fiber then it does not pass.
        def notIntercepted(temp):
            lst1 = coor(self.xCoordinates,temp[0])
            lst2 = coor(self.yCoordinates,temp[1])
            if temp[0]-lst1[0] > self.particleRadius and lst1[1]-self.threadWidth-temp[0] > self.particleRadius:
                if temp[1]-lst2[0] > self.particleRadius and lst2[1] - self.threadWidth-temp[1] >self.particleRadius:
                    return True
                else:
                    return False
            else:
                return False


        def f(a,x,z):
            return log(sqrt(abs(a**2 - x**2))+z**2)

        def FindT(lst,p):
            i = lst.index(coor(lst,p)[0])

            if i-50<0:
                low,high = 0,i+50
            elif i+50>len(lst)-1:
                low,high = i-50,len(lst)-1
            else:
                low,high = i-50,i+50

            def F(x):        #p is for position which means it is its coordinate.
                temp = 0
                for j in range(low,high):
                    temp = temp+ f(p-lst[j]-self.threadWidth/2,x,self.z)
                return temp

            dx = min(p-lst[i],lst[i+1]-p)     ############ Possiblity for list index going out of range

            const = F(dx)

            def integrate(x):
                return sqrt((high-low)*const - F(x))

            factor = (self.ParticleMass/(2*self.particleCharge*self.chargeDensity*18*1e18))**(1/2)

            return factor*quad(integrate,dx,0)[0]


        def FindTz(lst):

            i = self.xCoordinates.index(coor(self.xCoordinates,temp[0])[0])
            j = self.yCoordinates.index(coor(self.yCoordinates,temp[1])[0])

            if i-10<0:
                low,high = 0,i+10
            elif i+10>len(self.xCoordinates)-1:
                low,high = i-10,len(self.xCoordinates)-1
            else:
                low,high = i-10,i+10


            def Fx(z):
                temp =0
                for k in range(low,high):
                    temp = temp+ f(lst[0]-self.xCoordinates[k]-self.threadWidth/2,lst[0],z)
                return temp

            const1 = (high-low)*Fx(self.z)

            if j-10<0:
                low,high = 0,j+10
            elif j+10>len(self.yCoordinates)-1:
                low,high = j-10,len(self.yCoordinates)-1
            else:
                low,high = j-10,j+10


            def Fy(z):
                temp =0
                for k in range(low,high):
                    temp = temp+ f(lst[1]-self.yCoordinates[k]-self.threadWidth/2,lst[1],z)
                return temp


            const2 = (high-low)*Fy(self.z)

            factor = (2*self.particleCharge*self.chargeDensity*18*1e18)/self.ParticleMass
            #print(factor)

            def integrate(z):
                #print((self.particleSpeed**2 - factor*(Fx(z)+Fy(z)-const1-const2))**(-1/2))
                return (self.particleSpeed**2 - factor*(Fx(z)+Fy(z)-const1-const2))**(-1/2)

            return quad(integrate,self.z,0)[0]

        def notElectricallyIntercepted(temp):
            #print(FindTz(temp) ,FindT(self.xCoordinates,temp[0]) ,FindT(self.yCoordinates,temp[1]))
            if FindTz(temp) < min(FindT(self.xCoordinates,temp[0]),FindT(self.yCoordinates,temp[1])):
                return True
            else:
                return False

        for i in range(self.particleNum):
            temp =particle()
            if notIntercepted(temp) and notElectricallyIntercepted(temp):
                    self.passed+=1

    def filtrationEfficiency(self):
        print("The number of particle to pass ",self.passed)
        return (1-(self.passed/self.particleNum))*100


f = open("2.1 Mask.props.txt","r")
n = Mask(f)
n.throwParticles()
print("The filtration efficiency of the mask is ",n.filtrationEfficiency())
f.close()
