# This is code which includes both interception and electric effects.
# PLEASE GO THROUGH README FILE

from scipy.stats import *
from scipy.integrate import quad
from random import *
from math import *
import matplotlib.pyplot as plt


# The class Particle() expects a list of contents of file as an argument.
# That is the list is the list using readlines().
# Particle object has 5 attributes.
class Particle:

    def __init__(self,lst):
        PropList1 = lst

        prl = float(str(str(PropList1[7]).split("=")[1]).split('-')[0])  # prl - particle radius lower
        pru = float(str(str(PropList1[7]).split("=")[1]).split('-')[1])  # prl - particle radius upper
        l = float(str(PropList1[0]).split("=")[1])*1e4
        b = float(str(PropList1[1]).split("=")[1])*1e4
        h = float(str(PropList1[2]).split("=")[1])*1e4
        psl = float(str(str(PropList1[8]).split("=")[1]).split('-')[0])  # prl - particle spped lower
        psu = float(str(str(PropList1[8]).split("=")[1]).split('-')[1])  # prl - particle spped upper

# Particles are generated using normal distribution
        def normal(mean, sd, low, upp):
            return truncnorm(
                (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd).rvs()

        self.radius = normal((prl+pru)/2,(prl+pru)/3,prl,pru)
# The position of particles is generated using uniform distribution.
        self.position = [uniform(0,l),uniform(0,b),uniform(h/2,h)]
# Speed of particle along z-axis.
        self.speed = uniform(psl,psu)*1e6
# Charge of particle in nano coulombs
        self.charge = float(str(PropList1[9]).split("=")[1])
# Mass of particle in micro grams
        self.mass = float(str(PropList1[10]).split("=")[1])


class Mask:

    def __init__(self,f):
        self.PropList = f.readlines()

        self.length = float(str(self.PropList[0]).split("=")[1])*1e4
        self.breadth = float(str(self.PropList[1]).split("=")[1])*1e4
        self.twl = float(str(str(self.PropList[3]).split("=")[1]).split('-')[0])  # twl - thread width lower
        self.twu = float(str(str(self.PropList[3]).split("=")[1]).split('-')[1])  # twu - thread width upper
        self.tpl = float(str(str(self.PropList[4]).split("=")[1]).split('-')[0])  # tpl - thread pitch lower
        self.tpu = float(str(str(self.PropList[4]).split("=")[1]).split('-')[1])  # tpu - thread pitch upper
        self.chargeDensity = float(str(self.PropList[5]).split("=")[1])*1e-6
        self.scale = float(str(self.PropList[6]).split("=")[1])
        self.number = int(str(self.PropList[11]).split("=")[1])

# This is the distribution of thread widths and thread pitchs.
        def normal(mean, sd, low, upp):
            return truncnorm(
                (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd).rvs()

# This is a function to generate grid.
        def generateGrid():
            #list containing coordinates(tuple) of starting and ending of threads
            grid = []
            temp = normal((self.twu+self.twl)/2,(self.twu-self.twl)/3,self.twl,self.twu) # 1st thread generated
            grid.append((0,temp))
            while(grid[-1][1] < self.length):
                temp1 = normal((self.tpu+self.tpl)/2,(self.tpu-self.tpl)/3,self.tpl,self.tpu)
                temp2 = max(temp1,self.twu)
                temp = normal((temp2+self.twl)/2,(temp2-self.twl)/3,self.twl,temp2)
                grid.append((grid[-1][1]+temp1-temp,grid[-1][1]+temp1))
            return grid

        self.xCoordinates = generateGrid()
        self.yCoordinates = generateGrid()
        self.length = self.xCoordinates[-1][1]
        self.breadth = self.yCoordinates[-1][1]

    def throwParticles(self):
        self.particle = Particle(self.PropList)
        self.passed = 0
        self.electricallyPassed = 0

        poin = []
        for i in self.xCoordinates:
            poin.append((i[0]+i[1])/2)

        def threadIndex(lst,a):
            low,high = 0,len(lst)-1

            while(high-low>0):
                mid = int((low+high)/2)
                if lst[mid][0] < a and lst[mid+1][1] > a:
                    return mid
                else:
                    if lst[mid][0]<a:
                        low = mid
                    else:
                        high = mid

        def points(lst,lst1,a,b):
            return lst1[threadIndex(lst,a)+1:threadIndex(lst,b)]


        xInd = threadIndex(self.xCoordinates,self.particle.position[0])#∞∞∞∞
        yInd = threadIndex(self.yCoordinates,self.particle.position[1])

        def notIntercepted(p):
            if self.xCoordinates[xInd][1]+p.radius < p.position[0] and self.xCoordinates[xInd+1][0]-p.radius > p.position[0]:
                if self.yCoordinates[yInd][1]+p.radius < p.position[1] and self.yCoordinates[yInd+1][0]-p.radius > p.position[1]:
                    return True
                else:
                    return False
            else:
                return False

#potential due to grid at a point (x,y,z)
        def V(x,y,z):
            #print(x  ,y  ,z)
            #potential due to 'a' vertical threads
            def Vv(a,x,y,z):
                n1 = sqrt((self.breadth-y)**2+(a-x)**2+z**2)+ self.breadth-y
                d1 = sqrt(y**2 + (a-x)**2 +z**2) - y
                return log(n1/d1)

            #potential due to 'a' horizontal threads
            def Vh(a,x,y,z):
                n2 = sqrt((self.length-x)**2+(a-y)**2+z**2)+self.length-x
                d2 = sqrt(x**2+(a-y)**2+z**2) - x
                return log(n2/d2)


            low ,high = 0 , len(self.xCoordinates)-1
            temp1 =0
            for i in range(low,high):
                a = (self.xCoordinates[i][0]+self.xCoordinates[i][1])/2
                temp1 = temp1 + Vv(a,x,y,z)

            low ,high = 0 , len(self.yCoordinates)-1
            temp2 =0
            for i in range(low,high):
                a = (self.yCoordinates[i][0]+self.yCoordinates[i][1])/2
                temp2 = temp2 + Vh(a,x,y,z)

            return temp1+temp2

        const = V(self.particle.position[0],self.particle.position[1],self.particle.position[2])
        factor = 2*self.particle.charge*self.chargeDensity*9*1e18/self.particle.mass

        def FindTz(p):

            factor = 2*p.charge*self.chargeDensity*9*1e18/p.mass

            def integrate(z):
                return (p.speed**2 - factor*(V(p.position[0],p.position[1],z) - const))**(-1/2)

            return quad(integrate,0,p.position[2])[0]

        tz = FindTz(self.particle)

        def dx(p):

            def vel(p,pos):
                delV = p.charge*(V(p.position[0],p.position[1],p.position[2])-V(pos[0],pos[1],pos[2]))
                dV = V(p.position[0],p.position[1],p.position[2]) - V(p.position[0]+0.01,p.position[1],p.position[2])

                if delV < 0:
                    return 0
                else:
                    if dV>0:
                        sign = -1
                    else:
                        sign = 1
                    v = sign*(2*(delV)/p.mass)**(1/2)

                    return v

            dt = tz/50
            if V(p.position[0],p.position[1],p.position[2]) - V(p.position[0]+0.01,p.position[1],p.position[2]) < 0:
                pre = p.position[0]+0.01
                v = vel(p,[p.position[0]+0.01,p.position[1],p.position[2]])
            else:
                pre = p.position[0]-0.01
                v = vel(p,[p.position[0]-0.01,p.position[1],p.position[2]])

            t=0
            while(t<=tz):
                nxt = pre + v*dt
                t+=dt
                v = vel(p,[nxt,p.position[1],p.position[2]])
                pre = nxt

            return nxt

        def dy(p):
            def vel(p,pos):
                delV = p.charge*(V(p.position[0],p.position[1],p.position[2])-V(pos[0],pos[1],pos[2]))
                dV = V(p.position[0],p.position[1],p.position[2]) - V(p.position[0],p.position[1]+0.01,p.position[2])

                if delV < 0:
                    return 0
                else:
                    if dV>0:
                        sign = -1
                    else:
                        sign = 1
                    v = sign*(2*(delV)/p.mass)**(1/2)

                    return v

            dt = tz/50
            if V(p.position[0],p.position[1],p.position[2]) - V(p.position[0],p.position[1]+0.01,p.position[2]) < 0:
                pre = p.position[1]+0.01
                v = vel(p,[p.position[0],p.position[1]+0.01,p.position[2]])
            else:
                pre = p.position[1]-0.01
                v = vel(p,[p.position[0],p.position[1]-0.01,p.position[2]])

            t=0
            while(t<=tz):
                nxt = pre + v*dt
                t+=dt
                v = vel(p,[p.position[0],nxt,p.position[2]])
                pre = nxt

            return nxt



        def notElectricallyIntercepted(p,lst):
            if self.xCoordinates[xInd][1]+p.radius < lst[0] and self.xCoordinates[xInd+1][0]-p.radius > lst[0]:
                if self.yCoordinates[yInd][1]+p.radius < lst[1] and self.yCoordinates[yInd+1][0]-p.radius > lst[1]:
                    return True
                else:
                    return False
            else:
                return False



        for i in range(self.number):
            #print("The value of dx using Euler method is ",dx(self.particle))
            #print("The value of dy using Euler mathod is ",dy(self.particle))
            if notIntercepted(self.particle):
                self.passed+=1
            if notElectricallyIntercepted(self.particle,[dx(self.particle),dy(self.particle)]):
                self.electricallyPassed+=1
            self.particle = Particle(self.PropList)
            xInd = threadIndex(self.xCoordinates,self.particle.position[0])
            yInd = threadIndex(self.yCoordinates,self.particle.position[1])
            const = V(self.particle.position[0],self.particle.position[1],self.particle.position[2])

    def printAll(self):
        print(self.xCoordinates, len(self.xCoordinates))
        print('\n')
        print(self.yCoordinates,len(self.yCoordinates))

    def filtrationEfficiency(self):
        print("The interception filtration efficiency is ",(1-(self.passed/self.number))*100)
        print("The electric filtration efficiency is ",(1-(self.electricallyPassed/self.number))*100)

    def TheoreFiltraEffi(self):
# This function calculates the ratio of occupied area to total area of the mask.
        def TheoFilteffi():
            temp =0
            for i in range(len(self.xCoordinates)):
                temp = temp + self.xCoordinates[i][1] - self.xCoordinates[i][0]
            temp = temp*self.breadth

            temp1 = 0
            for i in range(len(self.yCoordinates)):
                temp2 = 0
                for j in range(len(self.xCoordinates)-1):
                    temp2 = temp2 + self.xCoordinates[j+1][0] - self.xCoordinates[j][1]
                temp1 = temp1 + temp2*(self.yCoordinates[i][1]-self.yCoordinates[i][0])

            return ((temp+temp1)/(self.length*self.breadth))*100

        print("The theoretical efficiency is ",TheoFilteffi())

    def visualise(self):
        self.potential = []
        self.x = []
        self.particle = Particle(self.PropList)

        def V(a,x,y,z):
            n1 = sqrt((self.breadth-y)**2+(a-x)**2)+ self.breadth-y
            d1 = sqrt(y**2 + (a-x)**2 ) - y
            n2 = sqrt((self.length-x)**2+(a-y)**2)+self.length-x
            d2 = sqrt(x**2+(a-y)**2) - x
            return log(n1/d1)+log(n2/d2)


        def Fx(x):
            low ,high = 0 , len(self.xCoordinates)-1
            temp = 0
            for i in range(low,high):
                temp = temp + V((self.xCoordinates[i][0]+self.xCoordinates[i][1])/2,x,self.breadth/2,5)
            return temp

        def Fy(x):
            low ,high = 0 , len(self.yCoordinates)-1
            temp = 0
            for i in range(low,high):
                temp = temp + V((self.yCoordinates[i][0]+self.yCoordinates[i][1])/2,x,self.breadth/2,5)
            return temp

        def V(x,y,z):
            #potential due to 'a' vertical threads
            def Vv(a,x,y,z):
                n1 = sqrt((self.breadth-y)**2+(a-x)**2+z**2)+ self.breadth-y
                d1 = sqrt(y**2 + (a-x)**2 +z**2) - y
                return log(n1/d1)

            #potential due to 'a' horizontal threads
            def Vh(a,x,y,z):
                n2 = sqrt((self.length-x)**2+(a-y)**2+z**2)+self.length-x
                d2 = sqrt(x**2+(a-y)**2+z**2) - x
                return log(n2/d2)


            low ,high = 0 , len(self.xCoordinates)-1
            temp1 =0
            for i in range(low,high):
                a = (self.xCoordinates[i][0]+self.xCoordinates[i][1])/2
                temp1 = temp1 + Vv(a,x,y,z)

            low ,high = 0 , len(self.yCoordinates)-1
            temp2 =0
            for i in range(low,high):
                a = (self.yCoordinates[i][0]+self.yCoordinates[i][1])/2
                temp2 = temp2 + Vh(a,x,y,z)

            return temp1+temp2


        x = 0
        dx = 10

        a = []
        Va = []
        for i in self.xCoordinates:
            a.append((i[0]+i[1])/2)
            Va.append(V(a[-1],self.breadth/2,5))
            #print(self.particle.charge,'\t',end = " ")
            #print(self.particle.charge*V(a[-1],self.breadth/2,5))

        while(x<self.length):
            self.potential.append(self.particle.charge*V(x,self.breadth/2,5))
            self.x.append(x)
            x= x+dx

        plt.plot(self.x,self.potential)
        #plt.plot(a,Va,marker='*')
        plt.xlabel("x")
        plt.ylabel("potential energy")

        plt.show()











f = open("2.2 Mask.props.txt",'r')
m = Mask(f)
#m.printAll()
#m.throwParticles()
#m.filtrationEfficiency()
#m.TheoreFiltraEffi()
m.visualise()
f.close()
