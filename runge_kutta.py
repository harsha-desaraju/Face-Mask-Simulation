# This is a program to solve a second order differential equation
# using Runge-Kutta method of 1st order. for ode of kind
# x" = px'+f(x)

import matplotlib.pyplot as plt
from math import *
from random import uniform, gauss

def func(x):
    temp = 0
    #print(x,end = " ")
    #print("The function value is ",temp)
    return temp

def RungeKutta(x0,v0,t0,t,N,p,m,f):

    x = []
    v = []
    T = []

    dt = (t-t0)/N

    Vn = v0
    Xn = x0
    t1 =0
    while(t1<t):
        x.append(Xn)
        v.append(Vn)
        T.append(t1)
        #print(Vn,'\t',Xn)
        temp = gauss(0,t1)
        Vn_1 = Vn + dt*(p*Vn+f(Xn)+temp)*m
        #print(gauss(0,t1))
        Xn_1 = Xn +dt*Vn
        t1+=dt
        Vn = Vn_1
        Xn = Xn_1

    return [x,v,T]

x,vx,tx = RungeKutta(0,0,0,1,100,-1,0.01,func)
y,vy,ty = RungeKutta(0,0,0,1,100,-1,0.01,func)
#print(x)
#print(v)
#print(t)

#plt.plot(tx,vx,label="velocity")
#plt.plot(t,x,label="position" )
plt.plot(x,y)
plt.plot([0],[0],marker = '*')
plt.xlabel("x")
plt.ylabel("y")

plt.show()
