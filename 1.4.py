# This code calculates the average efficiency of a mask with thread diameter
# and thread pitch given as ranges.

import random
from math import sqrt
import matplotlib.pyplot as plt
from scipy.stats import *

def normal(mean, sd, low, upp):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

def exponential(mean, sd, low, upp):
    return truncexpon(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)


def min_or_same(a,b):
    if b == diameter_range[0]:
        return b
    else:
        if a<b:
            return a
        else:
            return b

def norm_parameters(a,b):
    lst = []
    lst.append((a+b)/2)
    lst.append((b-a)/10)
    lst.append(a)
    lst.append(b)
    return lst

def generate_grid(lst1,lst2):
    temp_lst = norm_parameters(diameter_range[0],diameter_range[1])
    #print(temp_lst)
    temp1 = normal(temp_lst[0],temp_lst[1],temp_lst[2],temp_lst[3]).rvs()
    #print(temp1.rvs())
    lst1.append(temp1)
    while(lst1[-1]<=total_length):
        temp_lst1= norm_parameters(pitch_range[0],pitch_range[1])
        temp2 = exponential(temp_lst1[1],temp_lst1[2],temp_lst1[3]).rvs()
        #print(temp2.rvs())
        temp_lst = norm_parameters(diameter_range[0],min_or_same(temp2,diameter_range[1]))
        temp1 = normal(temp_lst[0],temp_lst[1],temp_lst[2],temp_lst[3]).rvs()
        #print(temp1.rvs())
        lst2.append(temp2-temp1)
        lst1.append(lst1[-1]+lst2[-1]+temp1)

def generate_particle():
    temp1 = x_coordinates[-1]-x_coordinates[-2]-x_pitch[-1]
    x = random.uniform(x_coordinates[0],x_coordinates[-1]-temp1/2)
    temp2 = y_coordinates[-1]-y_coordinates[-2]-y_pitch[-1]
    y = random.uniform(y_coordinates[0],y_coordinates[-1]-temp2/2)
    return [x,y]

def coordinates(lst,num):
    low,high = 0,len(lst)-1
    while(high-low>0):
        n = int((low+high)/2)
        #print(low,'\t',high,'\t',n,'\t',num)
        if lst[n]<num and lst[n+1]>num:
            return [n,n+1]
        else:
            if lst[n]<=num:
                low = n
            else:
                high = n

def condition(d,R,r):
    if d > (R+r)/sqrt(2):
        return True
    else:
        return False

def passes(lst1,lst2,num,rad):
    temp = coordinates(lst1,num)
    #print(temp,"\t",num)
    thread2 = lst1[temp[1]]-lst1[temp[0]]-lst2[temp[0]]
    thread1 = lst1[temp[0]]-lst1[temp[0]-1]- lst2[temp[0]-1]
    d1 = num - lst1[temp[0]] + (thread1/2)
    d2 = lst1[temp[1]] - num - (thread2/2)
    if condition(d1,thread1/2,rad) and condition(d2,thread2/2,rad):
        return True
    else:
        return False

def particle_passes(particle,rad):
    if passes(x_coordinates,x_pitch,particle[0],rad) and passes(y_coordinates,y_pitch,particle[1],rad):
        return True
    else:
        False

def filtration_efficiency(rad):
    passed = 0
    for i in range(N):
        particle = generate_particle()
        if particle_passes(particle,rad):
            passed+=1
    return ((N-passed)/N)*100

def porosity():
    temp = sum(y_pitch)
    free_area = 0
    for i in x_pitch:
        free_area = free_area + (i*temp)
    return (free_area/(x_coordinates[-1]*y_coordinates[-1]))*100

def min_max(lst,j):
    max,min = lst[0][j],lst[0][j]
    for i in range(len(final)):
        if lst[i][j]<min:
            min = lst[i][j]
        elif lst[i][j]>max:
            max = lst[i][j]
    ranges[0].append(min)
    ranges[1].append(max)

def calc_error():
    for i in range(len(efficiencies)):
        ranges[0][i] = efficiencies[i]-ranges[0][i]
        ranges[1][i] = ranges[1][i]-efficiencies[i]


total_length = float(input("Enter the total length of the mask (in cm)\n"))*1e-2
diameter_range = [float(ele)*1e-6 for ele in input("Enter a range for thread diameter (in um)\n").split('-')]
pitch_range = [float(ele)*1e-6 for ele in input("Enter a range for thread pitch (in um)\n").split("-")]
radii = list(float(ele)*1e-9 for ele in input("Enter some radii seperated by comma (in nm)\n").split(','))
N = int(input("Enter the number of particles for each radii\n"))
T = int(input("Enter the number of runs for averaging\n"))

x_coordinates = []
y_coordinates = []
x_pitch = []
y_pitch = []

generate_grid(x_coordinates,x_pitch)
generate_grid(y_coordinates,y_pitch)

final = []
for i in range(T):
    filtration_efficiencies =[]
    for radius in radii:
        filtration_efficiencies.append(filtration_efficiency(radius))
    final.append(filtration_efficiencies)

efficiencies = []
ranges = [[],[]]
for j in range(len(radii)):
    temp = 0
    min_max(final,j)
    for i in range(T):
        temp = temp+final[i][j]
    efficiencies.append(temp/T)

calc_error()

print("The porosity of the mask is ",porosity())
print(radii)
print(filtration_efficiencies)

plt.plot(radii,efficiencies)
plt.errorbar(radii,efficiencies,yerr=ranges,capsize = 5)

plt.xlabel("Radii")
plt.ylabel("filtration_efficiency")
plt.show()
