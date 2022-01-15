# This code calculates the average efficiency of filtration of a square grid mask
# with finite thread width.

import random
import matplotlib.pyplot as plt
from math import *

print("Running.....")

def generate_particle():
    x= random.uniform(thread_radius,total_length-thread_radius)
    y= random.uniform(thread_radius,total_length-thread_radius)
    return [x,y]

def coordinates(num):
    low,high = 0,len(thread_coordinates)-1
    while (high-low>1):
        n = round((low+high)/2)
        if thread_coordinates[n][1]<num:
            low =n
        else:
            if num<thread_coordinates[n][0]:
                high = n
            else:
                return thread_coordinates[n]

    if low==0:
        return thread_coordinates[low]
    else:
        return thread_coordinates[high]

def condition_satisfied(tup,num,particle_radius):
    d1 = num-tup[0]
    d2 = tup[1]-num
    if d1>num_cond and d2>num_cond:
        return True
    else:
        return False

def particle_passes(particle,particle_radius):
    x,y = False,False
    x_coordinates = coordinates(particle[0])
    if condition_satisfied(x_coordinates,particle[0],particle_radius):
        x = True
    y_coordinates = coordinates(particle[1])
    if condition_satisfied(y_coordinates,particle[1],particle_radius):
        y = True
    return x and y

def filtration_efficiency(rad):
    passed = 0
    for i in range(N):
        test_particle = generate_particle()
        if particle_passes(test_particle,rad):
            passed+=1
    return ((N-passed)/N)*100

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

total_length = float(input("Enter the total length of the mask (in cm)\n"))
thread_diameter = float(input("Enter the diameter of the thread (in um)\n"))
pore_size = float(input("Enter the size of the pore (in um)\n"))
N = int(input("Enter the number of particles of each radius\n"))

total_length = total_length*1e-2
thread_diameter = thread_diameter*1e-6
pore_size = pore_size*1e-6

n = round((total_length-thread_diameter)/(pore_size+thread_diameter))
total_length = n*pore_size+(n+1)*thread_diameter

# Constants for the program
thread_radius = thread_diameter/2

thread_coordinates = []
for i in range(n):
    temp_coordinates = (i*pore_size+(2*i+1)*thread_radius,(i+1)*pore_size+(2*(i+1)+1)*thread_radius)
    thread_coordinates.append((temp_coordinates))

radii = list(float(ele)*1e-9 for ele in input("Enter some radii seperated by comma (in um)\n").split(','))
T = int(input("Enter the number of runs for average\n"))

final = []
for i in range(T):
    filtration_efficiencies = []
    for radius in radii:
        num_cond = (thread_radius+radius)/sqrt(2)
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

porosity = (n*pore_size/total_length)**2*100

print("The porosity of the mask is ",porosity,"\n")
print("The radii are \n",radii,"\n")
print("The filtration efficiencies of radii are \n",efficiencies,"\n")
plt.plot(radii,efficiencies)
plt.errorbar(radii,efficiencies,yerr=ranges,capsize = 5)

plt.xlabel("Radii")
plt.ylabel("filtration_efficiency")
plt.show()
