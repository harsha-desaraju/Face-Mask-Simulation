import random
import matplotlib.pyplot as plt
from math import *

print("Running.....")

def generate_particle():
    return [random.uniform(0,total_length),random.uniform(0,total_length)]

def in_between(tup,num):
    if tup[0]<num and num<tup[1]:
        return True
    else:
        return False

def condition_satisfied(tup,num,particle_radius):
    d1 = num-tup[0]
    d2 = tup[1]-num
    if d1>(fibre_radius+particle_radius)/sqrt(2) and d2>(fibre_radius+particle_radius)/sqrt(2):
        return True
    else:
        return False

def particle_passes(particle,particle_radius):
    x,y = False,False
    for tup in fibre_coordinates:
        if in_between(tup,particle[0]):
            if condition_satisfied(tup,particle[0],particle_radius):
                x = True

    for tup in fibre_coordinates:
        if in_between(tup,particle[1]):
            if condition_satisfied(tup,particle[1],particle_radius):
                y = True
    return x and y

def filtration_efficiency(rad):
    passed = 0
    for i in range(N):
        test_particle = generate_particle()
        if particle_passes(test_particle,rad):
            passed+=1
    return ((N-passed)/N)*100



total_length = float(input("Enter the total length of the grid\n"))
fibre_width = float(input("Enter the width of the fibre\n"))
pore_size = float(input("Enter the size of the pore\n"))
N = int(input("Enter the number of simulations\n"))

n = round((total_length-fibre_width)/(pore_size+fibre_width))
total_length = n*pore_size+(n+1)*fibre_width

fibre_radius = fibre_width/2

fibre_coordinates = []
for i in range(n):
    coordinates = (i*pore_size+(2*i+1)*fibre_radius,(i+1)*pore_size+(2*(i+1)+1)*fibre_radius)
    fibre_coordinates.append((coordinates))

test_rad = float(input("Enter a radius for particles\n"))

print("The percentage of particles blocked is ",filtration_efficiency(test_rad))
