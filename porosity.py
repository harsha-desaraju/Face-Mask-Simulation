# With base unit of length as micro meter

import random
from scipy.stats import *
import matplotlib.pyplot as plt

def min_or_same(a,b):
    if b == thread_diameter[0]:
        return b
    else:
        if a<b:
            return a
        else:
            return b

def normal(mean, sd, low, upp):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd).rvs()

def norm_parameters(a,b):
    lst = []
    lst.append((a+b)/2)
    lst.append((b-a)/8)
    lst.append(a)
    lst.append(b)
    return lst

def exponential(lower,upper,scale):
    return truncexpon(b=(upper-lower)/scale,loc = lower,scale = scale).rvs()

def generate_grid(lst1,lst2):
    lst = norm_parameters(thread_pitch[0],thread_pitch[1])
    lst1.append(random.uniform(thread_diameter[0],thread_diameter[1]))
    while(lst1[-1]<total_length):
        temp1 = normal(lst[0],lst[1],lst[2],lst[3])
        lst1.append(lst1[-1]+temp1)
        temp2 = random.uniform(thread_diameter[0],min_or_same(temp1,thread_diameter[1]))
        lst2.append(temp1-temp2)

def generate_grid_1(lst1,lst2):
    lst1.append(random.uniform(thread_diameter[0],thread_diameter[1]))
    while(lst1[-1]<total_length):
        temp1 = exponential(thread_pitch[0],thread_pitch[1],10)
        pitch.append(temp1)
        lst1.append(lst1[-1]+temp1)
        temp2 = exponential(thread_diameter[0],min_or_same(temp1,thread_diameter[1]),15)*(-1)+thread_diameter[0]+thread_diameter[1]
        lst2.append(temp1-temp2)
        diameter.append(temp2)


def porosity(lst1,lst2):
    temp = sum(lst2)
    free_area = 0
    for i in lst1:
        free_area = free_area + (i*temp)
    return (free_area/(total_length**2))*100

total_length = float(input("Enter the total length of the mask (in cm) \n"))*1e4
thread_diameter = [float(ele) for ele in input("Enter thread diameter (in um)\n").split("-")]
thread_pitch = [float(ele) for ele in input("Enter thread pitch (in um)\n").split("-")]

x_coordinates = []
y_coordinates = []
x_pitch = []
y_pitch = []

x_coordinates1 = []
y_coordinates1 = []
x_pitch1 = []
y_pitch1 = []

pitch = []
diameter = []

generate_grid_1(x_coordinates,x_pitch)
generate_grid_1(y_coordinates,y_pitch)

#generate_grid_1(x_coordinates1,x_pitch1)
#generate_grid_1(y_coordinates1,y_pitch1)

#print(pitch)
print("The porosity of the mask is ",porosity(x_pitch,y_pitch))

plt.hist(pitch)
plt.show()
#print("The exponential porosity of the mask is ",porosity(x_pitch1,y_pitch1))
