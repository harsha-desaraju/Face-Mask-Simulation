import random
import matplotlib.pyplot as plt

def throw_ball():
    return [random.uniform(0,total_len),random.uniform(0,total_len)]

def find_grid(ball):
    for i in range(len(grid)):
        if grid[i]<ball[0] and grid[i+1]>ball[0]:
            break
    for j in range(len(grid)):
        if grid[j]<ball[1] and grid[j+1]>ball[1]:
            break
    return [grid[i],grid[j]]

def ball_passes(ball):
    grid_cord = find_grid(ball)
    if ball[0]-grid_cord[0]>ball_rad and ball[1]-grid_cord[1]>ball_rad:
        if (grid_cord[0]+pore_size)-ball[0]>ball_rad and (grid_cord[1]+pore_size)-ball[1]>ball_rad:
            return True
    else:
        return False

total_len = float(input("Enter the side length of the square grid\n"))
pore_size = float(input("Enter the size of square pore\n"))
ball_rad = float(input("Enter the radius of ball\n"))
N = int(input("Enter the number of balls\n"))

grid = []
n = int(total_len/pore_size)
for i in range(n+1):
    grid.append(round(pore_size*i,7))

success=0
prob = []
x = []

for i in range(1,N+1):
    ball = throw_ball()
    if ball_passes(ball):
        success+=1
    prob.append(success/i)
    x.append(i)

if ball_rad<pore_size/2:
    theo_prob = (1-2*(ball_rad/pore_size))*(1-2*(ball_rad/pore_size))
else:
    theo_prob = 0

print("The probability calculated from the simulation is ",prob[-1])
print("The theoretical probability is ",theo_prob)

plt.plot(x,prob,label = "simulated probability")
plt.plot([0,N],[theo_prob,theo_prob],label = "Theoretical probability")
plt.legend()
plt.xlabel("No. of balls")
plt.ylabel("probability")

plt.show()
