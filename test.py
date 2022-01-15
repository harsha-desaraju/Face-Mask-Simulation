
prob = []
l = []

for i in range(10):
    prob.append([])
    for j in range(10):
        prob[i].append(0)


for i in prob:
    print(i)

for i in range(10):
    for j in range(10):
        prob[i][j] = i+j

for i in prob:
    print(i)
