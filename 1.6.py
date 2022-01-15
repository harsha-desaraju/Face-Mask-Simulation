from scipy.stats import *
import matplotlib.pyplot as plt

r = truncexpon(b=(230-190)/1,loc=190,scale=1)
data = r.rvs(20000)*1e-6
print(data)
#print(r)
plt.hist(data)
plt.show()
