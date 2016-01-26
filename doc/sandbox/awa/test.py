"""
Compartimented optimization

LC 13/05/2015 
"""
import numpy as np
import matplotlib.pyplot as plt
import compmod, abapy, platform
from scipy import optimize, interpolate

eps = np.linspace(.05, 0.5, 100)
sigma = compmod.rheology.Bilinear(eps, E = 1.,sigmay = .005, 
                                                       n = .01, 
                                                       sigma_sat = 0.03)
fig = plt.figure(0)
plt.clf()
plt.plot(eps, sigma)
plt.grid()
plt.show()  
print sigma[:5]
