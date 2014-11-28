# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 13:22:39 2014

@author: issackm
"""
import numpy as np
import matplotlib.pyplot as plt
from compmod.distributions import Rectangular

N = 1000
mean, stdev = 5., 2.
rec = Rectangular(mean = mean, stdev = stdev)
data = rec.rvs(N)

x = np.linspace(0., 10., 1000)
y = rec.pdf(x)
plt.figure()
plt.clf()
plt.hist(data, bins = int(N**.5), histtype='step', normed = True, label = "Generated Random Numbers")

plt.plot(x,y, "g-", label = "Probability Density Function")
plt.grid()
plt.legend(loc = "best")
plt.show()