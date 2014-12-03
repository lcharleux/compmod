import numpy as np
import matplotlib.pyplot as plt
from compmod.distributions import Rayleigh

N = 1000
mean = 2
ray = Rayleigh(mean)
data = ray.rvs(N)

x =np.linspace(0., 10., 1000)
y = ray.pdf(x)
plt.figure()
plt.clf()
plt.hist(data, bins = int(N**.5), histtype='step', normed = True, label = "Generated Random Numbers")

plt.plot(x,y, "r-", label = "Probability Density Function")
plt.grid()
plt.legend(loc = "best")
plt.show()

