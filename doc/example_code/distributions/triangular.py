import numpy as np
import matplotlib.pyplot as plt
from compmod.distributions import Triangular

N = 1000
mean, stdev = 5., 2.
tri = Triangular(mean = mean, stdev = stdev)
data = tri.rvs(N)

x =np.linspace(0., 10., 1000)
y = tri.pdf(x)
plt.figure()
plt.clf()
plt.hist(data, bins = int(N**.5), histtype='step', normed = True, label = "Generated Random Numbers")

plt.plot(x,y, "r-", label = "Probability Density Function")
plt.grid()
plt.legend(loc = "best")
plt.show()
