# Distributions

from scipy import stats
import numpy as np

def Triangular(mean = 1., stdev = 1.):
  """
  A triangular symetric distribution function that returns a frozen distribution of the `scipy.stats.rv_continuous <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.rv_continuous.html>`_ class.
   
  :param mean: mean value
  :type mean: float
  :param stdev: standard deviation
  :type stdev: float
  :rtype: scipy.stats.rv_continuous instance
   
   
  >>> import compmod
  >>> tri = compmod.distributions.Triangular
  >>> tri = compmod.distributions.Triangular(mean = 1., stdev = .1)
  >>> tri.rvs(10)
  array([ 1.00410636,  1.05395898,  1.03192428,  1.01753651,  0.99951611,
          1.1718781 ,  0.94457269,  1.11406294,  1.08477038,  0.98861803])


  .. plot:: example_code/distributions/triangular.py
     :include-source:
  """
  width = np.sqrt(6) * stdev 
  return stats.triang(.5, loc = mean - width, scale = 2.*width)


def Rectangular(mean = 1., stdev = 1.):
   """
   A Rectangular symetric distribution
   
   ...
   """
   width = np.sqrt(3) * stdev 
   return stats.uniform(loc = mean - width, scale = 2.*width)
   
   
