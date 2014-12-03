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
   A Rectangular symetric distribution function that returns a frozen uniforn distribution of the 'scipy.stats.rv_continuous <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.uniform.html>'_class.
   
   :param mean: mean value
   :type mean: float
   :param stdev: standard deviation
   :type stdev: float
   :rtype: scipy.stats.rv_continuous instance


  >>> import compmod 
  >>> rec = compmod.distributions.Rectangular
  >>> rec = compmod.distributions.Rectangular(mean = 5. , stdev = 2.) 
  >>> rec.rvs(15)
  array([ 6.30703805,  5.55772119,  5.69890282,  5.41807602,  6.78339394,
        1.83640732,  3.50697054,  7.97707174,  4.54666157,  7.27897515,
        2.33288284,  2.62291176,  1.80274279,  3.39480096,  6.09699301])

  .. plot:: example_code/distributions/rectangular.py
     :include-source:
   """
   width = np.sqrt(3) * stdev 
   return stats.uniform(loc = mean - width, scale = 2.*width)
   
def Rayleigh(mean = 1):
   """
   A Rayleigh distribution
   
   ...
   """
   return stats.rayleigh(mean-1)
      
