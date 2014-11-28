# Distributions

from scipy import stats
import numpy as np

def Triangular(mean = 1., stdev = 1.):
   """
   A triangular symetric distribution
   
   ...
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
   
   
