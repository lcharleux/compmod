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
   
