# -*- coding: utf-8 -*-
# Rheology 

import numpy as np

class SaintVenant(object):
  """
  A class for parallel assembly of unit cells exhibiting time indepedent behavior. One parameter can be distributed using the cell fu
  
  .. plot:: example_code/rheology/demo.py
     :include-source:
  """
  def __init__(self, epsilon, cell, grid, dist):
    self.epsilon = epsilon
    self.cell = cell
    self.grid = grid
    self.dist = dist
    self.make_Sigma()
    self.make_Dist()
  
  def make_Sigma(self):
    """
    Computes the stress matrix
    """
    cell = self.cell
    eps = self.epsilon
    grid = self.grid
    Ne = len(eps)
    Ng = len(grid)
    Sigma = np.zeros([Ne, Ng])
    for i in xrange(Ng):
      Sigma[:,i] = cell(eps, grid[i])
    self.Sigma = Sigma
      
  def make_Dist(self):
    """
    Computes the probability density vector.
    """
    grid = self.grid
    self.Dist = self.dist.pdf(grid) 
  
  def sigma(self):
    """
    Computes the global stress vector.
    """
    Sigma = self.Sigma
    Dist = self.Dist
    grid = self.grid
    dg = grid[1] - grid[0]
    return (Sigma * Dist).sum(axis = 1) * dg 
    
    
def Bilinear(epsilon, E = 1., sigmay = .01, n = .1, sigma_sat =  None):
  """
  Loi de comportement bilinéaire saturée en traction monotone.
  """
  epsilon_y = sigmay / E
  
  if n == 0.:
    k = 0.
  else:
    k = (n**-1 + E**-1)**-1
  sigma = np.zeros_like(epsilon)    
  deps = epsilon[1:] - epsilon[:-1]
  for i in xrange(1, len(epsilon)):
    de = deps[i-1]
    sigma[i] = sigma[i-1] + E * de
    if abs(sigma[i]) > sigmay:
      de *= 1 - (sigmay - abs(sigma[i-1]))/E
      sigma[i] = np.sign(sigma[i-1]) * sigmay + k * de
      sigmay = abs(sigma[i]) 
  return sigma 
  """
  sigma = np.where( epsilon < epsilon_y, E * epsilon, sigmay + k * (epsilon - epsilon_y))
  if sigma_sat != None:
    sigma = np.where( sigma < sigma_sat, sigma, sigma_sat)
  """
  return sigma



class Bilinear2(object):
  def __init__(self, E = 1., sigma_y = 1., n = .1, sigma_sat = None, N = 1000):       
    self.E =  E
    self. sigma_y = sigma_y
    self.n = n
    self.sigma_sat = sigma_sat
    self._epsilon = np.zeros(N)
    self._sigma = np.zeros(N)
    self._pos = 1
    
  def set_epsilon(self, epsilon):
    pos = self._pos
    try:
      epsilon = float(epsilon)
      self._epsilon[pos] = epsilon
      pos += 1
    except:
      epsilon = np.array(epsilon)
      self._epsilon[pos:pos+len(epsilon)] = epsilon
  def get_epsilon(self):
    pass      
  
  def sigma(self):
    sigma = np.zeros_like(epsilon)
    sigma_y = self.sigma_y
    sigma_sat = self.sigma_sat
    n = self.n
    for e in epsilon:
      pass
      
    
    
    
  
