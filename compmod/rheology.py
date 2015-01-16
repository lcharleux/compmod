# -*- coding: utf-8 -*-
# Rheology 

import numpy as np

class SaintVenant(object):
  """
  A class for parallel assembly of unit cells.
  
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
  sigma = np.where( epsilon < epsilon_y, E * epsilon, sigmay + n * (epsilon - epsilon_y))
  if sigma_sat != None:
    sigma = np.where( sigma < sigma_sat, sigma, sigma_sat)
  return sigma
      