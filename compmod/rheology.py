"""
Rheology
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize as opt

class Cell(object):
  """
  A cell meta class.
  """
  def __init__(self, Nhist, weight):
    sk = self._state_keys
    self.state = pd.DataFrame(np.zeros((Nhist, len(sk))) , columns = sk)
    self.weight = weight
      
  def step_forward(self):
    """
    Makes a step forward in the history of the state variables
    """
    for row in self.state.iloc[:-1].iloc[::-1].itertuples(): # All but the last in reversed order...
      self.state.loc[row[0]+1] = row[1:]
      
  
  def step_backward(self):
    """
    Makes a step backward in the history of state variables
    """
    for row in self.state.iloc[1:].itertuples(): # All but the first ! 
      self.state.loc[row[0]-1] = row[1:]

class ElasticPlasticCell(Cell):
  """
  An Elastic (perfectly) plastic cell.
  """
  _state_keys = ["epse", "epsp"] # State variables
  
  def __init__(self, E = 1., sy = 1., weight = 1., Nhist = 2):       
    self.E =  E
    self.sy = sy
    Cell.__init__(self, Nhist, weight)
            
  def set_eps(self, eps):
    """
    Sets the value of the total strain eps.
    """
    state = self.state
    epse, epsp = state.epse.iloc[0], state.epsp.iloc[0]
    self.step_forward()
    sy, E = self.sy, self.E
    ey = sy/E
    deps = eps - epsp
    if abs(deps) <= ey:
      # Elastic domain
      state.epse.iloc[0]  = deps
      
    else:
      # Plastic loading
      state.epsp.iloc[0] += deps - np.sign(deps) * ey
      state.epse.iloc[0] = eps - state.epsp.iloc[0]
    #print "eps={0}".format(self.eps)
    
  def get_eps(self):
    return self.state.epse.iloc[0] + self.state.epsp.iloc[0]       
  
  eps = property(get_eps, set_eps)
  
  def sigma(self):
    E = self.E
    epse = self.state.epse.iloc[0]
    return E * epse

  
class CellBox(object):
  def __init__(self, cells):
    self.cells = cells
    
  def set_eps(self, eps):
    """
    Sets the value of the total strain eps.
    """
    for cell in self.cells:
      cell.eps = eps
    
  def get_eps(self):
    return self.cells[0].eps
  
  eps = property(get_eps, set_eps)  
  
  def sigma(self):
    sigma = np.array([cell.sigma() for cell in self.cells])
    weights = np.array([cell.weight for cell in self.cells])
    return (sigma * weights).sum() / weights.sum()
  
  def find_eps(self, sigma, epsa, epsb):
    """
    Finds the eps value corresponding to a given sigma.
    """
    def fonc(eps):
      """
      Target function
      """
      self.eps = eps
      error = self.sigma() - sigma
      self.step_backward()
      return error
    self.eps = opt.brentq(fonc, epsa, epsb)
  
  
    
  def step_backward(self):
    for cell in self.cells:
      cell.step_backward()
      
      
