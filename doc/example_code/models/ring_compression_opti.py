# SOME OPTIMIZATION WITH RING COMPRESSION

from abapy.materials import Hollomon
from compmod.models import RingCompression
from scipy import interpolate
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform
node = platform.node()


#FIXED PAREMETERS
settings = {}
settings['file_name'] = 'test_exp.txt'
settings['inner_radius'], settings['outer_radius'] = 95/2, 100/2
settings['Nt'], settings['Nr'] = 10, 3 
settings['displacement'] = 45.
settings['nFrames'] = 100
settings['E'] = 200.e9
settings['nu'] = .3

if node ==  'lcharleux':      
  abqlauncher   = '/opt/Abaqus/6.9/Commands/abaqus' # Local machine configuration
  workdir = "workdir/"
if node ==  'epua-pd47': 
  abqlauncher   = 'C:/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe' # Local machine configuration
  workdir = "D:/Simulations/Dossier_travail_Abaqus/"
label = "ringCompressionOpti"
elType = "CPS4"


def read_file(file_name):
  '''
  Read a two rows data file and converts it to numbers
  '''
  f = open(file_name, 'r') # Opening the file
  lignes = f.readlines() # Reads all lines one by one and stores them in a list
  f.close() # Closing the file
#    lignes.pop(0) # Delete le saut de ligne for each lines
  force_exp, disp_exp = [],[]

  for ligne in lignes:
      data = ligne.split() # Lines are splitted
      force_exp.append(float(data[0]))
      disp_exp.append(float(data[1]))
  return -np.array(disp_exp), -np.array(force_exp)



class Simulation(object):
  
  def __init__(self, sy, n, settings):
    self.sy = sy
    self.n = n
    self.settings = settings
    
    
  def Run(self):
    """
    Runs a simulation for a given couple (sy, n) and returns the (disp, force) couple.
    """
    #MODEL DEFINITION
    sy = self.sy
    n = self.n
  
    E = self.settings['E']
    nu = self.settings['nu']
    inner_radius = self.settings['inner_radius']
    outer_radius = self.settings['outer_radius']
    disp = self.settings['displacement']/2.
    nFrames = self.settings['nFrames']
    Nr = self.settings['Nr']
    Nt = self.settings['Nt']
    print E, nu, sy, n
    material = Hollomon(
      labels = "SAMPLE_MAT",
      E = E, nu = nu,
      sy = sy, n = n)
    m = RingCompression(
      material = material , 
      inner_radius = inner_radius, 
      outer_radius = outer_radius, 
      disp = disp, 
      nFrames = nFrames, 
      Nr = Nr, 
      Nt = Nt, 
      workdir = workdir,
      label = label, 
      elType = elType,
      abqlauncher = abqlauncher)
  
    # SIMULATION
    m.MakeMesh()
    m.MakeInp()
    m.Run()
    m.PostProc()
    outputs = m.outputs
    force = -2. * outputs['history']['force']
    disp = - 2 * outputs['history']['disp']
    
    self.disp = disp
    self.force = force
    
  
  
  def Interp(self):
    """
    Interpolate the curve Force-displacement on a known grid
    """
    disp, force = self.disp, self.force
    f = interpolate.interp1d(disp.data[0], force.data[0])
    return f

class Opti(object):
  
  def __init__(self, sy0, n0, settings):
    
    self.sy0 = sy0
    self.n0 = n0
    self.settings = settings
    self.sy = []
    self.n = []
    self.err = []
    disp_exp, force_exp = read_file(self.settings['file_name'])
    g = interpolate.interp1d(disp_exp, force_exp)
    self.disp_exp = disp_exp
    self.force_exp = force_exp
    self.g = g
    
  def Err(self, param):
    """
    Compute the residual error between experimental and simulated curve
    """
    sy = param[0]
    n =param[1]
   
    s = Simulation(sy, n ,self.settings)
    s.Run()
    f = s.Interp()
    d = self.settings['displacement']
    disp = np.linspace(0., d, 100)
    force_sim = f(disp)
    
    g = self.g
    force_exp = g(disp)
    
    err = np.sqrt(((force_exp - force_sim)**2).sum())
    self.sy.append(sy)
    self.n.append(n)
    self.err.append(err)
    return err
    
  def Optimize(self):
    p0 = [self.sy0, self.n0]
    
    result = minimize(self.Err, p0, method='nelder-mead', options={'disp':True, 'maxiter':5})
    self.result = result
    
O = Opti(200.e6, 0.3, settings)
O.Optimize()



#print s.force.data[0]
"""
f = s.Interp()
x = np.arange(0., 49., 0.1)
print f(x)

"""





