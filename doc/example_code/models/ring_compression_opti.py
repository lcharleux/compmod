# SOME OPTIMIZATION WITH RING COMPRESSION

from abapy.materials import Hollomon
from compmod.models import RingCompression
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform
node = platform.node()


#FIXED PAREMETERS
file_name = 'D:/Desktop/Article_compression_anneaux/Optimisation_compression_anneau/test_exp.txt'
inner_radius, outer_radius = 95/2, 100/2
Nt, Nr = 10, 3 
displacement = 50.
nFrames = 50
sy = 200.e6
E = 200.e9
nu = .3
n = .3
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
  return np.array(disp_exp), np.array(force_exp)



class Simulation(object):
  
  def __init__(self, sy, n):
    self.sy = sy
    self.n = n
    
    
  def Run(self):
    """
    Runs a simulation for a given couple (sy, n) and returns the (disp, force) couple.
    """
    #MODEL DEFINITION
    sy = self.sy
    n = self.n    
    material = Hollomon(
      labels = "SAMPLE_MAT",
      E = E, nu = nu,
      sy = sy, n = n)
    m = RingCompression(
      material = material , 
      inner_radius = inner_radius, 
      outer_radius = outer_radius, 
      disp = displacement/2, 
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



s = Simulation(200.e6, 0.3)
s.Run()
s.Interp()
#print s.force.data[0]
"""
f = s.Interp()
x = np.arange(0., 49., 0.1)
print f(x)

"""





