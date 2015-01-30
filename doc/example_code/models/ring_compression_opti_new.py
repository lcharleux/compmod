# SOME OPTIMIZATION WITH RING COMPRESSION

from abapy.materials import Hollomon
from compmod.models import RingCompression
from scipy import interpolate
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform



#FIXED PAREMETERS
settings = {}
settings['file_name'] = 'test_expD2.txt'
settings['inner_radius'], settings['outer_radius'] = 45.2 , 48.26
settings['Nt'], settings['Nr'], settings['Na'] = 20, 4, 5
settings['Ne'] =  settings['Nt']*settings['Nr']*settings['Na']
settings['displacement'] = 45.
settings['nFrames'] = 100
settings['E'] = 71413.
settings['nu'] = .3
settings['iteration'] = 15
settings['thickness'] = 50



workdir = "workdir/"
label = "ringCompression_opti"
elType = "CPE4"
cpus = 1
node = platform.node()
if node ==  'lcharleux':      abqlauncher   = '/opt/Abaqus/6.9/Commands/abaqus' # Ludovic
if node ==  'serv2-ms-symme': abqlauncher   = '/opt/abaqus/Commands/abaqus' # Linux
if node ==  'epua-pd47': 
  abqlauncher   = 'C:/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe' # Local machine configuration
if node ==  'SERV3-MS-SYMME': 
  abqlauncher   = '"C:/Program Files (x86)/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe"' # Local machine configuration
if node ==  'epua-pd45': 
  abqlauncher   = 'C:\SIMULIA/Abaqus/Commands/abaqus' 
  

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
      disp_exp.append(float(data[0]))
      force_exp.append(float(data[1]))
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
    Na = self.settings['Na']
    Ne = self.settings['Ne']
    thickness = self.settings['thickness']
    print E, nu, sy, n
    
    material = Hollomon(
      labels = "SAMPLE_MAT",
      E = E, nu = nu,
      sy = sy, n = n)
    m = RingCompression( material = material , 
      inner_radius = inner_radius, 
      outer_radius = outer_radius, 
      disp = disp,
      thickness = thickness,
      nFrames = nFrames, 
      Nr = Nr, 
      Nt = Nt, 
      Na = Na,
      workdir = workdir,
      label = label, 
      elType = elType,
      abqlauncher = abqlauncher,
      cpus = cpus,
      is_3D = False)
  
    # SIMULATION
    m.MakeMesh()
    m.MakeInp()
    m.Run()
    m.PostProc()
    outputs = m.outputs
    force = -2. * outputs['history']['force']
    disp = -2. * outputs['history']['disp']
    
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
    self.force_sim = []
    disp_exp, force_exp = read_file(self.settings['file_name'])
    g = interpolate.interp1d(disp_exp, force_exp)
    self.disp_exp = disp_exp
    self.force_exp = force_exp
    d = self.settings['displacement']
    self.disp_grid = np.linspace(0., d, 1000)
    self.force_exp_grid= g(self.disp_grid)

  def Err(self, param):
    """
    Compute the residual error between experimental and simulated curve
    """
    sy = param[0]
    n =param[1]
    disp_grid = self.disp_grid
    s = Simulation(sy, n ,self.settings)
    s.Run()
    f = s.Interp()
    force_sim = f(disp_grid)
    force_exp_grid = self.force_exp_grid
    
    err = np.sqrt(((force_exp_grid - force_sim)**2).sum())
    self.sy.append(sy)
    self.n.append(n)
    self.err.append(err)
    self.force_sim.append(force_sim)
    return err
    
  def Optimize(self):
    p0 = [self.sy0, self.n0]
    result = minimize(self.Err, p0, method='nelder-mead', options={'disp':True, 'maxiter':settings['iteration']})
    self.result = result
    
O = Opti(160., 0.1, settings)
O.Optimize()


fig = plt.figure('Load vs. disp')
plt.clf()

plt.plot(O.disp_grid, O.force_exp_grid, 'k-', label = 'experimental curve', linewidth = 2.)
plt.plot(O.disp_grid, O.force_sim[0], 'g-', label = 'initial curve', linewidth = 2.)
a = O.err
index = np.argmin(a)
plt.plot(O.disp_grid, O.force_sim[index], 'r-', label = 'optimized curve', linewidth = 2.)
for i in range(1, settings['iteration']):
  plt.plot(O.disp_grid, O.force_sim[i], 'b-', linewidth = .2)
#plt.plot(disp.data[1], force.data[1], 'b-', label = 'Unloading', linewidth = 2.)  
plt.legend(loc="lower right")
plt.grid()
plt.xlabel('Displacement, $U$')
plt.ylabel('Force, $F$')
plt.savefig(workdir + label + '_load-vs-disp.pdf')






