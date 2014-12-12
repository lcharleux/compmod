# SOME OPTIMIZATION WITH RING COMPRESSION

from abapy import materials
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
settings['inner_radius'], settings['outer_radius'] = 45.18 , 50.36
settings['Nt'], settings['Nr'], settings['Na'] = 40, 5, 10
settings['Ne'] =  settings['Nt']*settings['Nr']*settings['Na']
settings['displacement'] = 35.
settings['nFrames'] = 100
settings['E'] = 74.e3 * np.ones(settings['Ne'])
settings['nu'] = .3 * np.ones(settings['Ne'])
settings['iteration'] =10
settings['thickness'] = 20.02


if node ==  'lcharleux':      abqlauncher   = '/opt/Abaqus/6.9/Commands/abaqus' # Ludovic
if node ==  'serv2-ms-symme': abqlauncher   = '/opt/abaqus/Commands/abaqus' # Linux
if node ==  'epua-pd47': 
  abqlauncher   = 'C:/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe' # Local machine configuration
if node ==  'SERV3-MS-SYMME': 
  abqlauncher   = '"C:/Program Files (x86)/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe"' # Local machine configuration
if node ==  'epua-pd45': 
  abqlauncher   = 'C:\SIMULIA/Abaqus/Commands/abaqus'  
workdir = "D:\donnees_pyth/workdir/"
label = "ringCompressionOptiCompart"
elType = "CPS4"
cpus = 1




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
  
  def __init__(self, Ssat, n, sy_mean, settings):
    self.sy_mean = sy_mean
    self.n = n
    self.Ssat = Ssat
    self.settings = settings
    
    
  def Run(self):
    """
    Runs a simulation for a given couple (sy, n) and returns the (disp, force) couple.
    """
    #MODEL DEFINITION
    sy_mean = self.sy_mean * np.ones(settings['Ne'])
    n = self.n * np.ones(settings['Ne'])
    Ssat = self.Ssat * np.ones(settings['Ne'])
  
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
    
    #TASKS
    run_sim = True
    plot = True
    
    print E[0], nu[0], Ssat[0], n[0], sy_mean[0]
  
    ray_param = sy_mean/1.253314
    sy = np.random.rayleigh(ray_param, Ne)
    labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
    material = [materials.Bilinear(labels = labels[i], E = E[i], nu = nu[i], Ssat = Ssat[i], n=n[i], sy = sy[i]) for i in xrange(Ne)]

    m = RingCompression(
      material = material , 
      inner_radius = inner_radius, 
      outer_radius = outer_radius, 
      disp = disp, 
      nFrames = nFrames,
      thickness = thickness,
      Nr = Nr, 
      Nt = Nt, 
      workdir = workdir,
      label = label, 
      elType = elType,
      abqlauncher = abqlauncher,
      cpus = cpus,
      is_3D = True,
      compart = True)
  
    # SIMULATION
    m.MakeMesh()
    if run_sim:
      m.MakeInp()
      m.Run()
      m.PostProc()
      outputs = m.outputs
      force = -2. * outputs['history']['force']
      disp = -2 * outputs['history']['disp']
    
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
  
  def __init__(self, Ssat0, n0, sy_mean0, settings):
    
    self.sy_mean0 = sy_mean0
    self.n0 = n0
    self.Ssat0 = Ssat0
    self.settings = settings
    self.sy_mean = []
    self.n = []
    self.Ssat = []
    self.err = []
    self.force_sim = []
    disp_exp, force_exp = read_file(self.settings['file_name'])
    g = interpolate.interp1d(disp_exp, force_exp)
    self.disp_exp = disp_exp
    self.force_exp = force_exp
    self.g = g

  def Err(self, param):
    """
    Compute the residual error between experimental and simulated curve
    """
    n =param[1]
    Ssat = param[0]
    sy_mean = param[2]
   
    s = Simulation(Ssat, n , sy_mean, self.settings)
    s.Run()
    f = s.Interp()
    d = self.settings['displacement']
    disp = np.linspace(0., d, 100)
    force_sim = f(disp)
    
    g = self.g
    force_exp = g(disp)
    
    err = np.sqrt(((force_exp - force_sim)**2).sum())
    self.sy_mean.append(sy_mean)
    self.n.append(n)
    self.Ssat.append(Ssat)
    self.err.append(err)
    self.force_sim.append(force_sim)
    self.disp = disp
    self.force_exp = force_exp
    
    
    return err
    
  def Optimize(self):
    p0 = [self.Ssat0, self.n0, self.sy_mean0]
    
    result = minimize(self.Err, p0, method='nelder-mead', options={'disp':True, 'maxiter':settings['iteration']})
    self.result = result
    
O = Opti(800., 200., 200., settings)
O.Optimize()


fig = plt.figure('Load vs. disp')
plt.clf()

plt.plot(O.disp, O.force_exp, 'k-', label = 'experimental curve', linewidth = 2.)
plt.plot(O.disp, O.force_sim[0], 'g-', label = 'initial curve', linewidth = 2.)
a = O.err
index = np.argmin(a)
plt.plot(O.disp, O.force_sim[index], 'r-', label = 'optimized curve', linewidth = 2.)
for i in range(1, settings['iteration']):
  plt.plot(O.disp, O.force_sim[i], 'b-', linewidth = .2)
#plt.plot(disp.data[1], force.data[1], 'b-', label = 'Unloading', linewidth = 2.)  
plt.legend(loc="upper left")
plt.grid()
plt.xlabel('Displacement, $U$')
plt.ylabel('Force, $F$')
plt.savefig(workdir + label + '_load-vs-disp.pdf')


#print s.force.data[0]
"""
f = s.Interp()
x = np.arange(0., 49., 0.1)
print f(x)
"""




