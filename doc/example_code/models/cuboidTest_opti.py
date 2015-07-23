# SOME OPTIMIZATION WITH RING COMPRESSION

from abapy import materials
from compmod.models import CuboidTest_VER
from scipy import interpolate
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform
node = platform.node()

def read_file(file_name):
  '''
  Read a two rows data file and converts it to numbers
  '''
  f = open(file_name, 'r') # Opening the file
  lignes = f.readlines() # Reads all lines one by one and stores them in a list
  f.close() # Closing the file
#    lignes.pop(0) # Delete le saut de ligne for each lines
  strain_exp, stress_exp = [],[]

  for ligne in lignes:
      data = ligne.split() # Lines are splitted
      strain_exp.append(float(data[0]))
      stress_exp.append(float(data[1]))
  return np.array(strain_exp), np.array(stress_exp)


lateralbc = {"top":"pseudohomo"} # lateral boundary conditions : "pseudohomo"--> lateral nodes have the same displacement
is_3D = True
export_fields = False
label = "CuboidTestOpti"
cpus = 1
compart = False
if is_3D == False :
  elType = "CPS4"
else:
  elType = "C3D8"
#FIXED PAREMETERS
settings = {}

settings['file_name'] = 'Courbe_ref_alu1.txt'
strain_exp, stress_exp = read_file(settings['file_name'])


settings['lx'], settings['ly'], settings['lz']  = 1., 2., 1. #ly = tensile test direction
settings['Nx'], settings['Ny'], settings['Nz'] = 10, 20, 10
if is_3D == True :
    settings['Ne'] =  settings['Nx']*settings['Ny']*settings['Nz']
else :
    settings['Ne'] =  settings['Nx']*settings['Ny']
settings['displacement'] = strain_exp[-1]*settings['ly']
settings['nFrames'] = 100
settings['E'] = 64000.
settings['nu'] = .3 
settings['iteration'] = 20
#settings['thickness'] = 20.02


if node ==  'lcharleux':      
  abqlauncher   = '/opt/Abaqus/6.9/Commands/abaqus' # Local machine configuration
  workdir = "workdir/"
if node ==  'serv2-ms-symme': abqlauncher   = '/opt/abaqus/Commands/abaqus' # Linux
if node ==  'epua-pd47': 
  abqlauncher   = 'C:/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe' # Local machine configuration
  workdir = "workdir/"
if node ==  'SERV3-MS-SYMME': 
  abqlauncher   = '"C:/Program Files (x86)/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe"' # Local machine configuration
  workdir = "workdir/"
if node ==  'epua-pd45': 
  abqlauncher   = 'C:\SIMULIA/Abaqus/Commands/abaqus'  


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
    lx = self.settings['lx']
    ly = self.settings['ly']
    lz = self.settings['lz']
    disp = self.settings['displacement']
    nFrames = self.settings['nFrames']
    Nx = self.settings['Nx']
    Ny = self.settings['Ny']
    Nz = self.settings['Nz']
    Ne = self.settings['Ne']
    #thickness = self.settings['thickness']
    
    #TASKS
    run_sim = True
    plot = True

    print E, nu, sy, n

    labels = 'SAMPLE_MAT'
    material = materials.Hollomon(labels = labels, E = E, nu = nu, sy = sy, n=n)

    m = CuboidTest_VER(lx =lx, ly = ly, lz = lz, Nx = Nx, Ny = Ny, Nz = Nz, abqlauncher = abqlauncher, label = label, workdir = workdir, material = material, compart = compart, disp = disp, elType = elType, is_3D = is_3D, lateralbc = lateralbc, export_fields = export_fields, cpus = cpus)
    
    # SIMULATION
    m.MakeMesh()
    if run_sim:
      m.MakeInp()
      m.Run()
      m.PostProc()
      if m.outputs['completed']:
          # History Outputs
          disp =  np.array(m.outputs['history']['disp'].values()[0].data[0])
          force =  np.array(np.array(m.outputs['history']['force'].values()).sum().data[0])
          volume = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[0])
          length = ly + disp
          surface = volume / length
          logstrain = np.log10(1. + disp / ly)
          linstrain = disp/ly
          strain = linstrain
          stress = force / surface 
    
      self.disp = disp
      self.force = force
      self.logstrain = logstrain
      self.stress = stress
      self.strain = strain
    
  
  
  def Interp(self):
    """
    Interpolate the curve stress_vs_strain on a known grid
    """
    strain, stress = self.strain, self.stress
    f = interpolate.interp1d(strain, stress)
    return f

class Opti(object):
  
  def __init__(self, sy0, n0, settings):
    
    self.sy0 = sy0
    self.n0 = n0
    self.settings = settings
    self.sy = []
    self.n = []
    self.err = []
    self.stress_sim = []
    strain_exp, stress_exp = read_file(self.settings['file_name'])
    g = interpolate.interp1d(strain_exp, stress_exp)
    self.strain_exp = strain_exp
    self.stress_exp = stress_exp
    self.g = g

  def Err(self, param):
    """
    Compute the residual error between experimental and simulated curve
    """    
    n =param[1]
    sy = param[0]

    s = Simulation(sy, n , self.settings)

    s.Run()
    f = s.Interp()
    d = self.settings['displacement']*0.999
    ly = self.settings['ly']
    strain_grid = np.linspace(0., d, 100)/ly
    stress_sim = f(strain_grid)
    
    g = self.g
    stress_exp = g(strain_grid)
    
    err = np.sqrt(((stress_exp - stress_sim)**2).sum())
    self.sy.append(sy)
    self.n.append(n)
    self.err.append(err)
    self.stress_sim.append(stress_sim)
    self.stress_exp = stress_exp
    self.strain_grid = strain_grid
    
    return err
    
  def Optimize(self):
    p0 = [self.sy0, self.n0]   
    result = minimize(self.Err, p0, method='nelder-mead', options={'disp':True, 'maxiter':settings['iteration']})
    self.result = result

   
O = Opti(100., 0.1,  settings)

O.Optimize()


fig = plt.figure('Stress vs. strain')
plt.clf()

plt.plot(O.strain_grid, O.stress_exp, 'k-', label = 'experimental curve', linewidth = 2.)
plt.plot(O.strain_grid, O.stress_sim[0], 'g-', label = 'initial curve', linewidth = 2.)
a = O.err
index = np.argmin(a)
plt.plot(O.strain_grid, O.stress_sim[index], 'r-', label = 'optimized curve', linewidth = 2.)
for i in range(1, settings['iteration']):
  plt.plot(O.strain_grid, O.stress_sim[i], 'b-', linewidth = .2)
plt.legend(loc="lower right")
plt.grid()
plt.xlabel('Strain, $\epsilon$ (%)')
plt.ylabel('Stress, $\sigma$ (Mpa)')
plt.savefig(workdir + label + '_stress-vs-strain.pdf')


#print s.force.data[0]
"""
f = s.Interp()
x = np.arange(0., 49., 0.1)
print f(x)
"""




