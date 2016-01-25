from abapy import materials
from abapy.misc import load
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pickle, copy, platform, compmod


# GENERAL SETTINGS
settings = {}
settings['export_fields'] = True
settings['compart'] = True
settings['is_3D']   = True
settings['lx']      = 1.
settings['ly']      = 1. # test direction 
settings['lz']      = 1.
settings['Nx']      = 10
settings['Ny']      = 10
settings['Nz']      = 10
settings['disp']    = .01
settings['nFrames'] = 50
settings['workdir'] = "workdir/"
settings['label']   = "cuboidTest_3D"
settings['elType']  = "CPS4"
settings['cpus']    = 1

run_simulation = True # True to run it, False to just use existing results

E           = 1.
sy          = 0.001
nu          = 0.3
ratio       = .7

# ABAQUS PATH SETTINGS
node = platform.node()
if node ==  'lcharleux':      
  settings['abqlauncher']   = "/opt/Abaqus/6.9/Commands/abaqus" 
if node ==  'serv2-ms-symme': 
  settings['abqlauncher']   = "/opt/abaqus/Commands/abaqus"
if node ==  'epua-pd47': 
  settings['abqlauncher']   = "C:/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe" 
if node ==  'epua-pd45': 
  settings['abqlauncher']   = "C:\SIMULIA/Abaqus/Commands/abaqus"
if node ==  'SERV3-MS-SYMME': 
  settings['abqlauncher']   = "C:/Program Files (x86)/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe"


# MATERIALS CREATION
Ne = settings['Nx'] * settings['Ny'] 
if settings['is_3D']: Ne *= settings['Nz']
labels = ['mat_{0}'.format(i+1) for i in xrange(Ne)]
settings['material'] = []
for i in xrange(Ne):
  if i < ratio * Ne: 
    settings['material'].append(materials.VonMises(labels = labels[i], 
                                                   E = E, nu = nu, 
                                                   sy = sy))
  else:
    settings['material'].append(materials.Elastic(labels = labels[i], 
                                                 E = E, nu = nu )) 
np.random.shuffle(settings['material'])      
m = compmod.models.CuboidTest(**settings)
if run_simulation:
  m.MakeInp()
  m.Run()
  m.MakePostProc()
  m.RunPostProc()
m.LoadResults()
# Plotting results
if m.outputs['completed']:
  

  # History Outputs
  disp =  np.array(m.outputs['history']['disp'].values()[0].data[0])
  force =  np.array(np.array(m.outputs['history']['force'].values()).sum().data[0])
  volume = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[0])
  length = settings['ly'] + disp
  surface = volume / length
  logstrain = np.log10(1. + disp / settings['ly'])
  linstrain = disp/ settings['ly']
  strain = linstrain
  stress = force / surface 
   
  fig = plt.figure(0)
  plt.clf()
  sp1 = fig.add_subplot(2, 1, 1)
  plt.plot(disp, force, 'ro-')
  plt.xlabel('Displacement, $U$')
  plt.ylabel('Force, $F$')
  plt.grid()
  sp1 = fig.add_subplot(2, 1, 2)
  plt.plot(strain, stress, 'ro-', label = 'simulation curve', linewidth = 2.)
  plt.xlabel('Tensile Strain, $\epsilon$')
  plt.ylabel(' Tensile Stress $\sigma$')
  plt.grid()
  plt.savefig(settings['workdir'] + settings['label'] + 'history.pdf')
  
  
  # Field Outputs
  if settings["export_fields"]:
    m.mesh.dump2vtk(settings['workdir'] + settings['label'] + '.vtk')
  
