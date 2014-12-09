# Cuboid test simulation with multiple small simulations

from compmod.models import CuboidTest
from abapy import materials
from abapy.misc import load
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform
from scipy import interpolate



#PARAMETERS
Ns = 4 # Number of simulations to run
lx, ly = 1., 1.
Nx, Ny = 25, 25 
Ne = Nx * Ny
disp = .1
nFrames = 20
workdir = "D:\donnees_pyth/workdir/"
label = "cuboidTest"
elType = "CPE4"
cpus = 1
node = platform.node()
if node ==  'lcharleux':      
  abqlauncher   = '/opt/Abaqus/6.9/Commands/abaqus' # Local machine configuration
if node ==  'serv2-ms-symme': 
  abqlauncher   = '/opt/abaqus/Commands/abaqus' # Local machine configuration
if node ==  'epua-pd47': 
  abqlauncher   = 'C:\SIMULIA\Abaqus\6.13-1'
if node ==  'epua-pd45': 
  abqlauncher   = 'C:\SIMULIA/Abaqus/Commands/abaqus'  
compart = True

Strain = [] # Conteners to store data produced by multiple simulations
Stress = []
strain_grid = np.linspace(0., disp/lx, 100)

for s in xrange(Ns):
  if compart:
    E  = 1. * np.ones(Ne) # Young's modulus
    nu = .3 * np.ones(Ne) # Poisson's ratio
    sy_mean = .01
    sy = np.random.rayleigh(sy_mean, Ne)
    labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
    material = [materials.VonMises(labels = labels[i], E = E[i], nu = nu[i], sy = sy[i]) for i in xrange(Ne)]
  else:
    E = 1.
    nu =.3
    sy = .01
    labels = 'SAMPLE_MAT'
    material = materials.VonMises(labels = labels, E = E, nu = nu, sy = sy)

  m = CuboidTest(lx =lx, ly = ly, Nx = Nx, Ny = Ny, abqlauncher = abqlauncher, label = label, workdir = workdir, cpus = cpus, material = material, compart = compart, disp = disp, elType = elType)

  m.MakeInp()
  m.Run()
  m.MakePostProc()
  m.RunPostProc()


  # Plotting results
  if m.outputs['completed']:
    # History Outputs
    displacement =  np.array(m.outputs['history']['disp'].values()[0].data[0])
    force =  np.array(np.array(m.outputs['history']['force'].values()).sum().data[0])
    volume = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[0])
    length = ly + displacement
    surface = volume / length
    logstrain = np.log10(1. + displacement / ly)
    linstrain = displacement/ly
    strain = linstrain
    stress = force / surface
    f = interpolate.interp1d(strain, stress)
    Stress.append(f(strain_grid))
        
  
Stress = np.array(Stress)
Strain = np.array(Strain)
Average_Stress = np.average(Stress, axis = 0)
fig = plt.figure(0)
plt.clf()
sp1 = fig.add_subplot(1, 1, 1)
for s in xrange(len(Stress)):
  plt.plot(strain_grid, Stress[s], 'g-')
plt.plot(strain_grid, Average_Stress, 'or-', label = 'Average')  
plt.xlabel('Tensile Strain, $\epsilon$')
plt.ylabel('Tensile Stress $\sigma$')
plt.grid()
plt.legend()
plt.savefig(workdir + label + 'history_multiple.pdf')
  
  
