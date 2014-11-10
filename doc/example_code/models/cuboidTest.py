from compmod.models import CuboidTest
from abapy import materials
from abapy.misc import load
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform



#PAREMETERS
lx, ly = 1., 1.
Nx, Ny = 20, 20 
Ne = Nx * Ny
disp = .1
nFrames = 20
workdir = "workdir/"
label = "cuboidTest"
elType = "CPE4"
node = platform.node()
if node ==  'lcharleux':      abqlauncher   = '/opt/Abaqus/6.9/Commands/abaqus' # Local machine configuration
if node ==  'serv2-ms-symme': abqlauncher   = '/opt/abaqus/Commands/abaqus' # Local machine configuration



E  = 1. * np.ones(Ne) # Young's modulus
nu = .3 * np.ones(Ne) # Poisson's ratio
sy_mean = .01
sy = np.random.rayleigh(sy_mean, Ne)
labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
material = [materials.VonMises(labels = labels[i], E = E[i], nu = nu[i], sy = sy[i]) for i in xrange(Ne)]

m = CuboidTest(lx =lx, ly = ly, Nx = Nx, Ny = Ny, abqlauncher = abqlauncher, label = label, workdir = workdir, material = material, compart = True, disp = disp, elType = elType)

m.MakeInp()
m.Run()
m.MakePostProc()

m.RunPostProc()

# Plotting results
if m.outputs['completed']:
  disp =  np.array(m.outputs['history']['disp'].values()[0].data[0])
  force =  np.array(np.array(m.outputs['history']['force'].values()).sum().data[0])
  volume = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[0])
  length = ly + disp
  surface = volume / length
  logstrain = np.log10(1. + disp / ly)
  linstrain = disp/ly
  stress = force / surface


  fig = plt.figure(0)
  plt.clf()
  sp1 = fig.add_subplot(2, 1, 1)
  plt.plot(disp, force, 'ok-')
  plt.xlabel('Displacement, $U$')
  plt.ylabel('Force, $F$')
  plt.grid()
  sp1 = fig.add_subplot(2, 1, 2)
  plt.plot(linstrain, stress, 'ok-')
  plt.xlabel('Tensile Strain, $\epsilon$')
  plt.ylabel(' Tensile Stress $\sigma$')
  plt.grid()
  plt.savefig(label + 'history.pdf')




