from compmod.models import CuboidTest
from abapy import materials
from abapy.misc import load
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform


#PAREMETERS
lx, ly = 1., 1.
Nx, Ny = 4, 4 
disp = .1
nFrames = 100
workdir = "workdir/"
label = "cuboidTest"
elType = "CPS4"
node = platform.node()
if node ==  'lcharleux':      abqlauncher   = '/opt/Abaqus/6.9/Commands/abaqus' # Local machine configuration
if node ==  'serv2-ms-symme': abqlauncher   = '/opt/abaqus/Commands/abaqus' # Local machine configuration


Ne = Nx * Ny
E  = 1. * np.ones(Ne) # Young's modulus
nu = .3 * np.ones(Ne) # Poisson's ratio
sy_mean = .01
sy = np.random.rayleigh(sy_mean, Ne)
labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
material = [materials.VonMises(labels = labels[i], E = E[i], nu = nu[i], sy = sy[i]) for i in xrange(Ne)]

m = CuboidTest(lx =lx, ly = ly, Nx = Nx, Ny = Ny, abqlauncher = abqlauncher, label = label, workdir = workdir, material = material, compart = True)
m.MakeInp()
