# Voronoi cells for grain boundary generation ?
# A demo in 2D (3D is exaclty the same)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate
from compmod.models import CuboidTest_VER
from abapy import materials
from abapy.misc import load
from abapy.postproc import FieldOutput
import pickle, copy
import platform


# PLATFORM SETTINGS
abqlauncher = None
workdir = "workdir/"
is_3D = True
label = "cuboidTest_3D_Voronoi_Hall-Petch"
if is_3D:
  elType = "C3D8"
else:
  elType = "CPS4"
node = platform.node()
if node ==  'lcharleux':      
  abqlauncher   = '/opt/Abaqus/6.9/Commands/abaqus' # Local machine configuration
if node ==  'serv2-ms-symme': 
  abqlauncher   = '/opt/abaqus/Commands/abaqus' # Local machine configuration
if node ==  'epua-pd47': 
  abqlauncher   = 'C:/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe' # Local machine configuration
if node ==  'epua-pd45': 
  abqlauncher   = 'C:\SIMULIA/Abaqus/Commands/abaqus'
if node ==  'SERV3-MS-SYMME': 
  abqlauncher   = '"C:/Program Files (x86)/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe"' # Local machine configuration
if node == 'serv2-ms-symme':
  cpus = 6
else:
  cpus = 1
  
# SIMULATION SETTINGS
lx, ly, lz = 1., 1., 1.
Nx, Ny, Nz = 5, 5, 5
Ne = Nx * Ny * Nz

loading = {"displacement"} #"loading" : force or displacement
#disp = strain_exp[-1] * ly #if the simulation is comparated to experimental data, activate this line
disp = 0.1
force = 190.
force_fin = force + 20.
nFrames = 30
export_fields = True
compart = True
unloading_reloading = False #for one cycle of loading (F = force), unloding (F=0) and reloading (F = force_fin)
Nseed = 20
sigma_0_hp = 0.001
k_hp = .001

Run_simu = True
#lateralbc = {}
#lateralbc = {"right":"periodic","front":"periodic","top":"pseudohomo"}
lateralbc = {"top":"pseudohomo"}


if compart:
  E  = 64000. * np.ones(Ne) # Young's modulus
  nu = .3 * np.ones(Ne) # Poisson's ratio
  sy_mean = 158.6 * np.ones(Ne)
  Ssat = 506. * np.ones(Ne)
  n = 515 * np.ones(Ne)
  #n = 1. * np.ones(Ne)
  ray_param = sy_mean/1.253314
  sy = np.random.rayleigh(ray_param, Ne)
  labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
  material = [materials.Bilinear(labels = labels[i], E = E[i], nu = nu[i], Ssat = Ssat[i], n=n[i], sy = sy[i]) for i in xrange(Ne)]
else:
  E = 70000.
  nu =.3
  sy = 133.1
  n = .081
  labels = 'SAMPLE_MAT'
  material = materials.Hollomon(labels = labels, E = E, nu = nu, sy = sy, n=n)
  
     
model = CuboidTest_VER(lx =lx, ly = ly, lz = lz, Nx = Nx, Ny = Ny, Nz = Nz, abqlauncher = abqlauncher, label = label, workdir = workdir, material = material, compart = compart, force = force, force_fin = force_fin, disp = disp, loading = loading, elType = elType, is_3D = is_3D, cpus = cpus, export_fields = export_fields, unloading_reloading = unloading_reloading, lateralbc = lateralbc)


model.MakeMesh()
mesh = model.mesh
centroids = mesh.centroids()
nodes = mesh.nodes
nodes_postion = np.array([nodes.x, nodes.y, nodes.z]).transpose()
conn = np.array(mesh.connectivity)
labels = nodes.labels

# VORONOI CELL BUILDING

# Bounding box
xmin, xmax = nodes_postion[:,0].min(), nodes_postion[:,0].max()
ymin, ymax = nodes_postion[:,1].min(), nodes_postion[:,1].max() 
zmin, zmax = nodes_postion[:,2].min(), nodes_postion[:,2].max() 

# Seeds
seeds = np.random.rand(Nseed,3) * np.array([[xmax - xmin, ymax -ymin, zmax - zmin]]) + np.array([[xmin, ymin, zmin]])
seed_flags = np.arange(Nseed)
elem_flags = interpolate.griddata(seeds, seed_flags, centroids, method = "nearest")

# Building sets
elabels = np.array(mesh.labels)
grain_sets = []
for i in xrange(len(seed_flags)):
  flag = seed_flags[i]
  grain_sets.append(elabels[np.where(elem_flags == flag)[0]])
  mesh.add_set(label = "Kikinou{0}".format(i), elements = elabels[np.where(elem_flags == flag)[0]])
  




labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
material = [materials.Bilinear(labels = labels[i],E = E[i], nu = nu[i], n = n[i],sy = sy[i]) for i in xrange(Ne)]
model.material = material
sy_field = FieldOutput(labels = mesh.labels, data = sy, position = "element")


model.MakeInp()


