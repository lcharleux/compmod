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
Nx, Ny, Nz = 2, 2, 2
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

#model = CuboidTest(lx = lx, ly = ly, lz = lz, Nx = Nx, Ny = Ny, Nz = Nz, abqlauncher = abqlauncher, label = label, workdir = workdir, cpus = cpus, compart = compart, disp = disp, elType = elType, is_3D = True)
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

# Volumes and dimensions
elem_volume = mesh.volume()
grain_volume = np.zeros(len(seed_flags))
for i in xrange(len(seed_flags)):
  flag = seed_flags[i]
  grain_volume[i] = ((elem_flags == flag) * elem_volume).sum()
grain_diameter = (grain_volume * 6./np.pi)**(1. / 3.)   
elem_grain_diameter = grain_diameter[elem_flags]
  

# Hall-Petch
sy = sigma_0_hp + k_hp / elem_grain_diameter**.5

E  = E * np.ones(Ne) # Young's modulus
nu = nu * np.ones(Ne) # Poisson's ratio
n = n * np.ones(Ne)


labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
material = [materials.Bilinear(labels = labels[i],E = E[i], nu = nu[i], n = n[i],sy = sy[i]) for i in xrange(Ne)]
model.material = material
sy_field = FieldOutput(labels = mesh.labels, data = sy, position = "element")

if Run_simu:
  model.MakeInp()
  model.Run()
  model.MakePostProc()
  model.RunPostProc()
else:
  model.LoadResults()
# Plotting results
if model.outputs['completed']:
  U = model.outputs['field']['U'][0]
  mesh.nodes.apply_displacement(U)
  f = open("cuboidTest_voronoi.vtk", "w")
  f.write(mesh.dump2vtk())
  f.write(sy_field.dump2vtk(name = "Yield_Stress"))
  f.write( model.outputs['field']['S'][0].vonmises().dump2vtk(name = "Von_Mises_Stress"))
  f.close()
  
  # History Outputs
  disp =  np.array(model.outputs['history']['disp'].values()[0].data[0])
  force =  np.array(np.array(model.outputs['history']['force'].values()).sum().data[0])
  volume = np.array(np.array(model.outputs['history']['volume'].values()).sum().data[0])
  length = ly + disp
  surface = volume / length
  logstrain = np.log10(1. + disp / ly)
  linstrain = disp/ly
  strain = linstrain
  stress = force / surface 
   
  fig = plt.figure(0)
  plt.clf()
  sp1 = fig.add_subplot(2, 1, 1)
  plt.plot(disp, force, 'ok-')
  plt.xlabel('Displacement, $U$')
  plt.ylabel('Force, $F$')
  plt.grid()
  sp1 = fig.add_subplot(2, 1, 2)
  plt.plot(strain, stress, 'ok-')
  plt.xlabel('Tensile Strain, $\epsilon$')
  plt.ylabel(' Tensile Stress $\sigma$')
  plt.grid()
  plt.savefig(workdir + label + 'history.pdf')


