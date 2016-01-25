from compmod.models import RingCompression
from abapy import materials
from abapy.misc import load
from abapy.postproc import FieldOutput
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform
from scipy import interpolate

#PAREMETERS
is_3D = True
inner_radius, outer_radius = 30 , 40
Nt, Nr, Na = 10, 5, 2
Ne = Nt * Nr * Na
displacement = 5.
nFrames = 100
Nseed = 400 # Number of grain seeds
sigma_0_hp = .001
k_hp = .001
nu = 0.3
n = 0.001
E = 1.
thickness =10.
workdir = "workdir/"
label = "ringCompression_voronoi"
elType = "C3D8"
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
  
#TASKS
Run_simu = True



#MODEL DEFINITION
disp = displacement/2.
model = RingCompression( 
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
  is_3D = is_3D,
  compart = True)

# SIMULATION
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
  

# Hall-Petch (not mandatory)
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
  f = open(label + ".vtk", "w")
  f.write(mesh.dump2vtk())
  f.write(sy_field.dump2vtk(name = "Yield_Stress"))
  f.write( model.outputs['field']['S'][0].vonmises().dump2vtk(name = "Von_Mises_Stress"))
  f.close()
  
  # History Outputs
  force = -2. * model.outputs['history']['force']
  disp = -2. * model.outputs['history']['disp']
 
   
  fig = plt.figure(0)
  plt.clf()
  sp1 = fig.add_subplot(1, 1, 1)
  plt.plot(disp.data[0], force.data[0], 'ro-', label = 'Loading', linewidth = 2.)
  plt.plot(disp.data[1], force.data[1], 'bv-', label = 'Unloading', linewidth = 2.)
  plt.legend(loc="upper left")
  plt.xlabel('Displacement, $U$')
  plt.ylabel('Force, $F$')
  plt.grid()
  plt.savefig(workdir + label + 'history.pdf')
