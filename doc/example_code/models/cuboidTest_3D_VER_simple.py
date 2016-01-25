'''
This routine make a tensile test along Y axis on a cuboid with specific boundary conditions.
It is possible to pilot this test in force or displacement (loading = {force} or {displacement})

In case of force displacement, it is possible to have unloading and reloading conditions. In this case (unloading_reloading = True), "force" is the force of the beginning of the unloading step (end of step when F=0). Then "force_fin" is the force of reloading step (end of the step when F = force_fin)

At least, cuboid must have "pseudohomo" conditions on "top"  side. it involves that bottom nodes have "pseudhomo" conditions (i.e same displacement along Y axis)
It is possible for the other sides to have 3 types of boundary conditions:
* none boundary conditions, then lateralbc = {}
* "pseudohomo" conditions : same displacement along the normal vector of the side
* "periodic" conditions : the nodes of two opposide side have the same displacement (i.e right and left sides are linked and/or front and rear sides are linked)

In case of single periodic conditions, (only two sides linked): lateralbc = {"right":"periodic","top":"pseudohomo"}
In case of double periodic conditions, (right/left sides and front/rear sides): lateralbc = {"right":"periodic", "front":"periodic", "top":"pseudohomo"}

It is possible to mix boundary conditions : for example to have periodic conditions on right/left sides and pseudohomo conditions on front side. For this example, indicate : lateralbc = {"right":"periodic", "front":"pseudohomo", "top":"pseudohomo"}

It is also possible to launch a specify number of simulation (variable iteration). It can be used for determinate the varibility to given parameters for example. If just one simulation is needed, then specify : iteration = 1

'''

from compmod.models import CuboidTest_VER
from abapy import materials
import matplotlib.pyplot as plt
import numpy as np
import platform
from scipy import interpolate


#PARAMETERS
iteration = 1 #number of simulations
lx, ly, lz = 1., 2., 1.
Nx, Ny, Nz = 4,4,4
Ne = Nx * Ny * Nz
is_3D = True
steps =  []
steps.append({
  "name": "LOADING", 
  "control": "displacement",
  "frames": 100,
  "value": 0.04})
steps.append({
  "name": "UNLOADING", 
  "control": "force",
  "frames": 100,
  "value": 0.})

export_fields = False
compart = True

#lateralbc = {}
lateralbc = {"right":"periodic","front":"periodic","top":"pseudohomo"}
lateralbc = {"top":"pseudohomo"}
workdir = "workdir/"
label = "cuboidTest_3D_VER"
if is_3D:
  elType = "C3D8"
else:
  elType = "CPS4"
abqlauncher = "" 
cpus = 1 
node = platform.node()
if node ==  'lcharleux':      
  abqlauncher   = '/opt/Abaqus/6.9/Commands/abaqus' # Local machine configuration
if node ==  'serv2-ms-symme': 
  abqlauncher   = '/opt/abaqus/Commands/abaqus'# Local machine configuration
  cpus = 6
if node ==  'epua-pd47': 
  abqlauncher   = 'C:/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe' # Local machine configuration
  cpus = 1
if node ==  'epua-pd45': 
  abqlauncher   = 'C:\SIMULIA/Abaqus/Commands/abaqus'
if node ==  'SERV3-MS-SYMME': 
  abqlauncher   = '"C:/Program Files (x86)/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe"' # Local machine configuration
  cpus = 6



if compart:
  E  = 64000. * np.ones(Ne) # Young's modulus
  nu = .3 * np.ones(Ne) # Poisson's ratio
  Ssat = 673. * np.ones(Ne)
  n = 511.18 * np.ones(Ne)
  sy_mean = 150.
  ray_param = sy_mean/1.253314
  sy = np.random.rayleigh(ray_param, Ne)
  labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
  material = [materials.Bilinear(labels = labels[i], E = E[i], nu = nu[i], Ssat = Ssat[i], n=n[i], sy = sy[i]) for i in xrange(Ne)]
else:
  E = 64000.
  nu =.3
  sy = 145.8
  n = .081
  labels = 'SAMPLE_MAT'
  material = materials.Hollomon(labels = labels, E = E, nu = nu, sy = sy, n=n)

     
m =CuboidTest_VER(lx =lx, ly = ly, lz = lz, Nx = Nx, Ny = Ny, Nz = Nz, abqlauncher = abqlauncher, label = label, workdir = workdir, material = material, compart = compart, steps=steps, elType = elType, is_3D = is_3D, cpus = cpus, export_fields = export_fields,  lateralbc = lateralbc)
m.MakeInp()
#m.Run()
#m.MakePostProc()
#m.RunPostProc()
m.LoadResults()



# Plotting results
if m.outputs['completed']:
  U =  np.array(m.outputs['history']['U'].values()[0].data)
  RF =  np.array(m.outputs['history']['RF'].values()[0].data)
  CF =  np.array(m.outputs['history']['CF'].values()[0].data)
  volume = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[0])
  F = RF + CF
      
  length = ly + U
  surface = volume / length
  logstrain = np.log10(1. + U / ly)
  linstrain = U/ly
  strain = linstrain
  stress = F / surface
  
fig = plt.figure(0)
plt.clf()
for i in xrange(len(U)):
  plt.plot(strain[i], stress[i], label = "Step {0}".format(i))  
plt.grid()  
plt.legend()
plt.show()
  

