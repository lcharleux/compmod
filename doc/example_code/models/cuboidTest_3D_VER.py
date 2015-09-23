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

settings = {}
settings['file_name'] = 'Courbe_ref_alu_transv.txt' # experimental data
strain_exp, stress_exp = read_file(settings['file_name'])

#PARAMETERS
iteration = 1 #number of simulations
lx, ly, lz = 1., 2., 1.
Nx, Ny, Nz = 15,30, 15
Ne = Nx * Ny * Nz
is_3D = True
loading = {"displacement"} #"loading" : force or displacement
disp = strain_exp[-1] * ly #if the simulation is comparated to experimental data, activate this line
#disp = 0.1
force = 190.
force_fin = force + 20.
nFrames = 30
export_fields = False
compart = False
unloading_reloading = False #for one cycle of loading (F = force), unloding (F=0) and reloading (F = force_fin) 

#lateralbc = {}
#lateralbc = {"right":"periodic","front":"periodic","top":"pseudohomo"}
lateralbc = {"top":"pseudohomo"}
workdir = "workdir/"
label = "cuboidTest_3D_VER"
if is_3D:
  elType = "C3D8"
else:
  elType = "CPS4"
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


strain_tot, stress_tot =[], []
for i in xrange(iteration): 
  if compart:
    E  = 64000. * np.ones(Ne) # Young's modulus
    nu = .3 * np.ones(Ne) # Poisson's ratio
    sy_mean = 174.45 * np.ones(Ne)
    Ssat = 673.79 * np.ones(Ne)
    n = 511.18 * np.ones(Ne)
    #n = 1. * np.ones(Ne)
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
  
     
  m =CuboidTest_VER(lx =lx, ly = ly, lz = lz, Nx = Nx, Ny = Ny, Nz = Nz, abqlauncher = abqlauncher, label = label, workdir = workdir, material = material, compart = compart, force = force, force_fin = force_fin, disp = disp, loading = loading, elType = elType, is_3D = is_3D, cpus = cpus, export_fields = export_fields, unloading_reloading = unloading_reloading, lateralbc = lateralbc)
  m.MakeInp()
  m.Run()
  m.MakePostProc()
  m.RunPostProc()
  
  
  
  
  # Plotting results
  if m.outputs['completed']:
  
    # History Outputs
    
    if "displacement" in loading:
      displ =  np.array(m.outputs['history']['disp'].values()[0].data[0])
      force =  np.array(np.array(m.outputs['history']['force'].values()).sum().data[0])
      volume = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[0])
    if "force" in loading and unloading_reloading == False:
      displ =  np.array(m.outputs['history']['disp'].values()[0].data[0])
      force =  np.array(np.array(m.outputs['history']['load'].values()).sum().data[0])
      volume = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[0])
    if "force" in loading and unloading_reloading:
      volume_loading = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[0])
      volume_unloading = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[1])
      volume_reloading = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[2])
      displ_loading =  np.array(m.outputs['history']['disp'].values()[0].data[0])
      displ_unloading =  np.array(m.outputs['history']['disp'].values()[0].data[1])
      displ_reloading =  np.array(m.outputs['history']['disp'].values()[0].data[2])
      force_loading =  np.array(np.array(m.outputs['history']['load'].values()).sum().data[0])
      force_unloading =  np.array(np.array(m.outputs['history']['load'].values()).sum().data[1])
      force_reloading =  np.array(np.array(m.outputs['history']['load'].values()).sum().data[2])
      A, B, C = [],[], []
      for i in xrange(len(displ_loading)):
        A.append(displ_loading[i])
        B.append(force_loading[i])
        C.append(volume_loading[i])
      for i in xrange(len(displ_unloading)):
        A.append(displ_unloading[i])
        B.append(force_unloading[i])
        C.append(volume_unloading[i])
      for i in xrange(len(displ_reloading)):
        A.append(displ_reloading[i])
        B.append(force_reloading[i])
        C.append(volume_reloading[i])
      displ = np.array(A)
      force = np.array(B)
      volume = np.array(C)
        
    length = ly + displ
    surface = volume / length
    logstrain = np.log10(1. + displ / ly)
    linstrain = displ/ly
    strain = linstrain
    stress = force / surface
    
    strain_tot.append(strain)
    stress_tot.append(stress)


strain_min = min(strain_tot[0][-1], strain_exp[-1])
strain_grid = np.linspace(0., strain_min, 100)
g = interpolate.interp1d(strain_exp, stress_exp)
interp_stress_exp = g(strain_grid)
interp_stress_tot, err = [], []
for i in xrange(len(strain_tot)):      
    f = interpolate.interp1d(strain_tot[i], stress_tot[i])
    interp_stress_tot.append(f(strain_grid))
for i in xrange(len(strain_tot)):
    err.append(np.sqrt(((interp_stress_tot[i] - interp_stress_exp)**2).sum()))

     
fig = plt.figure(1)
plt.clf()
#  sp1 = fig.add_subplot(2, 1, 1)
#  plt.plot(strain, stress, 'k-')
#  plt.xlabel('Displacement, $U$')
#  plt.ylabel('Force, $F$')
#  plt.grid()
#  sp1 = fig.add_subplot(2, 1, 2)
for i in xrange(len(strain_tot)):
  plt.plot(strain_tot[i], stress_tot[i], 'r-', label = 'simulation curve', linewidth = 2.)
  plt.plot(strain_exp, stress_exp, 'ko', label = 'experimental curve', markersize = 5., markevery=10)
  plt.xlabel('Tensile Strain, $\epsilon \ (\%)$',fontsize=16)
  plt.ylabel(' Tensile Stress $\sigma \ (MPa)$',fontsize=16)
  plt.legend(loc="lower right")
  plt.grid()
plt.savefig(workdir + label + 'history.pdf')
  
'''
# Field Outputs
 
def field_func(outputs, step):
  
  A function that defines the scalar field you want to plot
  """
  epsilon = np.array(outputs['field']['LE'][step].get_component(22).data)
  return (epsilon - max_strain) / max_strain

def plot_mesh(ax, mesh, outputs, step, field_func =None, cbar = True, cbar_label = 'Z', cbar_orientation = 'horizontal', disp = True):
  """
  A function that plots the deformed mesh with a given field on it.
  """
  mesh2 = copy.deepcopy(mesh)
  if disp:
    U = outputs['field']['U'][step]
    mesh2.nodes.apply_displacement(U)
  X,Y,Z,tri = mesh2.dump2triplot()
  xb,yb,zb = mesh2.get_border() 
  xe, ye, ze = mesh2.get_edges()
  ax.plot(xb, yb,'k-', linewidth = 2.)
  ax.plot(xe, ye,'k-', linewidth = .5)
  if field_func != None:
    field = field_func(outputs, step)
    grad = ax.tricontourf( X, Y, tri, field)
    if cbar :
      bar = plt.colorbar(grad, orientation = cbar_orientation)
      bar.set_label(cbar_label)
    

outputs = m.outputs
mesh = outputs['mesh']
max_strain = strain.max()
fig = plt.figure("Fields")
plt.clf()
ax = fig.add_subplot(1, 1, 1)
ax.set_aspect('equal')
plt.grid()
plot_mesh(ax, mesh, outputs, 0, field_func, cbar_label = r'Relative Tensile Strain, $\frac{\epsilon - \epsilon_{av}}{\epsilon_{av}}$')
#plot_mesh(ax, mesh, outputs, 0, field_func = None, cbar = False, disp = False)
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.savefig(workdir + label + '_fields.pdf')
'''
