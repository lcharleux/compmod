from compmod.models import CuboidTest_VER
from compmod.distributions import Triangular 
from abapy import materials
from abapy.misc import load
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform

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
settings['file_name'] = 'cuivre_cufe2p_ANR.txt' # experimental data
strain_exp, stress_exp = read_file(settings['file_name'])

#PARAMETERS
lx, ly, lz = 1., 1., 1.
Nx, Ny, Nz = 10, 10, 10
Ne = Nx * Ny * Nz

loading = {"force"} #"loading" : force or displacement
#disp = strain_exp[-1] * ly
disp = 0.1
force = 190.
nFrames = 30
export_fields = False
compart = True
unloading_reloading = False #for one cycle of loading (F), unloding (F=0) and reloading (F+20N) 
lateralbc = {}
#lateralbc = { "right":"pseudohomo"}
workdir = "workdir/"
label = "cuboidTest_3D_VER"
elType = "C3D8"
cpus = 1
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



if compart:
  E  = 72000. * np.ones(Ne) # Young's modulus
  nu = .3 * np.ones(Ne) # Poisson's ratio
  sy_mean = 184.791 * np.ones(Ne)
  Ssat = 1031.394 * np.ones(Ne)
  n = 333.485 * np.ones(Ne)
  ray_param = sy_mean/1.253314
  sy = np.random.rayleigh(ray_param, Ne)
  labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
  material = [materials.Bilinear(labels = labels[i], E = E[i], nu = nu[i], Ssat = Ssat[i], n=n[i], sy = sy[i]) for i in xrange(Ne)]
else:
  E = 70000.
  nu =.3
  sy = 150
  n = .1
  labels = 'SAMPLE_MAT'
  material = materials.Hollomon(labels = labels, E = E, nu = nu, sy = sy, n=n)
      
m =CuboidTest_VER(lx =lx, ly = ly, lz = lz, Nx = Nx, Ny = Ny, Nz = Nz, abqlauncher = abqlauncher, label = label, workdir = workdir, material = material, compart = compart, force = force, disp = disp, loading = loading, elType = elType, is_3D = True, cpus = cpus, export_fields = export_fields, unloading_reloading = unloading_reloading, lateralbc = lateralbc)
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
   
  fig = plt.figure(0)
  plt.clf()
#  sp1 = fig.add_subplot(2, 1, 1)
#  plt.plot(strain, stress, 'k-')
#  plt.xlabel('Displacement, $U$')
#  plt.ylabel('Force, $F$')
#  plt.grid()
#  sp1 = fig.add_subplot(2, 1, 2)
  plt.plot(strain, stress, 'k-', label = 'simulation curve', linewidth = 2.)
#  plt.plot(strain_exp, stress_exp, 'r-', label = 'experimental curve', linewidth = 2.)
  plt.xlabel('Tensile Strain, $\epsilon$')
  plt.ylabel(' Tensile Stress $\sigma$')
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
