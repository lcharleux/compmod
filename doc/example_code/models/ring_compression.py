from compmod.models import RingCompression
from abapy.materials import Hollomon
from abapy.misc import load
from abapy import materials
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform


#PAREMETERS
cpus = 6
compart = False #True for a compartimentalized model
is_3D = True #True for a 3D simulation
unloading = False #True il the unloading part of the simulation is needed
export_fields = False #True if stress and strain fields are needed
inner_radius, outer_radius = 45.96 , 50
thickness = 15
Nt, Nr, Na = 80, 8, 15
if is_3D == False :
  Ne = Nt * Nr
  elType = "CPS4" #CPS4 for plane strain element, CPS4 for plane stress elements
else:
  Ne = Nt * Nr * Na
  elType = "C3D8"
displacement = 45.
nFrames = 100

E = 64000.
nu = .3

if compart == False: 
  sy = 200.
  n = 0.0875
  material = Hollomon(labels = "SAMPLE_MAT", E = E, nu = nu, sy = sy, n = n)
else:
  E_array = E * np.ones(Ne) # Young's modulus
  nu_array = nu * np.ones(Ne) # Poisson's ratio
  Ssat =1000 * np.ones(Ne)
  n = 200 * np.ones(Ne)
  sy_mean = 200.
  ray_param = sy_mean/1.253314
  sy = np.random.rayleigh(ray_param, Ne)
  labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
  material = [materials.Bilinear(labels = labels[i], E = E_array[i], nu = nu_array[i], Ssat = Ssat[i], n=n[i], sy = sy[i]) for i in xrange(Ne)]


workdir = "workdir/"
label = "ringCompression"

filename = 'force_vs_disp_ring1.txt'

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
run_sim = True
plot = True

def read_file(file_name):
  '''
  Read a two rows data file and converts it to numbers
  '''
  f = open(file_name, 'r') # Opening the file
  lignes = f.readlines() # Reads all lines one by one and stores them in a list
  f.close() # Closing the file
#    lignes.pop(0) # Delete le saut de ligne for each lines
  force_exp, disp_exp = [],[]

  for ligne in lignes:
      data = ligne.split() # Lines are splitted
      disp_exp.append(float(data[0]))
      force_exp.append(float(data[1]))
  return -np.array(disp_exp), -np.array(force_exp)


disp_exp, force_exp = read_file(filename)

#MODEL DEFINITION
disp = displacement/2
thick = thickness/2

m = RingCompression( material = material , 
      inner_radius = inner_radius, 
      outer_radius = outer_radius, 
      disp = disp,
      thickness = thick,
      nFrames = nFrames, 
      Nr = Nr, 
      Nt = Nt, 
      Na = Na,
      unloading = unloading,
      export_fields = export_fields,
      workdir = workdir,
      label = label, 
      elType = elType,
      abqlauncher = abqlauncher,
      cpus = cpus,
      is_3D = is_3D,
      compart = compart)

# SIMULATION
m.MakeMesh()
if run_sim:
  m.MakeInp()
  m.Run()
  m.PostProc()

# SOME PLOTS
mesh = m.mesh
outputs = load(workdir + label + '.pckl')

if outputs['completed']:
  # Fields
  if export_fields == True :
    def field_func(outputs, step):
          """
          A function that defines the scalar field you want to plot
          """
          return outputs['field']['S'][step].vonmises()
    if is_3D == False:
        def plot_mesh(ax, mesh, outputs, step, field_func =None, zone = 'upper right', cbar = True, cbar_label = 'Z', cbar_orientation = 'horizontal', disp = True):
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
          if zone == "upper right": kx, ky = 1., 1.
          if zone == "upper left": kx, ky = -1., 1.
          if zone == "lower right": kx, ky = 1., -1.
          if zone == "lower left": kx, ky = -1., -1.
          ax.plot(kx * xb, ky * yb,'k-', linewidth = 2.)
          ax.plot(kx * xe, ky * ye,'k-', linewidth = .5)
          if field_func != None:
            field = field_func(outputs, step)
            grad = ax.tricontourf(kx * X, ky * Y, tri, field.data)
            if cbar :
              bar = plt.colorbar(grad, orientation = cbar_orientation)
              bar.set_label(cbar_label)
        
        # Exp data
        disp_exp, force_exp = read_file("test_expD2.txt")    
        
        fig = plt.figure("Fields")
        plt.clf()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_aspect('equal')
        plt.grid()
        plot_mesh(ax, mesh, outputs, 0, field_func, cbar_label = '$\sigma_{eq}$')
        plot_mesh(ax, mesh, outputs, 0, field_func = None, cbar = False, disp = False)
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.savefig(workdir + label + '_fields.pdf')
        

  # Load vs disp
  force = -4. * outputs['history']['force']
  disp = -2. * outputs['history']['disp']
    
  fig = plt.figure('Load vs. disp')
  plt.clf()
  plt.plot(disp.data[0], force.data[0], 'ro-', label = 'Loading', linewidth = 2.)
  if unloading == True : plt.plot(disp.data[1], force.data[1], 'bv-', label = 'Unloading', linewidth = 2.)
  plt.plot(disp_exp, force_exp, 'k-', label = 'Exp', linewidth = 2.)
  plt.legend(loc="upper left")
  plt.grid()
  plt.xlabel('Displacement, $U$')
  plt.ylabel('Force, $F$')
  plt.savefig(workdir + label + '_load-vs-disp.pdf')
  
else: 
  print 'Simulation not completed'



