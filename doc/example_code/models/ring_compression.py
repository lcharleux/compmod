from compmod.models import RingCompression
from abapy.materials import Hollomon
from abapy.misc import load
from abapy import materials
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform
from scipy import interpolate
node = platform.node()

#PAREMETERS
compart = False #True for a compartimentalized model
is_3D = True #True for a 3D simulation
unloading = False #True il the unloading part of the simulation is needed
export_fields = False #True if stress and strain fields are needed
inner_radius, outer_radius = 45.96 , 50
thickness = 15.
Nt, Nr, Na = 80, 10, 12
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
  sy = 145.8
  n = 0.081
  material = Hollomon(labels = "SAMPLE_MAT", E = E, nu = nu, sy = sy, n = n)
else:
  E_array = E * np.ones(Ne) # Young's modulus
  nu_array = nu * np.ones(Ne) # Poisson's ratio
  Ssat = 673.79 * np.ones(Ne)
  n = 511.2 * np.ones(Ne)
  sy_mean = 174.45
  ray_param = sy_mean/1.253314
  sy = np.random.rayleigh(ray_param, Ne)
  labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
  material = [materials.Bilinear(labels = labels[i], E = E_array[i], nu = nu_array[i], Ssat = Ssat[i], n=n[i], sy = sy[i]) for i in xrange(Ne)]


workdir = "workdir/"
label = "ringCompression"

filename = 'force_vs_disp_ring1.txt'

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
  data1, data2 = [],[]

  for ligne in lignes:
      data = ligne.split() # Lines are splitted
      data1.append(float(data[0]))
      data2.append(float(data[1]))
  return -np.array(data1), -np.array(data2)


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
        X_exp, Y_exp = read_file("test_expD2.txt")    
        
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
  
    else:#for ploting
      
      mesh2 = copy.deepcopy(mesh)
      top_external_nodes = [val for val in m.mesh.nodes.sets['top'] if val in m.mesh.nodes.sets['external_nodes']]
      top_lateral_nodes = [val for val in m.mesh.nodes.sets['top'] if val in m.mesh.nodes.sets['lateral_nodes']]
      top_edge_nodes = []
      if unloading == True:
        step = 1
      else:
        step = 0
      U = outputs['field']['U'][step]
      mesh2.nodes.apply_displacement(U)
      x_coord = mesh2.nodes.x
      y_coord = mesh2.nodes.y
      X_lateral, Y_lateral, X_external, Y_external = [], [], [], []
      X,Y = [],[]
      for i in top_external_nodes:
        X_external.append(x_coord[i-1])
        Y_external.append(y_coord[i-1])
      for i in top_lateral_nodes:
        X_lateral.append(x_coord[i-1])
        Y_lateral.append(y_coord[i-1])
        
      for i in xrange(len(X_external)):
        if Y_lateral[i] <= Y_external[i]:
          X.append(X_external[i])
          Y.append(Y_external[i])
        else:
          X.append(X_lateral[i])
          Y.append(Y_lateral[i])
      
      filename1 = 'geometrie_anneau1.xyz'
      X_exp, Y_exp = read_file("geometrie_anneau1.xyz")    
      X_exp_fin = X_exp + 124.
      Y_exp_fin = Y_exp + 96.12343
      fig = plt.figure("Deformed shape")
      plt.clf()
      ax = fig.add_subplot(1, 1, 1)
      ax.set_aspect('equal')
      plt.grid()
      ax.set_ylim([-10,40])
      ax.set_xlim([-10,70])
      plt.plot(X, Y, 'b-', label = 'Simulated shape', linewidth = 2.)
      plt.plot(-X_exp_fin, Y_exp_fin, 'ro', label = 'Experimental shape', markersize = 5., markevery=2)
      plt.legend(loc="lower left")
      plt.xlabel('$x\ (mm)$',fontsize=18)
      plt.ylabel('$y\ (mm)$',fontsize=18)
      plt.savefig(workdir + label + '_deformed_shape.pdf') 
          
         
  # Load vs disp
  force = -4. * outputs['history']['force']
  disp = -2. * outputs['history']['disp']
    
  fig = plt.figure('Load vs. disp')
  plt.clf()
  plt.plot(disp.data[0], force.data[0], 'r-', label = 'Loading',  linewidth = 2.)
  if unloading == True : plt.plot(disp.data[1], force.data[1], 'b-', label = 'Unloading',  linewidth = 2.)
  plt.plot(disp_exp, force_exp, 'ko', label = 'Exp', markersize = 5., markevery=10)
  plt.legend(loc="upper left")
  plt.grid()
  plt.xlabel('Displacement, $U\ (mm)$',fontsize=16)
  plt.ylabel('Force, $F\ (N)$',fontsize=16)
  plt.savefig(workdir + label + '_load-vs-disp.pdf')
  
  file = open("force_deplacement_hollo_CPS4_200_24.txt", "w")
  for i in xrange(len(disp.data[0])):
      file.write(str(disp.data[0][i])+"\t"+str(force.data[0][i])+"\n")
  file.close()

  disp_min = min(disp.data[0][-1], disp_exp[-1])
  disp_grid = np.linspace(0., disp_min, 100)
  g = interpolate.interp1d(disp_exp, force_exp)
  interp_force_exp = g(disp_grid)    
  f = interpolate.interp1d(disp.data[0], force.data[0])
  interp_force = f(disp_grid)
  err = np.sqrt(((interp_force - interp_force_exp)**2).sum())
  
else: 
  print 'Simulation not completed'



