from compmod.models import RingCompression
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy

#PAREMETERS
inner_radius, outer_radius = 4.5, 5.
Nt, Nr = 50, 10 
disp = 3.
nFrames = 100
workdir = "workdir/"
label = "ringCompression"
elType = "CPS4"

#TASKS
run_sim = True
plot = True

#MODEL DEFINITION
m = RingCompression(
  inner_radius = inner_radius, 
  outer_radius = outer_radius, 
  disp = disp, 
  nFrames = nFrames, 
  Nr = Nr, 
  Nt = Nt, 
  workdir = workdir,
  label = label, 
  elType = elType)

# SIMULATION
m.MakeMesh()
if run_sim:
  m.MakeInp()
  m.Run()
  m.PostProc()

# SOME PLOTS
mesh = m.mesh
f = open(workdir + label + '.pckl')
outputs = pickle.load(f)
f.close()

# Fields
def field_func(outputs, step):
  """
  A function that defines the scalar field you want to plot
  """
  return outputs['field']['S'][step].vonmises()

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
force = -2. * outputs['history']['force']
disp = -2. * outputs['history']['disp']

fig = plt.figure('Load vs. disp')
plt.clf()
plt.plot(disp.data[0], force.data[0], 'or-', label = 'Loading', linewidth = 2.)
plt.plot(disp.data[1], force.data[1], 'vb-', label = 'Unloading', linewidth = 2.)
plt.legend(loc="upper left")
plt.grid()
plt.xlabel('Diplacement, $U$')
plt.ylabel('Force, $F$')
plt.savefig(workdir + label + '_load-vs-disp.pdf')




