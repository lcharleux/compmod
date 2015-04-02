from compmod.models import CuboidTest
from abapy import materials
from abapy.misc import load
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform


  
def field_func(outputs, step):
    """
    A function that defines the scalar field you want to plot
    """
    epsilon = np.array(outputs['field']['LE'][step].get_component(22).data)
    return epsilon
  
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


#PARAMETERS
lx, ly = 1., 1.
Nx, Ny = 20, 20 
Ne = Nx * Ny
disp = .2
nFrames = 20
workdir = "workdir/"
label0 = "cuboidTest_homo0"
label1 = "cuboidTest_homo1"

elType = "CPE4"
lateralbc = { "right":"pseudohomo", "left":"pseudohomo" }
cpus = 1
abqlauncher = None
node = platform.node()
if node ==  'lcharleux':      
  abqlauncher   = '/opt/Abaqus/6.9/Commands/abaqus' # Local machine configuration
if node ==  'serv2-ms-symme': 
  abqlauncher   = '/opt/abaqus/Commands/abaqus' # Local machine configuration
if node ==  'epua-pd47': 
  abqlauncher   = 'C:\SIMULIA\Abaqus\6.13-1'
if node ==  'epua-pd45': 
  abqlauncher   = 'C:\SIMULIA/Abaqus/Commands/abaqus'  
compart = True
run_simulation = True

if compart:
  E  = 1. * np.ones(Ne) # Young's modulus
  nu = .3 * np.ones(Ne) # Poisson's ratio
  sy_mean = .01
  sy = np.random.rayleigh(sy_mean, Ne)
  labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
  material = [materials.VonMises(labels = labels[i], E = E[i], nu = nu[i], sy = sy[i]) for i in xrange(Ne)]
else:
  E = 1.
  nu =.3
  sy = .01
  labels = 'SAMPLE_MAT'
  material = materials.VonMises(labels = labels, E = E, nu = nu, sy = sy)

m0 = CuboidTest(lx =lx, ly = ly, Nx = Nx, Ny = Ny, abqlauncher = abqlauncher, label = label0, workdir = workdir, cpus = cpus, material = material, compart = compart, disp = disp, elType = elType, lateralbc = lateralbc)
m1 = CuboidTest(lx =lx, ly = ly, Nx = Nx, Ny = Ny, abqlauncher = abqlauncher, label = label1, workdir = workdir, cpus = cpus, material = material, compart = compart, disp = disp, elType = elType)

if run_simulation:
  m0.MakeInp()
  m1.MakeInp()
  m0.Run()
  m1.Run()
  m0.PostProc()
  m1.PostProc()
else:
  m0.LoadResults()
  m1.LoadResults()

# Plotting results
if m0.outputs['completed'] and m1.outputs['completed']:
    
    disp0 =  np.array(m0.outputs['history']['disp'].values()[0].data[0])
    force0 =  np.array(np.array(m0.outputs['history']['force'].values()).sum().data[0])
    volume0 = np.array(np.array(m0.outputs['history']['volume'].values()).sum().data[0])
    length0 = ly + disp0
    surface0 = volume0 / length0
    logstrain0 = np.log10(1. + disp0 / ly)
    linstrain0 = disp0/ly
    strain0 = linstrain0
    stress0 = force0 / surface0 
    
    disp1 =  np.array(m1.outputs['history']['disp'].values()[0].data[0])
    force1 =  np.array(np.array(m1.outputs['history']['force'].values()).sum().data[0])
    volume1 = np.array(np.array(m1.outputs['history']['volume'].values()).sum().data[0])
    length1 = ly + disp1
    surface1 = volume1 / length1
    logstrain1 = np.log10(1. + disp1 / ly)
    linstrain1 = disp1/ly
    strain1 = linstrain1
    stress1 = force1 / surface1 
    
    fig = plt.figure(0)
    plt.clf()
    sp1 = fig.add_subplot(2, 1, 1)
    plt.plot(disp0, force0, 'or-', label = "Pseudo-homo")
    plt.plot(disp1, force1, 'ob-', label = "Free")
    plt.xlabel('Displacement, $U$')
    plt.ylabel('Force, $F$')
    plt.grid()
    plt.legend()
    sp1 = fig.add_subplot(2, 1, 2)
    plt.plot(strain0, stress0, 'or-', label = "Pseudo-homo")
    plt.plot(strain1, stress1, 'ob-', label = "Free")
    plt.xlabel('Tensile Strain, $\epsilon$')
    plt.ylabel(' Tensile Stress $\sigma$')
    plt.grid()
    plt.legend()
    #plt.savefig(workdir + label + 'history.pdf')
    
  # Field Outputs
     
          
      

    mesh0 = m0.outputs['mesh']
    mesh1 = m1.outputs['mesh']
    max_strain0 = strain0.max()
    max_strain1 = strain1.max()
    fig = plt.figure("Fields")
    plt.clf()
    ax = fig.add_subplot(1, 2, 1)
    ax.set_aspect('equal')
    ax.axis("off")
    plt.grid()
    plot_mesh(ax, mesh0, m0.outputs, 0, field_func, cbar_label = 'Tensile Strain, $\epsilon$')
      #plot_mesh(ax, mesh, outputs, 0, field_func = None, cbar = False, disp = False)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    ax = fig.add_subplot(1, 2, 2)
    ax.set_aspect('equal')
    ax.axis("off")
    plt.grid()
    plot_mesh(ax, mesh1, m1.outputs, 0, field_func, cbar_label = r'Tensile Strain, $\epsilon$')
      #plot_mesh(ax, mesh, outputs, 0, field_func = None, cbar = False, disp = False)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    #plt.savefig(workdir + label + '_fields.pdf')
    plt.show()
    
 
