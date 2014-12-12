# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 15:35:38 2014

@author: issackm
"""

from compmod.models import CuboidTest
from abapy import materials
from abapy.misc import load
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform



#PARAMETERS

lx, ly, lz = 1., 1., 1.
Nx, Ny, Nz = 10, 10, 25 
Ne = Nx * Ny * Nz
disp = .1
nFrames = 20
workdir = "D:\donnees_pyth/workdir/"
label = "cuboidTest"
elType = "CPE4"
cpus = 1
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

if compart:
      E  = 210000. * np.ones(Ne) # Young's modulus 210 GPa
      nu = .3 * np.ones(Ne) # Poisson's ratio
      n = 200. * np.ones(Ne)
      Ssat = 2100.* np.ones(Ne)
      sy_mean = 200
      ray_param = sy_mean/1.253314
      sy = np.random.rayleigh(ray_param, Ne)
      labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy))]
      material = [materials.Bilinear(labels = labels[i],E = E[i], nu = nu[i], n = n[i],sy = sy[i], Ssat = Ssat[i]) for i in xrange(Ne)]
    
m = CuboidTest(lx =lx, ly = ly, lz = lz, Nx = Nx, Ny = Ny, Nz = Nz, abqlauncher = abqlauncher, label = label, workdir = workdir, cpus =cpus, material = material, compart = compart, disp = disp, elType = elType, is_3D = True)
    
m.MakeInp()
m.Run()
m.MakePostProc()
    
m.RunPostProc()
    
# Plotting results
if m.outputs['completed']:


  # History Outputs
  disp =  np.array(m.outputs['history']['disp'].values()[0].data[0])
  force =  np.array(np.array(m.outputs['history']['force'].values()).sum().data[0])
  volume = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[0])
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
  
  """
  # Field Outputs
   
  def field_func(outputs, step):
    
    A function that defines the scalar field you want to plot
  
    epsilon = np.array(outputs['field']['LE'][step].get_component(22).data)
    return (epsilon - max_strain) / max_strain
  
  def plot_mesh(ax, mesh, outputs, step, field_func =None, cbar = True, cbar_label = 'Z', cbar_orientation = 'horizontal', disp = True):
 
    A function that plots the deformed mesh with a given field on it.

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
   """