# SOME OPTIMIZATION WITH RING COMPRESSION

from abapy.materials import Hollomon
from compmod.models import RingCompression
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy
import platform
node = platform.node()


#FIXED PAREMETERS
inner_radius, outer_radius = 95, 100
Nt, Nr = 10, 3 
displacement = 50.
nFrames = 50
sy = 200.e6
E = 200.e9
nu = .3
n = .3
if node ==  'lcharleux':      
  abqlauncher   = '/opt/Abaqus/6.9/Commands/abaqus' # Local machine configuration
  workdir = "workdir/"
if node ==  'epua-pd47': 
  abqlauncher   = 'C:/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe' # Local machine configuration
  workdir = "D:/Simulations/Dossier_travail_Abaqus/"
label = "ringCompressionOpti"
elType = "CPS4"


def RunSimulation(sy, n):
  """
  Runs a simulation for a given couple (sy, n) and returns the (disp, force) couple.
  """
  #MODEL DEFINITION
  material = Hollomon(
    labels = "SAMPLE_MAT",
    E = E, nu = nu,
    sy = sy, n = n)
  m = RingCompression(
    material = material , 
    inner_radius = inner_radius, 
    outer_radius = outer_radius, 
    disp = displacement, 
    nFrames = nFrames, 
    Nr = Nr, 
    Nt = Nt, 
    workdir = workdir,
    label = label, 
    elType = elType,
    abqlauncher = abqlauncher)

  # SIMULATION
  m.MakeMesh()
  m.MakeInp()
  m.Run()
  m.PostProc()
  outputs = m.outputs
  force = -2. * outputs['history']['force']
  disp = -2. * outputs['history']['disp']
  return disp, force

RunSimulation(200.e6, 0.3)




