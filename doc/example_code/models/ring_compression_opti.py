# SOME OPTIMIZATION WITH RING COMPRESSION

from compmod.models import RingCompression
import matplotlib.pyplot as plt
import numpy as np
import pickle, copy

#FIXED PAREMETERS
inner_radius, outer_radius = 95, 100
Nt, Nr = 10, 3 
disp = 50.
nFrames = 50
sy = 200.e6
E = 200.e9
nu = .3
n = .3
workdir = "workdir/"
label = "ringCompressionOpti"
elType = "CPS4"



#MODEL DEFINITION
material = Hollomon(
  labels = "SAMPLE_MAT",
  E = E, nu = nu,
  sy = sy, n = n)
m = RingCompression(
  material = material , 
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




