from abapy import materials
from abapy.misc import load
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import pickle, copy, platform, compmod, os




def Tensile_Test(settings):

  # MATERIALS CREATION
  Ne = settings['Nx'] * settings['Ny'] 
  if settings['is_3D']: Ne *= settings['Nz']
  if settings['compart']:
    E         = settings["E"] * np.ones(Ne) # Young's modulus
    nu        = settings["nu"]        * np.ones(Ne) # Poisson's ratio
    sy_mean   = settings["sy_mean"]   * np.ones(Ne)
    sigma_sat = settings["sigma_sat"] * np.ones(Ne)
    n         = settings["n"]         * np.ones(Ne)
    sy = compmod.distributions.Rayleigh(settings["sy_mean"]).rvs(Ne)
    labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy_mean))]
    settings['material'] = [materials.Bilinear(labels = labels[i], 
                                   E = E[i], nu = nu[i], Ssat = sigma_sat[i], 
                                   n=n[i], sy = sy[i]) for i in xrange(Ne)]
  else:
    labels = 'SAMPLE_MAT'
    settings['material'] = materials.Bilinear(labels = labels, 
                                  E = E, nu = nu, sy = sy_mean, Ssat = sigma_sat,
                                  n = n)

         
  m = compmod.models.CuboidTest(**settings)
  m.MakeInp()
  m.Run()
  m.MakePostProc()
  m.RunPostProc()
  m.LoadResults()
  # Plotting results
  if m.outputs['completed']:
    

    # History Outputs
    disp =  np.array(m.outputs['history']['disp'].values()[0].data[0])
    force =  np.array(np.array(m.outputs['history']['force'].values()).sum().data[0])
    volume = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[0])
    length = settings['ly'] + disp
    surface = volume / length
    logstrain = np.log10(1. + disp / settings['ly'])
    linstrain = disp/ settings['ly']
    strain = linstrain
    stress = force / surface 
    output = {}
    output["force"] = force 
    output["disp"]  = disp 
    output["stress"] = stress
    output["strain"] = strain
    df = pd.DataFrame(output)
    df.to_csv("{0}{1}.csv".format(settings["workdir"], settings["label"]))
    df.to_excel("{0}{1}.xls".format(settings["workdir"], settings["label"]))
    
    
   
  
