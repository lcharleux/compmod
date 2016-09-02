from abapy import materials
from abapy.misc import load
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import pickle, copy, platform, compmod, os
from scipy import optimize, interpolate



def Tensile_Test(settings):
  args = settings.copy()
  # MATERIALS CREATION
  Ne = settings['Nx'] * settings['Ny'] 
  if settings['is_3D']: Ne *= settings['Nz']
  if settings['compart']:
    E         = settings["E"] * np.ones(Ne) # Young's modulus
    nu        = settings["nu"]        * np.ones(Ne) # Poisson's ratio
    sy_mean   = settings["sy_mean"]   * np.ones(Ne)
    sigma_sat = settings["sigma_sat"] * np.ones(Ne)
    n_hol     = settings["n_hol"]     * np.ones(Ne)
    n_bil     = settings["n_bil"]     * np.ones(Ne)
    sy = compmod.distributions.Rayleigh(settings["sy_mean"]).rvs(Ne)
    labels = ['mat_{0}'.format(i+1) for i in xrange(len(sy_mean))]
    if args['material_type']  == "Bilinear":
      args['material'] = [materials.Bilinear(labels = labels[i], 
                                     E = E[i], nu = nu[i], Ssat = sigma_sat[i], 
                                     n=n_bil[i], sy = sy[i]) for i in xrange(Ne)]
    if args['material_type']  == "Hollomon":
      args['material'] = [materials.Hollomon(labels = labels[i], 
                                     E = E[i], nu = nu[i], n=n_hol[i], 
                                     sy = sy[i]) for i in xrange(Ne)]
  else:
    labels = 'SAMPLE_MAT'
    if args['material_type']  == "Bilinear":
      args['material'] = materials.Bilinear(labels = labels, 
                                    E = settings["E"], 
                                    nu = settings["nu"], 
                                    sy = settings["sy_mean"], 
                                    Ssat = settings["sigma_sat"],
                                    n = settings["n_bil"])
    if args['material_type']  == "Hollomon":
      args['material'] = materials.Hollomon(labels = labels, 
                                     E = settings["E"], 
                                    nu = settings["nu"], 
                                    sy = settings["sy_mean"], 
                                    n = settings["n_hol"])
         
  m = compmod.models.CuboidTest(**args)
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
    output["volume"] = volume
    output["length"] = length
    df = pd.DataFrame(output)
    df.to_csv("{0}{1}.csv".format(settings["workdir"], settings["label"]), index = False)
    df.to_excel("{0}{1}.xls".format(settings["workdir"], settings["label"]), index = False)
    inputs = pd.DataFrame(settings, index = [0])
    inputs.transpose().to_csv("{0}{1}_inputs.csv".format(settings["workdir"], settings["label"]))
    

   
class Optimize(object):
  def __init__(self, settings):
    settings = settings.copy()
    settings['compart']   = True # True: compartimented, False: homogeneous
    settings["material_type"]  = "Bilinear" # "Bilinear" or "Hollomon"
    self.settings = settings
    self.inputs = []
    self.sim_id = 1
    exp = pd.read_csv(settings["expdir"] + settings["experiment"], 
                    delim_whitespace = True)
    exp.stress *= settings['exp_stress_factor']
    exp.strain *= settings['exp_strain_factor']
    # STRAIN GRID FOR LEAST SQUARE OPTIMIZATION
    self.strain_grid = np.linspace(
             max(settings["eps_lim_min"], exp.strain.min()), 
             min(settings["eps_lim_max"], exp.strain.max()), 
             100) # Strain grid
    self.sigma_exp_grid = interpolate.interp1d(exp.strain, exp.stress)(self.strain_grid)
    
  def Cost_Function(self, X):
    settings = self.settings.copy()
    settings["sy_mean"]   = X[0] # [Pa] (only for bilinear)
    settings["n_bil"]     = X[1] # [Pa] Bilinear hardening
    settings["sigma_sat"] = X[2] # [Pa] 
    # SIMULATION LABEL RENAMING
    settings["label"] += "_{0}".format(self.sim_id)
    Tensile_Test(settings)
    sim = pd.read_csv(settings["workdir"] + settings["label"] + ".csv")
    strain_grid = self.strain_grid
    sigma_sim_grid = interpolate.interp1d(sim.strain, sim.stress)(strain_grid)
    err = ((sigma_sim_grid - self.sigma_exp_grid)**2).sum()/len(strain_grid)
    self.inputs.append([self.sim_id] + list(X))
    self.sim_id +=1
    return err
  
  def run(self):  
    # EXPERIMENTAL DATA
    settings = self.settings
    # OPTIMiZATION START POINT
    X0 = np.array([settings["sy_mean"], settings["n_bil"], settings['sigma_sat']])
    sol_compart = optimize.minimize(self.Cost_Function, X0, 
      method = "Nelder-Mead", 
      options = {"maxfev": settings["max_number_of_simulations"]})
    self.sol = sol_compart  
    inputs = np.array(self.inputs).transpose()
    df = pd.DataFrame({"sim_id": inputs[0],
                       "sy_mean":inputs[1],
                       "n_bil":  inputs[2],
                       "sigma_sat":inputs[3],})
    df.to_csv("{0}{1}_opti_results.csv".format(settings["workdir"], settings["label"]), index = False)                   
    
