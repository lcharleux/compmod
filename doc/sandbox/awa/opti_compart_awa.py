"""
Compartimented optimization

LC 13/05/2015 
"""
import numpy as np
import matplotlib.pyplot as plt
import compmod, abapy, platform
from scipy import optimize, interpolate

#-------------------------------------------------------------------------------
# FUNCTIONS AND CLASSES
def Saint_Venant(X):
  sy_mean = X[0]
  n = X[1]
  s_sat = X[2]
  cell= lambda eps, sy: compmod.rheology.Bilinear(eps, E = settings['E'], 
                                                       sigmay = sy, 
                                                       n = n, 
                                                       sigma_sat = s_sat)
  grid = np.linspace(0., s_sat*2, 500)
  ray = compmod.distributions.Rayleigh(sy_mean)
  sv = compmod.rheology.SaintVenant(eps, cell, grid, ray)
  return sv
  
def cost_function_sv(X):
  sv = Saint_Venant(X)
  sigma_sv = sv.sigma()
  err = ((sigma_sv - sigma_exp_grid)**2).sum()
  return err

def compart_test(X):
  '''
  Compartmentalized cost function.
  '''
  sigmay_mean = X[0]
  n = X[1]
  sigma_sat = X[2]
  if compart :
    Ec  = settings["E"]  * np.ones(Ne) * k # Young's modulus
    nuc = settings["nu"] * np.ones(Ne)    # Poisson's ratio
    nc  = settings["n"]  * np.ones(Ne) * k  # Hardening slope
    sigma_satc = sigma_sat * np.ones(Ne) * k
    syc = compmod.distributions.Rayleigh(sigmay_mean).rvs(Ne) * k
    labels = ['mat_{0}'.format(i+1) for i in xrange(len(syc))]
    material = [abapy.materials.Bilinear(labels = labels[i], E = Ec[i], nu = nuc[i], sy = syc[i], n = nc[i], Ssat = sigma_satc[i]) for i in xrange(Ne)]
  else:
    labels = 'SAMPLE_MAT'
    material = abapy.materials.Bilinear(labels = labels, E = E * k, nu = nu, sy = sigmay_mean * k, n = n * k , Ssat = sigma_sat * k)
  m = compmod.models.CuboidTest(lx = lx, ly = ly, lz = lz,
      Nx = Nx, Ny = Ny, Nz = Nz,
      abqlauncher = abqlauncher, 
      label = label, 
      workdir = workdir, 
      cpus = cpus, 
      material = material, 
      compart = compart, 
      disp = disp, 
      elType = elType, 
      lateralbc = lateralbc,
      is_3D = True,
      export_fields = False)
  if run_simulation:
    m.MakeInp()
    m.Run()
    m.PostProc()
  else:
    m.LoadResults()
  if m.outputs['completed']:
    disp_sim =  np.array(m.outputs['history']['disp'].values()[0].data[0])
    force_sim =  np.array(np.array(m.outputs['history']['force'].values()).sum().data[0])
    volume_sim = np.array(np.array(m.outputs['history']['volume'].values()).sum().data[0])
    length_sim = ly + disp_sim
    surface_sim = volume_sim / length_sim
    logstrain_sim = np.log10(1. + disp_sim / ly)
    linstrain_sim = disp_sim / ly
    eps_sim = linstrain_sim
    sigma_sim = force_sim / surface_sim / k
    return interpolate.interp1d(eps_sim, sigma_sim)(eps)
    
def cost_function_compart(X):
  sigma_sim_grid = compart_test(X)
  err = ((sigma_sim_grid - sigma_exp_grid)**2).sum()
  return err
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# SETTINGS
settings = {}
settings["exp_path"]     = "cuivre_cufe2p_ANR.txt"  # Path to exp file
settings["results_path"] = "results.txt"           # Path to results file
settings["E"]            = 123.e9                  # Young's modulus
settings["nu"]           = 0.3                     # Poisson's ratio
settings["sy_mean"]      = 500.e6                  # Initial mean yield stress
settings["n"]            = 20.e6                  # Initial hardening
settings["s_sat"]        = 2000.e6                 # Saturation stress 
settings["eps_lim"]      = [1.e-3, .1]             # Limits of the exp strain   
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# EXPERIMENTAL DATA
results = open(settings['results_path'], "wb")
results.write("Method\t sy_mean[MPa]\tn[MPa]\tsigma_sat[Mpa]\n")
data = abapy.misc.read_file(settings["exp_path"])
eps_exp = data[0]          # Experimental strain 
sigma_exp = data[1]        # Experimental stress in Pa

# UNIT CORRECTIONS : stress in Pa and strains in the right unit
sigma_exp *= 1.e6
eps_exp   *= 1.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# EXPERIMENTAL DATA PLOT
fig = plt.figure(0)
plt.clf()
plt.plot(eps_exp,  sigma_exp, "r-", linewidth = 2., label = "Experiment")
plt.grid()
plt.xlabel("Strain $\epsilon$")
plt.ylabel("Stress $\sigma$")
plt.legend(loc = "lower right")
plt.savefig(settings["exp_path"]+"_exp.pdf")
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# SAINT-VENANT PRE OPTIMISATION
# Starting parameters
# Grids definition
print "# Starting Saint-Venant pre optimization"
eps = np.linspace(settings["eps_lim"][0], settings["eps_lim"][1], 100) # Strain grid
sigma_exp_grid = interpolate.interp1d(eps_exp, sigma_exp)(eps)
X0 = [settings["sy_mean"], settings["n"], settings['s_sat']]
sol = optimize.minimize(cost_function_sv, X0, 
    method = "Nelder-Mead", 
    options = {"maxfev": 100, "maxiter":100})
sv = Saint_Venant(sol.x) 
print sol
fig = plt.figure(0)
plt.clf()
plt.plot(eps,  sigma_exp_grid, "r-", linewidth = 2., label = "Experiment")
plt.plot(eps, sv.sigma(), "b-", label = "Saint-Venant")
plt.grid()
plt.xlabel("Strain $\epsilon$")
plt.ylabel("Stress $\sigma$")
plt.legend(loc = "lower right")
plt.savefig(settings["exp_path"]+"_sv.pdf")
# New parameters from Saint-Venant
sigmay_mean = sol.x[0]
n = sol.x[1]
sigma_sat = sol.x[2]
results.write("Saint-Venant\t {0:.2f}\t{1:.2f}\t{2:.2f}\n".format(sol.x[0]*1.e-6, sol.x[1]*1.e-6,sol.x[2]*1.e-6))
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# COMPARTMENTALIZED OPTIMIZATION (start parameters = Saint Venant Solution).
lx, ly, lz = 1., 1., 1.
Nx, Ny, Nz = 20, 20, 20 
Ne = Nx * Ny * Nz
disp = settings["eps_lim"][1] * ly * 1.2
nFrames = 20
workdir = "workdir/"
label = "cuboidTest"
elType = "C3D8"
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
k = 1.e-6 # A factor to keep values reasonnably low in Abaqus, no impact on solution
X0 = sol.x

sol_compart = optimize.minimize(cost_function_compart, X0, 
    method = "Nelder-Mead", 
    options = {"maxfev": 100, "maxiter":100})
results.write("Compart\t {0:.2f}\t{1:.2f}\t{2:.2f}\n".format(sol_compart.x[0]*1.e-6, sol_compart.x[1]*1.e-6,sol_compart.x[2]*1.e-6))
results.close()

sigma_sim_grid = compart_test(sol_compart.x)  

# PLOTING
fig = plt.figure(0)
plt.clf()
plt.plot(eps,  sigma_exp_grid, "r-", linewidth = 2., label = "Experiment")
plt.plot(eps, sv.sigma(), "b-", label = "Saint-Venant")
plt.plot(eps, sigma_sim_grid, "g-", label = "Compart.")

plt.grid()
plt.xlabel("Strain $\epsilon$")
plt.ylabel("Stress $\sigma$")
plt.legend(loc = "lower right")
plt.savefig(settings["exp_path"]+"_sv-cm.pdf")
