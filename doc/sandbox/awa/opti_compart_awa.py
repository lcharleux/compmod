"""
Compartimented optimization

LC 13/05/2015 
"""
import numpy as np
import matplotlib.pyplot as plt
import compmod, abapy, platform
from scipy import optimize, interpolate

# EXPERIMENTAL DATA
path = "cuivre_cufe2p_ANR.txt"
results = open("results.txt", "wb")
results.write("Method\t sy_mean[MPa]\tn[MPa]\tsigma_sat[Mpa]\n")
data = abapy.misc.read_file(path)
eps_exp = data[0]          # Experimental strain 
sigma_exp = data[1] * 1.e6 # Experimental stress in Pa

# FIXED  PARAMETERS
E = 123.e9 
nu = .3 

# SAINT-VENANT PRE OPTIMISATION
# Starting parameters
sigmay_mean = 400.e6
n = 200.e6
sigma_sat = 5000.e6
# Grids definition
eps = np.linspace(1.e-3, 0.1, 100) # Strain grid
sigma_exp_grid = interpolate.interp1d(eps_exp, sigma_exp)(eps)

def Saint_Venant(X):
  sigmay_mean = X[0]
  n = X[1]
  sigma_sat = X[2]
  cell= lambda eps, sy: compmod.rheology.Bilinear(eps, E, sy, n, sigma_sat)
  grid = np.linspace(0., 5000.e6, 1000)
  ray = compmod.distributions.Rayleigh(sigmay_mean)
  sv = compmod.rheology.SaintVenant(eps, cell, grid, ray)
  return sv
  
def cost_function(X):
  sv = Saint_Venant(X)
  sigma_sv = sv.sigma()
  err = ((sigma_sv - sigma_exp_grid)**2).sum()
  return err
  
X0 = [sigmay_mean, n, sigma_sat]
sol = optimize.minimize(cost_function, X0, 
    method = "Nelder-Mead", 
    options = {"maxfev": 5000, "maxiter":5000})
sv = Saint_Venant(sol.x) 

# New parameters from Saint-Venant
sigmay_mean = sol.x[0]
n = sol.x[1]
sigma_sat = sol.x[2]
results.write("Saint-Venant\t {0:.2f}\t{1:.2f}\t{2:.2f}\n".format(sol.x[0]*1.e-6, sol.x[1]*1.e-6,sol.x[2]*1.e-6))


# COMPARTMENTALIZED OPTIMIZATION (start parameters = Saint Venant Solution).

#PARAMETERS
lx, ly, lz = 1., 1., 1.
Nx, Ny, Nz = 5, 5, 5 
Ne = Nx * Ny * Nz
disp = .12
nFrames = 20
workdir = "workdir/"
label = "cuboidTest_homo"
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

k = 1.e-6
E = E
nu = nu

def compart_test(X):
  """
  Compartmentalized cost function.
  """
  sigmay_mean = X[0]
  n = X[1]
  sigma_sat = X[2]
  if compart :
    Ec  = E * np.ones(Ne) * k # Young's modulus
    nuc = nu * np.ones(Ne)    # Poisson's ratio
    nc = n * np.ones(Ne) * k  # Hardening slope
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

X0 = sol.x

sol_compart = optimize.minimize(cost_function_compart, X0, 
    method = "Nelder-Mead", 
    options = {"maxfev": 1, "maxiter":1})
results.write("Compart\t {0:.2f}\t{1:.2f}\t{2:.2f}\n".format(sol_compart.x[0]*1.e-6, sol_compart.x[1]*1.e-6,sol_compart.x[2]*1.e-6))
results.close()

sigma_sim_grid = compart_test(sol_compart.x)  
 
# PLOTING
fig = plt.figure(0)
plt.plot(eps_exp,  sigma_exp, "r-", linewidth = 2., label = "Experiment")
plt.plot(eps, sv.sigma(), "b-", label = "Saint-Venant")
plt.plot(eps, sigma_sim_grid, "g-", label = "Compart.")

plt.grid()
plt.xlabel("Strain $\epsilon$")
plt.ylabel("Stress $\sigma$")
plt.legend(loc = "lower right")
plt.show()
