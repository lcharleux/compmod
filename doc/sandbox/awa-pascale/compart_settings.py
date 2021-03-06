from compart_classes import Tensile_Test
import platform

# GENERAL SETTINGS
settings = {}
settings['compart']   = True # True: compartimented, False: homogeneous
settings['is_3D']     = True  # leave True....
# MATERIAL
settings["material_type"]  = "Bilinear" # "Bilinear" or "Hollomon"
settings["E"]         = 102.e9 # [Pa]
settings["sy_mean"]   = 380.e6 # [Pa] (only for bilinear)
settings["nu"]        = 0.3    # Poisson's ratio
settings["sigma_sat"] = 450.e6 # [Pa] 
settings["n_hol"]     = 0.1     # Hollomon hardening
settings["n_bil"]     = 2.e9    # [Pa] Bilinear hardening

# SAMPLE DIMENSIONS
settings['lx']        = 1.
settings['ly']        = 1. # test direction 
settings['lz']        = 1.
# MESH
settings['Nx']        = 10
settings['Ny']        = 10
settings['Nz']        = 10
settings['elType']    = "CPS4" # Only for 2D 
# LOADING
settings['disp']      = .061
# NUMBER OF FRAMES
settings['nFrames']   = 100
# NAMES & DIRECTORIES
settings['workdir']   = "workdir/" # Simulation work directory
settings['label']     = "tensile_test"
settings['expdir']    = "experiments/" # Experiments directory 
settings['experiment']= "Ess4_bis_Titane_0.csv" 
# EXPERIMENTS
settings['exp_stress_factor'] = 1.e6 # A factor to adjust units between exp and simulations
settings['exp_strain_factor'] = 1.   # A factor to adjust units between exp and simulations
settings["eps_lim_min"]           = 5.e-4  # Limits of the exp strain
settings["eps_lim_max"]           = .06  # Limits of the exp strain   
# NUMBER OF CORES/CPUS
settings['cpus']      = 1          #Number of CPUs/cores
# OPTIMIZATION SETTINGS
settings["max_number_of_simulations"] = 8 # 4 or more, 100 is ok
settings["max_number_of_iterations"] = 3 # Number of optimization iterations.
# ABAQUS PATH SETTINGS
node = platform.node()
if node ==  'lcharleux-HP':      
  settings['abqlauncher']   = "/opt/abaqus/scratch/Commands/abaqus" 
if node ==  'serv2-ms-symme': 
  settings['abqlauncher']   = "/opt/abaqus/Commands/abaqus"
if node ==  'epua-pd47': 
  settings['abqlauncher']   = "C:/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe" 
if node ==  'epua-pd45': 
  settings['abqlauncher']   = "C:\SIMULIA/Abaqus/Commands/abaqus"
if node ==  'SERV3-MS-SYMME': 
  settings['abqlauncher']   = "C:/Program Files (x86)/SIMULIA/Abaqus/6.11-2/exec/abq6112.exe"
if node == "administrateur": # Awa
  settings['abqlauncher']   = "C:/SIMULIA/Abaqus/6.10-1/exec/abq6101"
if node == 'epua-pb211': # Pascale
  settings['abqlauncher'] = "C:/SIMULIA/Abaqus/Commands/abaqus"

if "abqlauncher" not in settings.keys(): 
  print "Warning: abqlauncher not defined"

if __name__ == '__main__':
  Tensile_Test(settings)

