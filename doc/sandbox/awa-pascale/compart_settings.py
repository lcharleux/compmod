from compart_classes import Tensile_Test
import platform

# GENERAL SETTINGS
settings = {}
settings['compart']   = True # True: compartimented, False: homogeneous
settings['is_3D']     = True  # leave True....
# MATERIAL
settings["material_type"]  = "Hollomon" # "Bilinear" or "Hollomon"
settings["E"]         = 210.e9 # [Pa]
settings["sy_mean"]   = 400.e6 # [Pa] (only for bilinear)
settings["nu"]        = 0.3    # Poisson's ratio
settings["sigma_sat"] = 600.e6 # [Pa] 
settings["n_hol"]     = 0.1     # Hollomon hardening
settings["n_bil"]     = 1.e9    # [Pa] Bilinear hardening

# SAMPLE DIMENSIONS
settings['lx']        = 1.
settings['ly']        = 1. # test direction 
settings['lz']        = 1.
# MESH
settings['Nx']        = 3
settings['Ny']        = 3
settings['Nz']        = 3
settings['elType']    = "CPS4" # Only for 2D 
# LOADING
settings['disp']      = .2
# NUMBER OF FRAMES
settings['nFrames']   = 100
# NAMES & DIRECTORIES
settings['workdir']   = "workdir/"
settings['label']     = "tensile_test" 
# NUMBER OF CORES/CPUS
settings['cpus']      = 1          #Number of CPUs/cores


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

