from compart_classes import Tensile_Test
import platform

# GENERAL SETTINGS
settings = {}
settings['compart'] = True
settings['is_3D']   = True
settings['lx']      = 5.
settings['ly']      = 5. # test direction 
settings['lz']      = 5.
settings['Nx']      = 3
settings['Ny']      = 3
settings['Nz']      = 3
settings['disp']    = 1.2
settings['nFrames'] = 50
settings['workdir'] = "workdir/"
settings['label']   = "tensile_test"
settings['elType']  = "CPS4"
settings['cpus']    = 1
settings["E"]       = 1.
settings["sy_mean"]     = .001
settings["nu"]          = 0.3
settings["sigma_sat"]   = .005
settings["n"]           = 0.05


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


Tensile_Test(settings)

