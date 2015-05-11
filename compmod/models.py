'''
Models
=========
'''
# A repository for compartmentalized models compatible simulation models.
from abapy.mesh import RegularQuadMesh
from abapy.materials import VonMises, Hollomon
from abapy.misc import load
import numpy as np
import os, time, subprocess, pickle, copy


class Simulation(object):
  """
  :param abqlauncher: path to the abaqus executable.
  :type abqlauncher: string
  :param material: material instance from `abapy.materials`
  :param label: label of the simulation (default: 'simulation')
  :type label: string
  :param workdir: path to the simulation work directory (default: '')
  :type workdir: string
  :param compart: indicated if the simulation homogeneous or compartimented (default: False)
  :type compart: boolean 
  :param nFrames: number or frames per step (default: 100)
  :type nFrames: integer
  :param elType: element type (default: 'CPS4')
  :type elType: string
  :param is_3D: indicates if the model is 2D or 3D (default: False)
  :type is_3D: boolean 
  :param cpus: number of CPUs to use (default: 1)
  :type compart: integer
  """
  def __init__(self, **kwargs):
    defaultArgs = {"abqlauncher":None, "material": VonMises(), "label": "simulation",  "workdir": "", "compart":False, "nFrames": 100, "elType": "CPS4", "is_3D": False, "cpus" : 1}
    for key, value in defaultArgs.iteritems(): setattr(self, key, value)
    for key, value in kwargs.iteritems(): setattr(self, key, value)
    if self.is_3D: 
      self.elType = 'C3D8'
    
  def Run(self, deleteOldFiles = True):
    '''
    Runs the simulation.
    
    :param deleteOldFiles: indicates if existing simulation files are deleted before the simulation starts.
    :type deleteOlfFiles: boolean
    '''
    if deleteOldFiles: self.DeleteOldFiles()
    t0 = time.time()
    print '< Running simulation {0} in Abaqus>'.format(self.label) 
    p = subprocess.Popen( '{0} job={1} input={1}.inp cpus={2} interactive'.format(self.abqlauncher, self.label, self.cpus), cwd = self.workdir, shell=True, stdout = subprocess.PIPE)
    trash = p.communicate()
    print trash[0]
    t1 = time.time()
    self.duration = t1 - t0
    print '< Ran {0} in Abaqus: duration {1:.2f}s>'.format(self.label, t1 - t0)   
  
  def DeleteOldFiles(self):
    """
    Deletes old job files.
    """
    suffixes = ['.odb', '.lck', '.log', '.dat', '.com', '.sim', '.sta', '.prt', '.msg']
    for s in suffixes:
      path = self.workdir + self.label + s
      try:
        os.remove(path)
      except OSError:
        pass
  
  def RunPostProc(self):
    """
    Runs the post processing script.
    """
    t0 = time.time()
    p = subprocess.Popen( "{0} viewer -noGUI {1}".format(self.abqlauncher, self.label + '_abqpostproc.py'), cwd = self.workdir,stdout = subprocess.PIPE, shell=True)
    trash = p.communicate()
    print trash[0]
    t1 = time.time()
    print '< Post Processed {0} in Abaqus: duration {1:.2f}s>'.format(self.label, t1 - t0) 
    self.outputs = load(self.workdir + self.label + ".pckl")
  
  def PostProc(self):
    """
    Makes the post proc script and runs it.
    """
    self.MakePostProc()
    self.RunPostProc()
  
  def LoadResults(self):
    """
    Loads the results from a pickle file.
    """
    self.outputs = load(self.workdir + self.label + ".pckl")
    
    
class CuboidTest(Simulation):
  """
  Performs a tensile or compressive test on an cuboid rectangular cuboid. 

  
  :param lx: length of the box along the x axis (default = 1.)
  :type lx: float
  :param ly: length of the box along the y axis (default = 1.)
  :type ly: float
  :param lz: length of the box along the z axis (default = 1.). Only used in 3D simulations.
  :type lz: float
  :param Nx: number of elements along the x direction.
  :type Nx: int
  :param Ny: number of elements along the y direction.
  :type Ny: int
  :param Nz: number of elements along the z direction.
  :type Nz: int
  :param disp: imposed displacement along the y direction (default = .25)
  :type disp: float 
  :param export_fields: indicates if the field outputs are exported (default = True). Can be set to False to speed up post processing.
  :type export: boolean
  :param lateral_bc: indicates the type of lateral boundary conditions to be used.
  :type lateral_bc: dict
  {0}
  
    
  CuboidTest with microstructure generated using Voronoi cells : 
   
  * Source: :download:`cuboidTest_voronoi <example_code/models/cuboidTest_voronoi.py>`.
  * VTK output: :download:`cuboidTest_voronoi <example_code/models/cuboidTest_voronoi.vtk>`.
 
  """
  __doc__ = __doc__.format(Simulation.__doc__)
  
  
  def __init__(self, **kwargs):
    
    
    defaultArgs = {"Nx":10, "Ny":10, "Nz":10, "lx":1., "ly":1., "lz":1., "disp":.25, "export_fields": True, "lateralbc":{}}
    for key, value in defaultArgs.iteritems(): setattr(self, key, value)
    for key, value in kwargs.iteritems(): setattr(self, key, value)
    super(CuboidTest, self).__init__(**kwargs)
  
  def MakeMesh(self):
    """
    Builds the mesh.
    """
    Nx , Ny, Nz = self.Nx, self.Ny, self.Nz
    lx, ly, lz = self.lx, self.ly, self.lz
    elType = self.elType
    if self.is_3D:
      Ne = Nx * Ny * Nz
    else:
      Ne = Nx * Ny
    m = RegularQuadMesh(Nx, Ny, l1= lx, l2 = ly, name = elType)
    m.add_set(label = "AllElements", elements = m.labels)
    nsets = copy.copy(m.nodes.sets) 
    if self.is_3D: 
       m = m.extrude(N = Nz, l = lz)
       m.nodes.sets['bottomleft'] = nsets['bottomleft']
       m.nodes.sets['bottomright'] = nsets['bottomright']
    self.mesh = m
    
  def MakeInp(self):
    """
    Writes the Abaqus INP file in the workdir.
    """
    pattern = """**----------------------------------
**DISTRIBUTED MECHANICAL PROPERTIES
**----------------------------------
**HEADER
*Preprint, echo=NO, model=NO, history=NO, contact=NO
**----------------------------------
** PART "pSAMPLE" DEFINITION
*Part, name = pSample
#MESH
#SECTIONS
#LATERALBC
*End part
**----------------------------------
** ASSEMBLY
*Assembly, name = Assembly
*Instance, name=iSample, part=pSample
*End Instance
*End Assembly
**----------------------------------
** MATERIALS
#MATERIALS
**----------------------------------
** STEPS
*Step, Name=Loading0, Nlgeom=YES, Inc=1000000
*Static
#FRAME_DURATION, 1, 1e-08, #FRAME_DURATION
** BOUNDARY CONDITIONS
*Boundary
iSample.Bottom, 2, 2
iSample.BottomLeft, 1, 1#3DBOUNDARY
iSample.Top,    2, 2, #DISP
** RESTART OPTIONS 
*Restart, write, frequency=0
** FIELD OUTPUTS
*Output, field, frequency=999999
*Node Output
U
*Element Output, directions=YES
E, PE, EE, PEEQ, S
** HYSTORY OUTPUTS
*Output, history
*Energy Output
ALLPD, ALLSE, ALLWK
*Node Output, nset=iSample.Top
RF2
*Node Output, nset=iSample.TopLeft
U2
*Node Output, nset=iSample.Left
COOR1
*Node Output, nset=iSample.Right
COOR1
*Element Output, elset=iSample.allElements, directions=NO
EVOL
*End Step
  """
    Nx , Ny, Nz = self.Nx, self.Ny, self.Nz
    lx, ly, lz = self.lx, self.ly, self.lz
    elType = self.elType
    material = self.material
    disp = self.disp
    nFrames = self.nFrames
    if self.is_3D:
      Ne = Nx * Ny * Nz
    else:
      Ne = Nx * Ny
    sections = ""
    matinp = ""
    if self.compart:
      section_pattern = "*Solid Section, elset=Elset{0}, material={1}\n*Elset, Elset=Elset{0}\n{0},\n"
      labels = [mat.labels[0] for mat in material]
      for i in xrange(Ne):
        sections += section_pattern.format(i+1, labels[i]) 
        matinp += material[i].dump2inp() + '\n'
    else:
      section_pattern = "*SOLID SECTION, ELSET = ALLELEMENTS, MATERIAL = {0}\n{1}"  
      label = material.labels[0]
      sections = section_pattern.format(label, self.lz)
      matinp = material.dump2inp() 
    if hasattr(self, "mesh") == False:
      self.MakeMesh()
    m = self.mesh
    lateralbc = ""
    if len(self.lateralbc.keys()) != 0:
      lateralbc += "*EQUATION\n"
      lateralbc_keys = self.lateralbc.keys()
      for lbck in lateralbc_keys:
        if lbck == "right": 
          direction = 1
          nset = m.nodes.sets['right'] 
        if lbck == "left": 
          direction = 1
          nset = m.nodes.sets['left']
        if self.lateralbc[lbck] == "pseudohomo":
          for nodelabel in nset[1:]:
            lateralbc += "2\n{0}, 1, 1, {1}, 1, -1\n".format(nodelabel, nset[0])
        """
        if self.lateralbc[lbck] == "periodic":
          left_nodes = m.nodes.sets['left']
          right_nodes = m.nodes.sets['right']
          xl = np.array([])
          side_pairs = []
          top_pair = []
          bottom_pair = []
          
          for nl in left_nodes:
            nlx = 
        """       
    pattern = pattern.replace("#LATERALBC", lateralbc[:-1])  
    pattern = pattern.replace("#MESH", m.dump2inp())
    pattern = pattern.replace("#SECTIONS", sections[:-1])
    pattern = pattern.replace("#MATERIALS", matinp[:-1])
    pattern = pattern.replace("#DISP", str(disp))
    pattern = pattern.replace("#FRAME_DURATION", str(1./nFrames))
    if self.is_3D:
      pattern = pattern.replace("#3DBOUNDARY", "\niSample.BottomLeft, 3, 3\niSample.BottomRight, 3, 3")
    else:  
      pattern = pattern.replace("#3DBOUNDARY", "")
    f = open(self.workdir + self.label + '.inp', 'wb')
    f.write(pattern)
    f.close()
  def MakePostProc(self):
    """
    Makes the post-proc script
    """
    pattern = """# ABQPOSTPROC.PY
# Warning: executable only in abaqus abaqus viewer -noGUI,... not regular python.
import sys
from abapy.postproc import GetFieldOutput_byRpt as gfo
from abapy.postproc import GetVectorFieldOutput_byRpt as gvfo
from abapy.postproc import GetTensorFieldOutput_byRpt as gtfo
from abapy.postproc import GetHistoryOutputByKey as gho
from abapy.postproc import GetMesh
from abapy.indentation import Get_ContactData
from abapy.misc import dump
from odbAccess import openOdb
from abaqusConstants import JOB_STATUS_COMPLETED_SUCCESSFULLY



# Odb opening  
file_name = '#FILE_NAME'
odb = openOdb(file_name + '.odb')
data = {}

# Check job status:
job_status = odb.diagnosticData.jobStatus

if job_status == JOB_STATUS_COMPLETED_SUCCESSFULLY:
  data['completed'] = True"""
    if self.export_fields : pattern += """
  # Field Outputs
  data['field'] = {}
  fo = data['field']
  fo['U'] = [
    gvfo(odb = odb, 
      instance = 'ISAMPLE', 
      step = 0,
      frame = -1,
      original_position = 'NODAL', 
      new_position = 'NODAL', 
      position = 'node',
      field = 'U', 
      delete_report = True)
      ]
 
      
  fo['S'] = [
    gtfo(odb = odb, 
      instance = 'ISAMPLE', 
      step = 0,
      frame = -1,
      original_position = 'INTEGRATION_POINT', 
      new_position = 'NODAL', 
      position = 'node',
      field = 'S', 
      sub_set_type = 'element', 
      delete_report = True),
    ]
   
  fo['LE'] = [
    gtfo(odb = odb, 
      instance = 'ISAMPLE', 
      step = 0,
      frame = -1,
      original_position = 'INTEGRATION_POINT', 
      new_position = 'NODAL', 
      position = 'node',
      field = 'LE', 
      sub_set_type = 'element', 
      delete_report = True),
    ] 
      
  fo['EE'] = [
    gtfo(odb = odb, 
      instance = 'ISAMPLE', 
      step = 0,
      frame = -1,
      original_position = 'INTEGRATION_POINT', 
      new_position = 'NODAL', 
      position = 'node',
      field = 'EE', 
      sub_set_type = 'element', 
      delete_report = True),
    ]     
  
  fo['PE'] = [
    gtfo(odb = odb, 
      instance = 'ISAMPLE', 
      step = 0,
      frame = -1,
      original_position = 'INTEGRATION_POINT', 
      new_position = 'NODAL', 
      position = 'node',
      field = 'PE', 
      sub_set_type = 'element', 
      delete_report = True),
    ] 
  
  fo['PEEQ'] = [
    gfo(odb = odb, 
      instance = 'ISAMPLE', 
      step = 0,
      frame = -1,
      original_position = 'INTEGRATION_POINT', 
      new_position = 'NODAL', 
      position = 'node',
      field = 'PEEQ', 
      sub_set_type = 'element', 
      delete_report = True),
    ] """
    
    pattern += """# History Outputs
  data['history'] = {} 
  ho = data['history']
  ho['disp'] =   gho(odb,'U2')
  ho['force'] =   gho(odb,'RF2')
  ho['allse'] =   gho(odb,'ALLSE').values()[0]
  ho['allpd'] =   gho(odb,'ALLPD').values()[0]
  ho['allwk'] =   gho(odb,'ALLWK').values()[0]
  ho['volume'] =  gho(odb,'EVOL')
  
  # Mesh 
  data['mesh'] = GetMesh(odb, "ISAMPLE")
  
else:
  data['completed'] = False
# Closing and dumping
odb.close()
dump(data, file_name+'.pckl')"""
    pattern = pattern.replace("#FILE_NAME", self.label)
    f = open(self.workdir + self.label + '_abqpostproc.py', 'w')
    f.write(pattern)
    f.close()
  
  
class RingCompression(Simulation):
  """
 let see 2 kind of RingCompression , one homogenous and the second compartmentalized
 
 RingCompression_3D
    
 .. plot:: example_code/models/ring_compression_3D.py
     :include-source:
     
 RingCompression_3D compartimentalized

 .. plot:: example_code/models/ring_compression_3D_compart.py
     :include-source:
  """
  def __init__(self, **kwargs):
    """
    :param inner_radius: inner radius of the ring
    :type inner_radius: float
    :param outer_radius: outer radius of the ring
    :type outer_radius: float
    :param Nr: Number of elements in the radial direction
    :type Nr: int
    :param Nt: Number of elements in the orthoradial direction
    :type Nt: int
    :param Na: Number of elements in the axial direction. Used only if the is_3D option is active.
    :type Na: int
    """
    defaultArgs = {
      "inner_radius": 1., 
      "inner_radius": 2.,  
      "thickness": 1.,
      "Nr":10, 
      "Nt":10,
      "Na":10, 
      "disp": .5,
      "unloading": True,
      "export_fields": True
      }
    for key, value in defaultArgs.iteritems(): setattr(self, key, value)
    for key, value in kwargs.iteritems(): setattr(self, key, value)
    super(RingCompression, self).__init__(**kwargs)
    
  def MakeMesh(self):
    """
    Builds the mesh
    """
    Ri = self.inner_radius
    Ro = self.outer_radius
    thickness = self.thickness
    Nr, Nt, Na = self.Nr, self.Nt, self.Na
    mesh = RegularQuadMesh(Nt, Nr, .25, Ro - Ri, name = self.elType)
    mesh.nodes.add_set_by_func('left_nodes', lambda x, y, z, labels: x == 0.)
    mesh.nodes.add_set_by_func('right_nodes', lambda x, y, z, labels: x == .25)
       
    def function(x, y, z, labels):
      theta = 2 * np.pi * (.25 - x)
      r = y + Ri
      ux = -x + r * np.cos(theta)
      uy = -y + r * np.sin(theta)
      uz = 0. * z
      return ux, uy, uz
    vectorField = mesh.nodes.eval_vectorFunction(function)
    mesh.nodes.apply_displacement(vectorField)
    nodes = mesh.nodes
    for i in xrange(len(nodes.labels)):
      if nodes.x[i] < 0.: 
        nodes.x[i] = 0. 
    
    mesh.add_set('all_elements', mesh.labels)
    mesh.add_set('surface_elements',range( Nt * (Nr-1)+1, Nt*Nr+1  ))
    mesh.add_surface('surface_faces',[ ('surface_elements',3) ])
    if self.is_3D:
       mesh = mesh.extrude(N = Na, l = thickness, mapping = {self.elType: self.elType})  
    self.mesh = mesh
  
  def MakeInp(self):
    pattern = """**----------------------------------
**RING COMPRESSION SIMULATION
**----------------------------------
**HEADER
*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO, CONTACT=NO
**----------------------------------
** SAMPLE DEFINITION
*PART, NAME = P_SAMPLE
#RING_MESH
#SECTIONS
*END PART
**----------------------------------
** INDENTER DEFINITION
**----------------------------------
*PART, NAME = P_PLATE
*END PART
**----------------------------------
** ASSEMBLY
**----------------------------------
*ASSEMBLY, NAME = ASSEMBLY
*INSTANCE, NAME = I_SAMPLE, PART = P_SAMPLE
*END INSTANCE
*INSTANCE, NAME = I_PLATE, PART= P_PLATE
*NODE, NSET=REFNODE
  1, 0., 0., 0.
*SURFACE, TYPE=#SURFTYPE, NAME=SURFACE
  START, #OUTER_RADIUS, #OUTER_RADIUS
  LINE,  0., #OUTER_RADIUS         
*RIGID BODY, REF NODE=REFNODE, ANALYTICAL SURFACE=SURFACE
*END INSTANCE
*END ASSEMBLY
**----------------------------------
** SURFACE INTERACTIONS
**----------------------------------
*SURFACE INTERACTION, NAME = SURF_INT
*FRICTION
0.0,
*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE = HARD
*CONTACT PAIR, INTERACTION = SURF_INT, SUPPLEMENTARY CONSTRAINTS = NO
I_SAMPLE.SURFACE_FACES, I_PLATE.SURFACE
**----------------------------------
** MATERIALS
**----------------------------------
** SAMPLE MATERIAL
#MATERIALS
**----------------------------------
** STEPS
**----------------------------------
*STEP, NAME = LOADING, NLGEOM = YES, INC=1000
*Static
#FRAME_DURATION, 1, 1e-08, #FRAME_DURATION
*BOUNDARY
I_SAMPLE.LEFT_NODES, 1, 1
I_SAMPLE.RIGHT_NODES, 2, 2#3DBOUNDARY
I_PLATE.REFNODE, 2, 2, #DISP
I_PLATE.REFNODE, 1, 1
I_PLATE.REFNODE, 3, 6
*RESTART, WRITE, FREQUENCY = 0
*OUTPUT, FIELD, FREQUENCY = 1
*NODE OUTPUT
COORD, U, 
*ELEMENT OUTPUT, ELSET=I_SAMPLE.ALL_ELEMENTS, DIRECTIONS = YES
LE, EE, PE, PEEQ, S, 
*OUTPUT, HISTORY
*ENERGY OUTPUT
ALLFD, ALLWK
*ENERGY OUTPUT, ELSET=I_SAMPLE.ALL_ELEMENTS
ALLPD
*ENERGY OUTPUT, ELSET=I_SAMPLE.ALL_ELEMENTS
ALLSE
*NODE OUTPUT, NSET=I_PLATE.REFNODE
RF2, U2
*END STEP
"""
    
    if self.unloading:
      pattern += """*STEP, NAME = UNLOADING, NLGEOM = YES, INC=1000
*Static
#FRAME_DURATION, 1, 1e-08, #FRAME_DURATION
*BOUNDARY
I_SAMPLE.LEFT_NODES, 1, 1
I_SAMPLE.RIGHT_NODES, 2, 2
I_PLATE.REFNODE, 2, 2, 0.
I_PLATE.REFNODE, 1, 1
I_PLATE.REFNODE, 3, 6
*RESTART, WRITE, FREQUENCY = 0
*OUTPUT, FIELD, FREQUENCY = 1
*NODE OUTPUT
COORD, U, 
*ELEMENT OUTPUT, ELSET=I_SAMPLE.ALL_ELEMENTS, DIRECTIONS = YES
LE, EE, PE, PEEQ, S, 
*OUTPUT, HISTORY
*ENERGY OUTPUT
ALLFD, ALLWK
*ENERGY OUTPUT, ELSET=I_SAMPLE.ALL_ELEMENTS
ALLPD
*ENERGY OUTPUT, ELSET=I_SAMPLE.ALL_ELEMENTS
ALLSE
*NODE OUTPUT, NSET=I_PLATE.REFNODE
RF2, U2
*END STEP"""
    Nr , Nt, Na = self.Nr, self.Nt, self.Na
    if self.is_3D:
      Ne = Nr * Nt * Na
    else:  
      Ne = Nr * Nt
    material = self.material
    sections = ""
    matinp = ""
    if self.compart:
      section_pattern = "*Solid Section, elset=Elset{0}, material={1}\n{2}\n*Elset, Elset=Elset{0}\n{0},\n"
      labels = [mat.labels[0] for mat in material]
      for i in xrange(Ne):
        sections += section_pattern.format(i+1, labels[i], self.thickness) 
        matinp += material[i].dump2inp() + '\n'
      matinp = matinp[:-1]  
    else:
      section_pattern = "*SOLID SECTION, ELSET = ALL_ELEMENTS, MATERIAL = {0}\n{1}"  
      label = material.labels[0]
      sections = section_pattern.format(label, self.thickness)
      matinp = material.dump2inp() 
    
    if hasattr(self, "mesh") == False:
      self.MakeMesh()    
    
    pattern = pattern.replace('#RING_MESH', self.mesh.dump2inp())
    pattern = pattern.replace('#OUTER_RADIUS', str(self.outer_radius))
    pattern = pattern.replace('#DISP', str(-self.disp))
    pattern = pattern.replace('#FRAME_DURATION', str(1. / self.nFrames))
    pattern = pattern.replace('#SECTIONS', sections)
    pattern = pattern.replace('#MATERIALS', matinp)
    if self.is_3D:
      labels = self.mesh.nodes.sets['topleft']
      nl = len(labels)
      label = labels[(nl-1)/2]
      pattern = pattern.replace("#3DBOUNDARY", "\nI_SAMPLE.{0}, 3, 3".format(label))
      pattern = pattern.replace('#SURFTYPE', "CYLINDER")
    else:  
      pattern = pattern.replace("#3DBOUNDARY", "")
      pattern = pattern.replace('#SURFTYPE', "SEGMENTS")  
    f =open(self.workdir + self.label + ".inp", 'w')
    f.write(pattern)
    f.close()

  def MakePostProc(self):
    import os, subprocess, time, pickle
    if self.unloading :
        step_number = 2
    else: step_number = 1
    pattern = """# ABQPOSTPROC.PY
# Warning: executable only in abaqus abaqus viewer -noGUI,... not regular python.
import sys
from abapy.postproc import GetFieldOutput_byRpt as gfo
from abapy.postproc import GetVectorFieldOutput_byRpt as gvfo
from abapy.postproc import GetTensorFieldOutput_byRpt as gtfo
from abapy.postproc import GetHistoryOutputByKey as gho
from abapy.indentation import Get_ContactData
from abapy.misc import dump
from odbAccess import openOdb
from abaqusConstants import JOB_STATUS_COMPLETED_SUCCESSFULLY



# Odb opening  
file_name = '#FILE_NAME'
odb = openOdb(file_name + '.odb')
data = {}
step_number= #STEP_NUMBER

# Check job status:
job_status = odb.diagnosticData.jobStatus

if job_status == JOB_STATUS_COMPLETED_SUCCESSFULLY:
  data['completed'] = True
  """
    
    if self.export_fields : pattern += """
  # Field Outputs
  data['field'] = {"U":[], "S":[], "LE":[], "EE":[], "PE":[], "PEEQ":[]}
  fo = data['field']
  for i in xrange(step_number):
      
      fo['U'].append(
        gvfo(odb = odb, 
          instance = 'I_SAMPLE', 
          step = i,
          frame = -1,
          original_position = 'NODAL', 
          new_position = 'NODAL', 
          position = 'node',
          field = 'U', 
          sub_set_type = 'element', 
          delete_report = True))
        
          
      fo['S'].append(
        gtfo(odb = odb, 
          instance = 'I_SAMPLE', 
          step = i,
          frame = -1,
          original_position = 'INTEGRATION_POINT', 
          new_position = 'NODAL', 
          position = 'node',
          field = 'S', 
          sub_set_type = 'element', 
          delete_report = True))
       
      fo['LE'].append(
        gtfo(odb = odb, 
          instance = 'I_SAMPLE', 
          step = i,
          frame = -1,
          original_position = 'INTEGRATION_POINT', 
          new_position = 'NODAL', 
          position = 'node',
          field = 'LE', 
          sub_set_type = 'element', 
          delete_report = True)) 
          
      fo['EE'].append(
        gtfo(odb = odb, 
          instance = 'I_SAMPLE', 
          step = i,
          frame = -1,
          original_position = 'INTEGRATION_POINT', 
          new_position = 'NODAL', 
          position = 'node',
          field = 'EE', 
          sub_set_type = 'element', 
          delete_report = True))  
      
      fo['PE'].append(
        gtfo(odb = odb, 
          instance = 'I_SAMPLE', 
          step = i,
          frame = -1,
          original_position = 'INTEGRATION_POINT', 
          new_position = 'NODAL', 
          position = 'node',
          field = 'PE', 
          sub_set_type = 'element', 
          delete_report = True))
      
      fo['PEEQ'].append(
        gfo(odb = odb, 
          instance = 'I_SAMPLE', 
          step = i,
          frame = -1,
          original_position = 'INTEGRATION_POINT', 
          new_position = 'NODAL', 
          position = 'node',
          field = 'PEEQ', 
          sub_set_type = 'element', 
          delete_report = True))
      """ 
    pattern += """# History Outputs
  data['history'] = {} 
  ho = data['history']
  ref_node = odb.rootAssembly.instances['I_PLATE'].nodeSets['REFNODE'].nodes[0].label
  ho['force'] =  gho(odb,'RF2')['Node I_PLATE.'+str(ref_node)] # GetFieldOutputByKey returns all the occurences of the required output (here 'RF2') and stores it in a dict. Each dict key refers to a location. Here we have to specify the location ('Node I_INDENTER.1') mainly for displacement which has been requested at several locations.
  ho['disp'] =   gho(odb,'U2')['Node I_PLATE.'+str(ref_node)]
  ho['allse'] =   gho(odb,'ALLSE').values()[0]
  ho['allpd'] =   gho(odb,'ALLPD').values()[0]
  ho['allfd'] =   gho(odb,'ALLFD').values()[0]
  ho['allwk'] =   gho(odb,'ALLWK').values()[0]
 
 
else:
  data['completed'] = False
# Closing and dumping
odb.close()
dump(data, file_name+'.pckl')"""  
    pattern = pattern.replace('#FILE_NAME', self.label)  
    pattern = pattern.replace('#STEP_NUMBER', str(step_number))  
    f = open(self.workdir + self.label + '_abqpostproc.py', 'w')
    f.write(pattern)
    f.close()
    
     
    
  
  
  
        
