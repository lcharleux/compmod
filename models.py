# A repository for compartmentalized models compatible simulation models.
from abapy.mesh import RegularQuadMesh
from abapy.materials import VonMises
import numpy as np

class RingCompression(object):
  """
 Ring compression test.
  """
  def __init__(self, inner_radius = 1., outer_radius = 2., Nr = 8, Nt = 8, disp = .5,  elType = 'CPS4', material =VonMises(labels = 'SAMPLE_MAT') ):
    """
    :param inner_radius: inner radius of the ring
    :type inner_radius: float
    :param outer_radius: outer radius of the ring
    :type outer_radius: float
    :param Nr: Number of elements in the radial direction
    :type Nr: int
    :param Nt: Number of elements in the orthoradial direction
    :type Nt: int
    """
    self.inner_radius = inner_radius
    self.outer_radius = outer_radius
    self.Nr = Nr
    self.Nt = Nt
    self.elType = elType
    self.material = material
    self.disp = disp
    
  def MakeMesh(self):
    """
    """
    Ri = self.inner_radius
    Ro = self.outer_radius
    mesh = RegularQuadMesh(self.Nt, self.Nr, .25, Ro - Ri, name = self.elType)
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
      if nodes.x[i] < 0.: nodes.x[i] = 0. 
    Nr, Nt = self.Nr, self.Nt
    Ne = Nt * Nr
    mesh.add_set('all_elements', mesh.labels)
    mesh.add_set('surface_elements',range( Nt * (Nr-1)+1, Nt*Nr+1  ))
    mesh.add_surface('surface_faces',[ ('surface_elements',3) ])
    
    self.mesh = mesh
  
  def MakeInp(self, path):
    pattern = """**----------------------------------
**RING COMPRESSION SIMULATION
**----------------------------------
**HEADER
*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO, CONTACT=NO
**----------------------------------
** SAMPLE DEFINITION
*PART, NAME = P_SAMPLE
#RING_MESH
*SOLID SECTION, ELSET = ALL_ELEMENTS, MATERIAL = SAMPLE_MAT
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
*SURFACE, TYPE=SEGMENTS, NAME=SURFACE
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
#SAMPLE_MAT
**----------------------------------
** STEPS
**----------------------------------
*STEP, NAME = LOADING0, NLGEOM = YES, INC=1000000
*STATIC, DIRECT
0.001, 1.
*BOUNDARY
I_SAMPLE.LEFT_NODES, 1, 1
I_SAMPLE.RIGHT_NODES, 2, 2
I_PLATE.REFNODE, 2, 2, #DISP
I_PLATE.REFNODE, 1, 1
I_PLATE.REFNODE, 3, 6
*RESTART, WRITE, FREQUENCY = 0
*OUTPUT, FIELD, FREQUENCY = 1
*NODE OUTPUT
COORD, U, 
*ELEMENT OUTPUT, ELSET=I_SAMPLE.ALL_ELEMENTS, DIRECTIONS = YES, POSITION=NODES
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
    pattern = pattern.replace('#RING_MESH', self.mesh.dump2inp())
    pattern = pattern.replace('#SAMPLE_MAT', self.material.dump2inp())
    pattern = pattern.replace('#OUTER_RADIUS', str(self.outer_radius))
    pattern = pattern.replace('#DISP', str(-self.disp))
    f =open(path, 'w')
    f.write(pattern)
    f.close()
  
  
  
  
  
  
        
