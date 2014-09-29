# A repository for compartmentalized models compatible simulation models.
from abapy.mesh import RegularQuadMesh
import numpy as np

class RingCompression(object):
  """
 Ring compression test.
  """
  def __init__(self, inner_radius = 1., outer_radius = 2., Nr = 8, Nt = 8, elType = 'CPS4'):
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
    
  def MakeMesh(self):
    """
    """
    Ri = self.inner_radius
    Ro = self.outer_radius
    mesh = RegularQuadMesh(self.Nt, self.Nr, .251, Ro - Ri, name = self.elType)
    def function(x, y, z, labels):
      theta = 2 * np.pi * x
      r = y + Ri
      ux = -x + r * np.cos(theta)
      uy = -y + r * np.sin(theta)
      uz = 0. * z
      return ux, uy, uz
    vectorField = mesh.nodes.eval_vectorFunction(function)
    mesh.nodes.apply_displacement(vectorField)
    self.mesh = mesh
      
