# A repository for compartmentalized models compatible simulation models.


class RingCompression(object):
  """
 Ring compression test.
  """
  def __init__(self, inner_radius, outer_radius, Nr, Nt):
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
