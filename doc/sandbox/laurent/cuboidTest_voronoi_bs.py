# Voronoi cells for grain boundary generation ?
# A demo in 2D (3D is exaclty the same)

import numpy as np
import abapy
from scipy import interpolate

  
# SIMULATION SETTINGS
lx, ly, lz = 1., 1., 1.
Nx, Ny, Nz = 5, 5, 5
Ne = Nx * Ny * Nz


Nseed = 200



 
mesh = abapy.mesh.RegularQuadMesh(l1 = lx,
                                  l2 = ly,
                                  N1 = Nx,
                                  N2 = Ny)
mesh = mesh.extrude(l= lz, N = Nz)                                  
                                  
centroids = mesh.centroids()
nodes = mesh.nodes
nodes_postion = np.array([nodes.x, nodes.y, nodes.z]).transpose()
conn = np.array(mesh.connectivity)
labels = nodes.labels

# VORONOI CELL BUILDING

# Bounding box
xmin, xmax = nodes_postion[:,0].min(), nodes_postion[:,0].max()
ymin, ymax = nodes_postion[:,1].min(), nodes_postion[:,1].max() 
zmin, zmax = nodes_postion[:,2].min(), nodes_postion[:,2].max() 

# Seeds
seeds = np.random.rand(Nseed,3) * np.array([[xmax - xmin, ymax -ymin, zmax - zmin]]) + np.array([[xmin, ymin, zmin]])
seed_flags = np.arange(Nseed)
elem_flags = interpolate.griddata(seeds, seed_flags, centroids, method = "nearest")

# Building sets
elabels = np.array(mesh.labels)


for i in xrange(len(seed_flags)):
  flag = seed_flags[i]
  eset = elabels[np.where(elem_flags == flag)[0]]
  if len(eset) != 0 : 
    mesh.add_set(label = "RAW_SET{0:06d}".format(i), elements = eset)
  
inp = mesh.dump2inp().split("RAW_SET")
inp2 = inp[0]
for i in xrange(1, len(inp)):
  inp2 += "G{0}".format(i)
  inp2 += inp[i][6:]
  
print inp2  
  



