from compmod.models import RingCompression
import matplotlib.pyplot as plt


model = RingCompression()
model.MakeMesh()
mesh = model.mesh
x,y,z = mesh.get_edges() # Mesh edges
X,Y,Z,tri = mesh.dump2triplot()
xb,yb,zb = mesh.get_border() # mesh borders
xe, ye, ze = mesh.get_edges()
fig = plt.figure(figsize=(10,10))
fig.gca().set_aspect('equal')
#plt.axis('off')
plt.grid()
plt.plot(xb,yb,'k-', linewidth = 2.)
plt.plot(xe, ye,'k-', linewidth = .5)
plt.triplot(X,Y,tri)
plt.show()

