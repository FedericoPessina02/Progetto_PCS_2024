from stl import mesh
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

figure = plt.figure()
axes = figure.add_subplot(projection='3d')

your_mesh = mesh.Mesh.from_file(r'../Release/polygonal_mesh.stl')
poly_collection = mplot3d.art3d.Poly3DCollection(your_mesh.vectors,edgecolor='black',linewidths=(1,))
poly_collection.set_color((0.7,0.7,0.7))
axes.add_collection3d(poly_collection)

plt.show()