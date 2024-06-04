from stl import mesh
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

figure = plt.figure()
axes = figure.add_subplot(projection='3d')

data = mesh.Mesh.from_file(r'../Release/polygonal_mesh.stl')
poly_collection = mplot3d.art3d.Poly3DCollection(data.vectors)
poly_collection.set_edgecolor("red")
poly_collection.set_facecolor("grey")
axes.add_collection3d(poly_collection)

scale = data.points.flatten()
axes.auto_scale_xyz(scale, scale, scale)
plt.show()