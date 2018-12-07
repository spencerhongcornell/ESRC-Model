from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pickle
import numpy as np
realCOPl = pickle.load(open("COPl","rb"))

x = np.arange(0, 1, 0.01)
y = np.arange(370, 385,1)
z = realCOPl


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(x, y, z)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
