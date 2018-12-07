import pickle

from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt

realCOPl = pickle.load(open("COPl","rb"))
x = np.arange(0, 1, 0.01)
y = np.arange(370, 385,1)
z = np.array(realCOPl)

xi= np.linspace(0, 1, len(z))
yi = np.linspace(370, 385, len(z))

X,Y= np.meshgrid(xi,yi)
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
plt.contourf(xi, yi, zi)

plt.show() 