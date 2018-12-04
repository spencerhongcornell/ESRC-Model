import numpy as np
import matplotlib.pyplot as plt
import pandas

colnames = ['pressure', 'ammoniacomp', 'temperature', 'vaporwater', 'vaporammonia', 'liquidwater', 'liquidammonia']
filename = 'txy-ammonia'
data = pandas.read_csv('%s.csv' %filename, names = colnames)

ammoniacomp = data.ammoniacomp.tolist()
temperature = data.temperature.tolist()
vaporammonia = data.vaporammonia.tolist()
liquidammonia = data.liquidammonia.tolist()

plt.plot(vaporammonia, temperature, color = 'black')
plt.plot(liquidammonia, temperature, color = "#1775f0")
plt.savefig("test.svg", format="svg")