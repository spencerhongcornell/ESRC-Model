import numpy as np
import matplotlib.pyplot as plt
from sympy import solve, Eq, symbols
import sys 
import pandas
import math
# This function calculates the degradation of the COP based on each entropy generation term (S)
def degradeCOP(Tevap, Tcond, Qall, S):
	degraded = ((Tevap * Tcond)/(Tcond - Tevap)) * (S/Qall)
	return degraded

# This function calculates the reversible COP of a ES refrigerator based on thermal reservoirs
def reversibleCOP(Tevap, Tcond, Tgen):
	revCOP = (((Tgen - Tcond)/(Tgen))/((Tcond-Tevap)/Tevap))
	return revCOP

def absorberevaporator(m6, xa4, ya3, xa6):
	m4, m3, m5= symbols(['m4', 'm3', 'm5'])
	system = [
	Eq((xa4*m4)+ (ya3 * m3) - (ya3 * m5) - (xa6 * m6), 0),
	Eq(m4 + m3 - m5 - m6, 0),
	Eq(m4 - (ya3*m3) - (m5) + (ya3 * m6), 0)
	]
	soln = solve(system, [m4, m3, m5, m6])
	return soln

	#a = np.array([[xa4, ya3, -ya3, -xa6], [1, 1, -1, -1], [1, -ya3, -1, ya3]])
	#b = np.array([0, 0, 0])

	#return np.linalg.solve(a, b)

def generator(m1):
	m6 = m1
	return m6

def interpolate(filename, targetcomp):
	# done at 4 bar
	# must use the entropy-ammonia csv
	# WORKING!
	colnames = ['pressure', 'ammoniacomp', 'entropy']
	data = pandas.read_csv('%s.csv' %filename, names=colnames)

	ammoniacomp = data.ammoniacomp.tolist()
	entropy = data.entropy.tolist()

	#print(ammoniacomp)
	#print(entropy)

	#print(int(math.floor(targetcomp /0.05))*0.05 )
	#print((int((math.floor(targetcomp /0.05))+1))*0.05)
	lowerentropy = entropy[int(math.floor(targetcomp /0.05))]
	higherentropy = entropy[(int((math.floor(targetcomp /0.05))+1))]

	theta = (targetcomp - int(math.floor(targetcomp /0.05))*0.05 ) / ((int((math.floor(targetcomp /0.05))+1))*0.05 - int(math.floor(targetcomp /0.05))*0.05 )
	return (theta * higherentropy) + (1-theta)*lowerentropy

def leverrule(inputflow, temp, inputcomp):
	#t-xy of ammonia-water
	colnames = ['pressure', 'ammoniacomp', 'temperature', 'vaporwater', 'vaporammonia', 'liquidwater', 'liquidammonia']
	filename = 'txy-ammonia'
	data = pandas.read_csv('%s.csv' %filename, names = colnames)

	ammoniacomp = data.ammoniacomp.tolist()
	temperature = data.temperature.tolist()
	vaporammonia = data.vaporammonia.tolist()
	liquidammonia = data.liquidammonia.tolist()


	index, valuetemp =  min(enumerate(temperature), key=lambda x: abs(x[1]-temp))

	liquiddistance = inputcomp - liquidammonia[index]
	vapordistance = vaporammonia[index] - inputcomp

	vaporflow = symbols('vaporflow')
	system = [
	Eq((vapordistance * vaporflow) + (-1.0*liquiddistance*(float(inputflow) - vaporflow)), 0)
	]
	soln = solve(system, [vaporflow])
	return soln, (inputflow - soln[vaporflow]) ,liquidammonia[index], vaporammonia[index]

	# load T-xy diagram
	# do leverrule to calculate the two "arms"
print(leverrule(50, 345, 0.5))

