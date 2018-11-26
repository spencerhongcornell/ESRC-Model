import numpy as np
import matplotlib.pyplot as plt
from sympy import solve, Eq, symbols
import sys 
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

def entropyammoniawater(targetcomp):
	# done at 4 bar

	entropy = np.array([1.52,
	1.22,
	0.91,
	0.59,
	0.28,
	0.02,
	-0.17,
	-0.28,
	-0.29,
	-0.3,
	-0.31])

	lowerentropy = entropy[int(math.floor(targetcomp *10))]
	higherentropy = entropy[int((math.floor(targetcomp * 10)))+1]
	return lowerentropy, higherentropy
print(entropyammoniawater(0.78))

#print(absorberevaporator(10, 0.45, 0.8, 0.7))


