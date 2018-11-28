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

def massabsorberevaporator(m6, xa4, ya3, xa6):
	m4, m3, m5= symbols(['m4', 'm3', 'm5'])
	system = [
	Eq((xa4*m4)+ (ya3 * m3) - (ya3 * m5) - (xa6 * m6), 0),
	Eq(m4 + m3 - m5 - m6, 0),
	Eq(m4 - (ya3*m3) - (m5) + (ya3 * m6), 0)
	]
	soln = solve(system, [m4, m3, m5, m6])
	return soln

def interpolate(filename, targetcomp):
	# done at 4 bar
	# must use the entropy-ammonia-water csv, entropy-ammonia-butane csv, or enthalpy-ammonia-water csv
	# WORKING!
	colnames = ['pressure', 'ammoniacomp', 'prop']
	data = pandas.read_csv('%s.csv' %filename, names=colnames)

	ammoniacomp = data.ammoniacomp.tolist()
	prop = data.prop.tolist()

	lower = prop[int(math.floor(targetcomp /0.05))]
	higher = prop[(int((math.floor(targetcomp /0.05))+1))]

	theta = (targetcomp - int(math.floor(targetcomp /0.05))*0.05 ) / ((int((math.floor(targetcomp /0.05))+1))*0.05 - int(math.floor(targetcomp /0.05))*0.05 )
	return (theta * higher) + (1-theta)*lower

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

def Qgenerator(massin, compin):
	massout = massin

	enthalpyin = interpolate('enthalpy-kjmol-266K-ammoniawater', compin)

	enthalpyout = interpolate('enthalpy-kjmol-375K-ammoniawater', compin)

	Qgen = -1*(enthalpyin - enthalpyout)

	return Qgen

def Sgenerator(massin, compin, Qgen):
	massout = massin

	entropyin = interpolate('entropy-kjmol-266K-ammoniawater', compin)

	entropyout = interpolate('entropy-kjmol-375K-ammoniawater', compin)

	Sgen = symbols('Sgen')
	system = [
	Eq((-1 * massin * entropyin ) + (massout*entropyout) + (Qgen/375) - Sgen, 0)

	]
	soln = solve(system, [Sgen])

	return soln[Sgen]

Qgen= Qgenerator(float(sys.argv[1]), float(sys.argv[2]))
Sgen = Sgenerator(float(sys.argv[1]), float(sys.argv[2]), Qgen)

print(Qgen, Sgen)


