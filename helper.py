import numpy as np
import matplotlib.pyplot as plt
from sympy import solve, Eq, symbols
import sys 
import pandas
import math
import os

def degradeCOP(Tevap, Tcond, Qall, S):
    degraded = ((Tevap * Tcond)/(Tcond - Tevap)) * (S/Qall)
    return degraded

# This function calculates the reversible COP of a ES refrigerator based on thermal reservoirs
def reversibleCOP(Tevap, Tcond, Tgen):
    revCOP = (((Tgen - Tcond)/(Tgen))/((Tcond-Tevap)/Tevap))
    return revCOP

# This function solves the system of equations to calculate the mass flowrates of
# the combined absorber-evaporator system
def massabsorberevaporator(m6, m4, xa4, ya3, xa6):
    m3, m5= symbols(['m3', 'm5'])
    system = [
    #Eq((xa4*m4)+ (ya3 * m3) - (0 * m5) - (xa6 * m6), 0),
    Eq(m5 - (1-ya3)*m3,0),
    Eq(m4 + m3 - m5 - m6, 0),
    #Eq((1-ya3)*m3-m5, 0)
    #Eq(m4 - (ya3*m3) - (m5) + (ya3 * m6), 0)
    ]
    soln = solve(system, [m4, m3, m5, m6])

    return float(m4), float(soln[m3]), float(soln[m5]), float(m6)

# This is an interpolate helper function to be used in other functions.
# targetcomp refers to ammonia composition. All CSV files are in ammonia composition.
def interpolateAW(targetcomp, T, var):
    # must use the entropy-ammonia-water csv, entropy-ammonia-butane csv, or enthalpy-ammonia-water csv

    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    rel_path = 'data/mixed/'
    abs_file_path = os.path.join(script_dir, rel_path)

    colnames = ['pressure', 'ammoniacomp', 'prop']
    propertyname = ''
    if var is 1:
        propertyname = 'enthalpy'
    else:
        propertyname = 'entropy'

    filename = propertyname + '-kjmol-' + str(T) + 'K-ammoniawater'

    data = pandas.read_csv(str(abs_file_path) + '%s.csv' %filename, names=colnames)

    ammoniacomp = data.ammoniacomp.tolist()
    prop = data.prop.tolist()

    lower = prop[int(math.floor(targetcomp /0.05))]
    higher = prop[(int((math.floor(targetcomp /0.05))+1))]

    theta = (targetcomp - int(math.floor(targetcomp /0.05))*0.05 ) / ((int((math.floor(targetcomp /0.05))+1))*0.05 - int(math.floor(targetcomp /0.05))*0.05 )
    return (theta * higher) + (1-theta)*lower

def interpolateAB(targetcomp, T, var):
    # must use the entropy-ammonia-water csv, entropy-ammonia-butane csv, or enthalpy-ammonia-water csv

    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    rel_path = 'data/mixed/'
    abs_file_path = os.path.join(script_dir, rel_path)

    colnames = ['pressure', 'ammoniacomp', 'prop']
    propertyname = ''
    if var is 1:
        propertyname = 'enthalpy'
    else:
        propertyname = 'entropy'

    filename = propertyname + '-kjmol-' + str(T) + 'K-ammoniabutane'

    data = pandas.read_csv(str(abs_file_path) + '%s.csv' %filename, names=colnames)

    ammoniacomp = data.ammoniacomp.tolist()
    prop = data.prop.tolist()

    lower = prop[int(math.floor(targetcomp /0.05))]
    higher = prop[(int((math.floor(targetcomp /0.05))+1))]

    theta = (targetcomp - int(math.floor(targetcomp /0.05))*0.05 ) / ((int((math.floor(targetcomp /0.05))+1))*0.05 - int(math.floor(targetcomp /0.05))*0.05 )
    return (theta * higher) + (1-theta)*lower
# This calculates the two mass flowrates and the compositions coming out of the flash drum 
# given a mass flowrate and composition of ammonia coming in
# inputcomp is the ammonia composition
# temp is the temperature that the flash drum flashes at
def leverrule(inputflow, temp, inputcomp):
    #t-xy of ammonia-water
    #input composition of ammonia
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    rel_path = "data/txy/"
    abs_file_path = os.path.join(script_dir, rel_path)
    colnames = ['pressure', 'ammoniacomp', 'temperature', 'vaporwater', 'vaporammonia', 'liquidwater', 'liquidammonia']
    filename = 'txy-ammonia-4bar'
    data = pandas.read_csv( str(abs_file_path) + '%s.csv' %filename, names = colnames)

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

    # the order is: vapor flow, liquid flow, liquid ammonia composition. vapor ammonia composition
    return soln[vaporflow], (inputflow - soln[vaporflow]) ,liquidammonia[index], vaporammonia[index]

# This calculates the Q of the generator
# compin is the ammonia composition
def Qgenerator(massin, compin, Tgen):
    massout = massin

    enthalpyin = interpolateAW(compin, 325, 1 )

    enthalpyout = interpolateAW(compin, Tgen, 1)

    Qgen = -1*(massin*enthalpyin - massout*enthalpyout)

    return Qgen
# This calculates the S of the flash
# compin is the ammonia flash
def Sgenerator(massin, compin, Qgen, Tgen):
    massout = massin

    entropyin = interpolateAW(compin, 325, 2)
    #RAHUL fixed Line 95 - wrong entropy values read in
    entropyout = interpolateAW(compin, Tgen, 2)

    Sgen = symbols('Sgen')
    system = Eq((-1 * massin * entropyin ) + (massout*entropyout) + (Qgen/Tgen) - Sgen, 0)

    soln = solve([system], Sgen)

    return soln[Sgen]

def Qflash(massin, massvapor, massliquid, compin, vaporammonia, liquidammonia, Tgen):

    enthalpyin = interpolateAW( compin, Tgen, 1)

    enthalpyoutvapor = interpolateAW(vaporammonia, Tgen, 1)

    enthalpyoutliquid = interpolateAW( liquidammonia, Tgen, 1)

    Qflash = symbols('Qflash')
    system = Eq(( massin* enthalpyin ) + (-1* massvapor*enthalpyoutvapor) + (-1* massliquid*enthalpyoutliquid) + Qflash, 0)

    
    soln = solve([system], [Qflash])

    return soln[Qflash]
# This calculates the S of the flash
# compin is the ammonia flash
def Sflash(massin, massvapor, massliquid, compin, vaporammonia, liquidammonia, Qflash, Tgen):

    entropyin = interpolateAW( compin, Tgen, 2)

    entropyoutvapor = interpolateAW(vaporammonia, Tgen, 2)

    entropyoutliquid = interpolateAW(liquidammonia, Tgen, 2)

    Sflash = symbols('Sflash')
    system = Eq(( massin* entropyin ) + (-1* massvapor*entropyoutvapor) + (-1* massliquid*entropyoutliquid) + (Sflash/Tgen) - Qflash, 0)

    
    soln = solve([system], Sflash)

    return soln[Sflash]
def Qevaporator(m2, m3, m5, ya3, ya2, xa5, Tgen):

    enthalpym2 = interpolateAW(ya2, Tgen, 1)

    enthalpym3 = interpolateAB( ya3, 266, 1)
    #print(enthalpym3)
    enthalpym5 = interpolateAB( xa5, 325, 1)
    #print(enthalpym5)


   # print(enthalpym2 + enthalpym3 + enthalpym5)
    Qevap = symbols('Qevap')
    system = Eq(( m2 * enthalpym2 ) + (-1* m3*enthalpym3) + (m5*enthalpym5) + Qevap, 0)
    soln = solve([system], Qevap)

    #print(type(soln))
    return soln[Qevap]

def Sevaporator(m2, m3, m5, ya3, ya2, xa5, Qevap, Tgen):

    entropym2 = interpolateAW(ya2, Tgen, 2)

    entropym3 = interpolateAB( ya3, 266, 2)

    entropym5 = interpolateAB( xa5, 325, 2)

    Sevap = symbols('Sevap')
    system = Eq(( m2 * entropym2 ) + (-1* m3*entropym3) + (m5*entropym5) + (Qevap/266) - Sevap, 0)
    

    soln = solve([system], Sevap)

    return soln[Sevap]

def Qabsorber(m3, m4, m5, m6, ya3, xa4, xa5, xa6, Tgen):

    enthalpym3 = interpolateAB( ya3, 266, 1)
    
    enthalpym4 = interpolateAW( xa4, Tgen, 1)
    
    enthalpym5 = interpolateAB(xa5, 325, 1)
    
    enthalpym6 = interpolateAW( xa6, 325, 1)
    
    Qabs = symbols('Qabs')
    system =  (m3 * enthalpym3 ) + (m4 * enthalpym4) + (-1*m5 * enthalpym5) + (-1 * m6 * enthalpym6) + Qabs
    

    soln = solve([system], Qabs)

    return soln[Qabs]


def Sabsorber(m3, m4, m5, m6, ya3, xa4, xa5, xa6, Qabs, Tgen):

    entropym3 = interpolateAB( ya3, 266, 2)

    entropym4 = interpolateAW( xa4, Tgen, 2)

    entropym5 = interpolateAB( xa5, 325, 2)

    entropym6 = interpolateAW(xa6, 325, 2)

    Sabs = symbols('Sabs')
    system =  Eq((m3*entropym3) + (m4 * entropym4) + (-1*m5 * entropym5) + (-1*m6*entropym6) + (Qabs/325)- Sabs, 0)
    

    soln = solve([system], Sabs)

    return soln[Sabs]