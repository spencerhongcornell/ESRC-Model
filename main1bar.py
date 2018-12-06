from helper1bar import *
import numpy as np 
from sympy import solve, Eq, symbols
import csv 
import pandas
import math
import os
import pickle

def run(ammoniacompin, Tgen):
	massintoflash = 10
	massvaporoutflash, massliquidoutflash, liquidammoniacomp, vaporammoniacomp = leverrule(massintoflash, Tgen, ammoniacompin)

	m2 = massvaporoutflash
	ya2 = vaporammoniacomp

	m4 = massliquidoutflash

	m6old = massintoflash
	xa4 =  liquidammoniacomp
	ya3 = 0.75    #AZEOTROPE
	xa6 = ammoniacompin
	xa5 = 0        #ASSUME TO BE PURE STREAM - RAHUL


	m4, m3, m5, m6new = massabsorberevaporator(m6old, massliquidoutflash, xa4, ya3, xa6)
	#print("m2: " +  str(m2) + "\nm3: " + str(m3) +  "\n" + "m4: " + str(m4)+ "\n" + "m5: " + str(m5) + "\n" + "m6: " + str(m6new))

	Qevap = Qevaporator(m2, m3, m5, ya3, ya2, xa5, Tgen)
	Sevap = Sevaporator(m2, m3, m5, ya3, ya2, xa5, Qevap, Tgen)

	Qflas = Qflash(massintoflash, m2, m4, ammoniacompin, ya2, xa4, Tgen)
	Sflas = Sflash(massintoflash, m2, m4, ammoniacompin, ya2, xa4, Qflas, Tgen)


	#print("Qevap: " + str(Qevap) + "\nSevap: " + str(Sevap))


	Qabs = Qabsorber(m3, m4, m5, m6new, ya3, xa4, xa5, xa6, Tgen)
	Sabs = Sabsorber(m3, m4, m5, m6new, ya3, xa4, xa5, xa6, Qabs, Tgen)

	#print("Qabs: " + str(Qabs) + "\nSabs: " + str(Sabs))

	Qgen = Qgenerator(m6new, xa6,Tgen)
	Sgen = Sgenerator(m6new, xa6,Qgen, Tgen)

	#print("Qgen: " + str(Qgen) + "\nSgen: " + str(Sgen))

	COPreal = (Qevap/Qgen)
	#Sall = Sgen + Sevap + Sabs
	COPrev = reversibleCOP(266, 325, Tgen)

	#print(COPrev)
	#Sall = - Qevap/266 - Qabs/325 + Qgen/375

	# Sall = [Sgen, Sevap, Sabs]
	# COPdegrade = 0
	# for S in Sall:
	#     COPdegrade += degradeCOP(266, 325, Qgen, S)

	# print(COPdegrade)
	#COPdegrade = degradeCOP(266, 325, 4000000, Sall)
	return Qevap, m4, m3, Qgen, COPreal, Qflas, Qabs, Sabs, Sflas


m4l = []
Qevapl = []
m3l = []
Qgenl = []

Qflasl = []

with open("output.csv", "w") as f:
    writer = csv.writer(f)
    overall = [np.arange(369, 385, 1)]
    for i in np.arange(0, 1, 0.01):
        print('read at %.2f'%i)
        realCOPl = [i]
        for T in np.arange(370, 375, 1):
            Qevap, m4, m3, Qgen, COPreal, Qflas, Qabs, Sabs, Sflas= run(i, T)
            print('______')
            if Qevap > 0 and m4 > 0 and m3 > 0 and Qgen > 0 and Qabs < 0:
                Qevapl.append(Qevap)
                m4l.append(m4)
                m3l.append(m3)
                Qgenl.append(Qgen)
                Qflasl.append(Qflas)
                realCOPl.append(COPreal)
                print('success')
                print('________')
        overall.append(realCOPl)



    writer.writerows(overall)





#pickle.dump(realCOPl,open("COPl","wb"))

 