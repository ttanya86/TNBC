#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 16:28:37 2021

@author: tatianamiti
"""

import os
from matplotlib import pylab
import random
import itertools
import json
import numpy as np
import statistics
from statistics import mean

N = 2
Nsteps = 500
Stroma = []
my_file = open("QuadratStrOn.txt", "r")
StromaStr = my_file.read()
StromaStr = StromaStr.split(", ")
#print(len(StromaStr))
myZeros = 0

#=============================================================================
averNumDiv = Nsteps*[0]
averNumDie = Nsteps*[0]
averTotalGrid = Nsteps*[0]
sdTotalGrid = []
SDaver = []

for ss in range(0,1):
    my_file = open("ContDeath" + str(ss) + ".txt", "r")
    content = my_file.read()
    content = content[1:-1]
    content = content.split("],")
    
    newL = []
    for t in content:
        myS = ''
        for w in t:
            if w != '[':
                if w != ']':
                    myS = myS + w
        newL.append(myS)
        
    ValueList = []
    for ee in newL:
        ValueList.append(ee.split (","))
    
    numDiv = []
    numDie = []
    totalGrid = []
    
    for ff in ValueList:
        numDiv.append(int(ff[0]))
        numDie.append(int(ff[1]))
        totalGrid.append(int(ff[2]))
    if totalGrid[499] > 1094:
        myZeros = myZeros + 1


    
    for gf in range(Nsteps):    
        averNumDiv[gf] = averNumDiv[gf] + numDiv[gf]
        averNumDie[gf] = averNumDie[gf] + numDie[gf]
        averTotalGrid[gf] = averTotalGrid[gf] + totalGrid[gf]
        
         

for ee in range(Nsteps):
    averNumDiv[ee] = averNumDiv[ee]/(N-1) 
    averNumDie[ee] = averNumDie[ee]/(N-1) 
    averTotalGrid[ee] = averTotalGrid[ee]/(N-1) - 1094# len(StromaStr)/2
   
print(mean(averTotalGrid[50:90]) + 1094)
fractTissue = [mean(averTotalGrid[50:90]) + 1094]

print(averTotalGrid[99])
print(averTotalGrid[112])



x = [i/3 for i in range(0,500)]
pylab.figure()
pylab.plot(x, averNumDiv, "r", label = "Dividing Cells")
pylab.plot(x, averNumDie, "b", label = "Dying Cells")
pylab.xlabel("Time (days)")
pylab.ylabel("# Tumor Cells")
pylab.legend()
pylab.title("Chemotherapy & 0% Proliferation Excess & 0 Lingering Death", fontsize=8)
#pylab.text(100, mean(averNumDiv[20:70])+10, "ALectinib",fontsize= 8 )
#pylab.text(109, mean(averNumDiv[50:90])+8, "ALectinib + Reduced Stroma",fontsize=8  )
pylab.savefig("DivDieAver0perc3rad.png", dpi = 300)


NumDivFract = []
NumDieFract = []
for gg in range(0,500):
    if averTotalGrid[gg] != 0:
    	NumDivFract.append(averNumDiv[gg]/averTotalGrid[gg])
    	NumDieFract.append(averNumDie[gg]/averTotalGrid[gg])
    else:
    	NumDivFract.append(0)
    	NumDieFract.append(0)

print(mean(NumDivFract[50:90]))
print(mean(NumDieFract[50:90]))
FractZeros = [myZeros]


pylab.figure()
pylab.plot(x[10:], NumDivFract[10:], "r", label = "Dividing Cells")
pylab.plot(x[10:], NumDieFract[10:], "b", label = "Dying Cells")
pylab.xlabel("Time(days")
pylab.ylabel("# Tumor Cells")
pylab.legend()
pylab.ylim(0, 1)
pylab.title("Chemotherapy & 0% Proliferation Excess & 0 Lingering Death", fontsize=8)
#pylab.text(100, mean(averNumDiv[50:70])+10, "Untreated",fontsize= 8 )
#pylab.text(109, mean(averNumDie[50:90])+8, "Chemotherapy",fontsize=8  )
pylab.savefig("DivDieAverfrct0perc3rad.png", dpi = 300)



pylab.figure()
pylab.plot(x[10:], averTotalGrid[10:], "b", label = "Total Cells")
#pylab.errorbar(x[10:], averTotalGrid[10:], yerr = SDaver[10:], label = "SD")
pylab.xlabel("Time (days)")
pylab.ylabel("# Tumor Cells")
pylab.legend()
pylab.ylim(0, 10000)
pylab.title("Chemotherapy & 0% Proliferation Excess & 0 Lingering Death", fontsize=8)
#pylab.text(200, mean(averTotalGrid[20:90])+100, "Untreated",fontsize= 8 )
#pylab.text(79, mean(averTotalGrid[200:290])+50, "Chemotherapy",fontsize=8  )
pylab.savefig("TotalAver0perc3rad.png", dpi = 300)


if os.path.isfile('averTotalGrid.txt'):
     with open('averTotalGrid.txt', 'w') as f_ffma:
         f_ffma.write(json.dumps(averTotalGrid))
else:
     with open('averTotalGrid.txt', 'w') as f_ffma:
         f_ffma.write(json.dumps(averTotalGrid))


if os.path.isfile('FractTissue.txt'):
     with open('FractTissue.txt', 'w') as f_ff:
         f_ff.write(json.dumps(fractTissue))
else:
     with open('FractTissue.txt', 'w') as f_ff:
         f_ff.write(json.dumps(fractTissue))

if os.path.isfile('FractZeros.txt'):
     with open('FractZeros.txt', 'w') as f_ffd:
         f_ffd.write(json.dumps(FractZeros))
else:
     with open('FractZeros.txt', 'w') as f_ffd:
         f_ffd.write(json.dumps(FractZeros))

from statistics import mean
#FractSD = [statistics.stdev(FractZeros)]
     
#if os.path.isfile('FractSD.txt'):
 #    with open('FractSD.txt', 'w') as f_ffmaa:
  #       f_ffmaa.write(json.dumps(FractSD))
#else:
 #    with open('FractSD.txt', 'w') as f_ffmaa:
  #       f_ffmaa.write(json.dumps(FractSD))
         
## =============================================================================
# if os.path.isfile('averNumDie.txt'):
#     with open('averNumDie.txt', 'w') as f_ffma:
#         f_ffma.write(json.dumps(averNumDie))
# else:
#     with open('averNumDie.txt', 'w') as f_ffma:
#         f_ffma.write(json.dumps(averNumDie))
#         
#         
# if os.path.isfile('averTotalGrid.txt'):
#     with open('averTotalGrid.txt', 'w') as f_ffma:
#         f_ffma.write(json.dumps(averTotalGrid))
# else:
#     with open('averTotalGrid.txt', 'w') as f_ffma:
#         f_ffma.write(json.dumps(averTotalGrid))
#   
#     
# if os.path.isfile('AverMean.txt'):
#     with open('AverMean.txt', 'w') as f_ffma:
#         f_ffma.write(json.dumps(AverMean))
# else:
#     with open('AverMean.txt', 'w') as f_ffma:
#         f_ffma.write(json.dumps(AverMean))
# =============================================================================
# =============================================================================
