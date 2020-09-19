# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 17:43:14 2020

@author: Pauline
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

data = np.genfromtxt('lakeshore.txt')

T0,V0,_ = np.transpose(data)

V=np.linspace(V0[0],V0[-1],1001)

TofV=interp1d(V0,T0,kind='cubic')

maxT = interp1d(V0, T0, kind='previous')
minT = interp1d(V0, T0, kind='next')

plt.clf()

plt.plot(V,TofV(V),'-', label='Interpolation', zorder=10,linewidth=2.5)
plt.plot(V,minT(V),'-', label='Min')
plt.plot(V,maxT(V),'-', label='Max')
plt.plot(V0,T0,'.', label='Data', zorder=10)

plt.xlabel('Voltage (V)')
plt.ylabel('Temperature (K)')
plt.legend()
#plt.savefig('TofV.png')

#%%

#Try a voltage value here!
yourvalue=1.2

if yourvalue>np.max(V0) or yourvalue<np.min(V0): print("Hmm that's not right")
else: print("At your value %s V, the temperature is %.1f K"%(yourvalue, TofV(yourvalue)))

