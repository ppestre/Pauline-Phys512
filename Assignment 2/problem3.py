# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 12:20:10 2020

@author: Pauline
"""

import numpy as np
from scipy import integrate
import time
from matplotlib import pyplot as plt

half_life=[1.409e17,2.082e6,2.41e4,7.742e12,2.377e12,5.046e10,330350,186,1608,1194,164.3e-6,7.033e8,1.5815e8, 1.1956e7]

def fun(x,y,half_life=half_life):
    dydx=np.zeros(len(half_life)+1)
    dydx[0]=-y[0]/half_life[0]
    for i in range(1,len(half_life)): 
        dydx[i]=y[i-1]/half_life[i-1]-y[i]/half_life[i]
    dydx[-1]=y[-2]/half_life[-2]
    return dydx

#Half lives, start with just uranium
y0=np.zeros(len(half_life)+1)
y0[0]=1

#time interval to look at
ti=0
tf=1e18
#t=np.linspace(ti,tf,100)


t1=time.time();
ans_rk4=integrate.solve_ivp(fun,[ti,tf],y0);
t2=time.time();
print('took ',ans_rk4.nfev,' evaluations and ',t2-t1,' seconds to solve with RK4.')

t1=time.time()
ans_stiff=integrate.solve_ivp(fun,[ti,tf],y0,method='Radau')
t2=time.time()
print('took ',ans_stiff.nfev,' evaluations and ',t2-t1,' seconds to solve implicitly')

# The half life of Uranium is much longer than the others, so the stiff solver is better

ans=ans_stiff

plt.plot(ans.t, ans.y[0])


    
