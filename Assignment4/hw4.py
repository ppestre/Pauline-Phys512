# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 11:58:47 2020

@author: Pauline
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft,irfft,fft,ifft


def gaussian(x, mu, sig):
    return np.exp(-(x - mu)**2 / (2 * sig**2))

def roll(array,n):
    N=len(array)
    k=np.arange(N)
    return np.real(ifft(fft(array)*np.exp(-2j*np.pi*n*k/N)))

def corr(f,g):
    return irfft(rfft(f)*np.conjugate(rfft(g)))

def corr_pad(f,g):
    N=min(len(f),len(g))
    f_pad, g_pad=np.zeros(2*len(f)),np.zeros(2*len(g))
    f_pad[:len(f)],g_pad[:len(g)]=f,g
    return irfft(rfft(f_pad)*np.conjugate(rfft(g_pad)))[:N]

#%%
'''Question 1'''

x=np.linspace(0,100,100)
y=gaussian(x,50,10)
yshift=roll(y,50)

#plt.plot(x,y, label="Initial Gaussian")
#plt.plot(x,yshift, label="Shifted Gaussian")
##plt.plot(x,y*yshift, label="Multiplied")
#plt.legend()

#plt.savefig('shifted')

#%%
'''Question 2'''

#plt.plot(corr(y,y))
#plt.xlabel('Shift')
#plt.ylabel('$f*g$')
#plt.legend()
#
#plt.savefig('corr')

#%%
'''Question 3'''

for n in [0,20,50,90]:plt.plot(corr(y,roll(y,n)), label="n = %i" % n)
plt.xlabel('Relative Shift')
plt.ylabel('$f*g$')
plt.legend()

plt.savefig('corrshift')

#%%
'''Question 4'''


plt.plot(corr_pad(y,y))
plt.xlabel('Shift')
plt.ylabel('$f*g$')


#plt.savefig('corr')

#%%
'''Question 5'''

T=1. #acquisition time
F=100. #sampling frequency
N=T*F

t=np.arange(0,T,1/F)
f=np.arange(-F/2, F/2, 1/T)

mycos = lambda cycles: 2*np.cos(2*np.pi*(cycles*F)*t/N)
window = 0.5-0.5*np.cos(2*np.pi*t*F/N)

#plt.plot(t,window)
plt.plot(f,np.abs(fft(mycos(4.5))), label="4.5 cycles")
plt.plot(f,np.abs(fft(mycos(9))), label="9 cycles")
plt.plot(f,np.abs(fft(mycos(4.5)*window)), label="4.5 cycles windowed")
plt.legend()

plt.savefig('dft')

plt.plot(np.abs(fft(window)))
plt.savefig('window')
