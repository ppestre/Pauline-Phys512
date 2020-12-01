# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 19:08:37 2020

@author: Pauline
"""
import matplotlib.pyplot as plt
import numpy as np


def lorentz_dev(n):
    '''Random variable with lorentzian distribution'''
    x=np.random.rand(n)
    return np.tan(np.pi*(x/2))


def exp_from_lorentz(x):
    '''Rejection'''
    accept_prob=np.exp(-x)/(1/(1+x**2))
    assert(np.max(accept_prob)<=1)
    accept = np.random.rand(len(x))<accept_prob
    return x[accept]


n=100000

# Sample the bounding function
lor=lorentz_dev(n)

# Do the rejection
rej = exp_from_lorentz(lor)
print('Accepted fraction was %f'  %(len(rej)/len(lor)) )

# Compare to analytic 
transform = -np.log(np.random.rand(len(rej)))

# Count distribution and compare
counts_l,bins=np.histogram(lor[lor<20],200)
counts_r,bins=np.histogram(rej[rej<20],bins=bins)

area=np.sum(counts_r*np.diff(bins))

x=0.5*(bins[1:]+bins[:-1])
# x=bins[:-1]

# Normalize
y_l = counts_l#/area
y_r = counts_r#/area


# Plot it!
fig,[ax1,ax2]=plt.subplots(1,2, figsize=(10,4))

ax1.plot(x,y_l,'.', label='Lorentzian')
ax1.plot(x,y_r,'.', label='Accepted')
# counts,bins,hist = ax1.hist(rej[rej<=20],200,density=True, label='Accepted')
# counts,bins,hist = ax1.hist(lor[lor<=20],200,density=True, label='Lorentzian')

x=0.5*(bins[1:]+bins[:-1])
ax1.plot(x,np.exp(-x)*area,'C3', label='$e^{-t}$')
ax1.set_title('PDF')

ax2.plot(np.sort(transform), 'c',label='transformation')
ax2.plot(np.sort(rej), label='rejection')
ax2.set_title('CDF')

ax1.legend()
ax2.legend()

plt.tight_layout()
plt.savefig('hist.png')