# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 23:12:26 2020

@author: Pauline
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def rat_eval(p,q,x):
    '''
    Evaluates a rational function of x from lists of coefficients p and q.
    Copied from example code provided
    '''
    top=0
    for i in range(len(p)):
        top=top+p[i]*x**i
    bot=1
    for i in range(len(q)):
        bot=bot+q[i]*x**(i+1)
    return top/bot

def rat_fit(x,y,n,m):
    '''
    Fits a rational function of order n/m to points x,y
    Copied from example code provided
    '''
    assert(len(x)==n+m-1)
    assert(len(y)==len(x))
    mat=np.zeros([n+m-1,n+m-1])
    for i in range(n):
        mat[:,i]=x**i
    for i in range(1,m):
        mat[:,i-1+n]=-y*x**i
    # fails for Lorentzian with .inv, succeeds with .pinv    
    pars=np.dot(np.linalg.pinv(mat),y)
    p=pars[:n]
    q=pars[n:]
    return p,q

#%%
 
    
### Cosine
#f  = np.cos
#x0 = np.linspace(-np.pi/2,np.pi/2,11)

## Lorentzian
f  = lambda x: 1/(1+x**2)
x0 = np.linspace(-1,1,11)

y0 = f(x0)

x = np.linspace(x0[0],x0[-1],1001)
y = f(x)

# Order of rational function
n,m=5,7
# Order of polyfit
pl = 6

# Do the fits
y_rat  = rat_eval(*rat_fit(x0,y0,n,m),x)
y_poly = np.polyval(np.polyfit(x0,y0,pl),x)
y_int  = interp1d(x0,y0,kind='cubic')(x)

# Print errors
print('Rational:',   np.std(y_rat -y))
print('Polynomial:', np.std(y_poly-y))
print('Spline:',     np.std(y_int -y))

# Plot points and fits
plt.clf()
plt.plot(x0,y0,'o', label='Data')
plt.plot(x ,y ,':', label='Analytic')
plt.plot(x ,y_rat , label='Rational')
plt.plot(x ,y_poly, label='Polynomial')
plt.plot(x ,y_int , label='Spline')

plt.xlabel('x')
plt.ylabel('y')

plt.legend()

#%%
