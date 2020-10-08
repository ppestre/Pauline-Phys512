# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 14:01:37 2020

@author: Pauline
"""

import numpy as np
import camb
from matplotlib import pyplot as plt


def A(pars,lmax=2000):
    # Gets the spectrum for a set of parameters
    H0,ombh2,omch2,tau,As,ns=pars[:6]
    pars=camb.CAMBparams()
    pars.set_cosmology(H0=H0,ombh2=ombh2,omch2=omch2,mnu=0.06,omk=0,tau=tau)
    pars.InitPower.set_params(As=As,ns=ns,r=0)
    pars.set_for_lmax(lmax,lens_potential_accuracy=0)
    results=camb.get_results(pars)
    powers=results.get_cmb_power_spectra(pars,CMB_unit='muK')
    cmb=powers['total']
    tt=cmb[:,0]    #you could return the full power spectrum here if you wanted to do say EE
    return tt


def get_grad(x,pars, parsteps):
    grad = np.zeros((len(x),len(parsteps)))
    
    # i = index of output, idx = index of gradient, d = delta of parameter
    for i,step in enumerate(parsteps):

        # Generate vector [0,...,0,ds,0,...0] to add/sub to parameter vector
        delta = np.zeros(len(pars))
        delta[i] = step
        
        # Calculate the spectrum at two points, and calculate derivative
        plus = A(pars + delta)
        mins = A(pars - delta)
        # Set output for this index
        grad[:,i] = (plus - mins)/(2 * step)

    return grad


#def chop(y):
#    if len(y)==2051: return y[2:1201]
#    elif len(y)==1199: return y
#    else: print('Hmmm')

# Load data
wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
xdata,ydata,yerror=wmap[:,0],wmap[:,1],wmap[:,2]





chisq=np.sum((ydata-ymodel[2:1201])**2/yerror**2)


N=np.diag(yerror**2)
Ninv=np.diag(1/yerror**2)

pars_guess=np.array([65,0.02,0.1,0.05,2e-9,0.96])

parsstep=np.array([65,0.02,0.1,0,2e-9,0.96])*0.01

pars=pars_guess.copy()


ymodel=A(pars)
x=np.arange(len(ymodel))


resid=ydata-ymodel[2:1201] #data minus current model
chisq=np.sum(resid**2/yerror**2)

lam=0
for iter in range(10):
    
    dAdm=get_grad(x,pars, 0.01*pars)
    grad=dAdm.T@(Ninv@resid)
    curv=dAdm.T@Ninv@dAdm
    
    chisqnew=chisq
    while chisqnew>=chisq:
        step=np.linalg.inv((1+lam)*curv)@grad
        parsnew=pars+step
        ynew=A(parsnew)
        chisqnew=np.sum((ydata-ynew[2:1201])**2/yerror**2)
        #increase and check again
        if lam==0: lam=1
        lam*=2
        
    #decrease lam
    if lam<1: lam=0
    lam*=np.sqrt(2)
    
    ymodel=ynew
    pars=parsnew
    chisq=chisq
 
        #increase lambda 2 (or +2 if 0)
    
    
    
    print('iteration ',iter,' has step ',step)
    
par_errs=np.sqrt(np.diag(np.linalg.inv(lhs)))
print('final parameters are ',pars_cur,' with errors ',par_errs)


#x=range(len(ymodel))

plt.clf();

#plt.errorbar(wmap[:,0],wmap[:,1],wmap[:,2],fmt='*')
plt.plot(xdata,ydata,'.',label='Data')
plt.plot(x, A(pars_cur), label='Spectrum model')

