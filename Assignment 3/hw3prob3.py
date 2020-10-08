# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 14:01:37 2020

@author: Pauline
"""

import numpy as np
import camb
from matplotlib import pyplot as plt


def A(x,pars,lmax=2000):
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
    return tt[x.astype(int)]


def get_grad(x, pars, ds):
    # for ds a list of parameters to vary and their indices
    grad = np.zeros((len(x),len(pars)))
    
    # i = index of output, idx = index of gradient, d = delta of parameter
    for i,pair in enumerate(ds):
        # Get values
        idx = pair[0]
        d = pair[1]
        
        # Generate vector [0,...,0,ds,0,...0] to add/sub to parameter vector
        delta = np.zeros(len(pars))
        delta[idx] = d
        
        # Calculate the spectrum at two points, and calculate derivative
        plus = A(x, pars + delta)
        mins = A(x, pars - delta)
        # Set output for this index
        grad[:,i] = (plus - mins)/(2 * d)
        print('grad row')

    return grad


# Load data
wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
xdata,ydata,yerror=wmap[:,0],wmap[:,1],wmap[:,2]

N=np.diag(yerror**2)
Ninv=np.diag(1/yerror**2)

pars=np.array([65,0.02,0.1,0.05,2e-9,0.96])
ymodel=A(xdata, pars)

resid=ydata-ymodel #data minus current model
chisq=np.sum(resid**2/yerror**2)

# parameters to vary
ds = list(zip([0,1,2,4,5], np.array([65,0.02,0.1,2e-9,0.96])/70))

lam=0
for iter in range(10):
    print('start iter')
    
    dAdm=get_grad(xdata,pars, ds)
    grad=dAdm.T@(Ninv@resid)
    curv=dAdm.T@Ninv@dAdm
    
    chisqnew=chisq
    while chisqnew>=chisq:
        print(chisq)
        print('start step')
        step=np.linalg.pinv((1+lam)*curv)@grad
        parsnew=pars+step
        ynew=A(xdata, parsnew)
        chisqnew=np.sum((ydata-ynew)**2/yerror**2)
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
    
#par_errs=np.sqrt(np.diag(np.linalg.inv(lhs)))
#print('final parameters are ',pars_cur,' with errors ',par_errs)


x=range(2000)

plt.clf();

#plt.errorbar(wmap[:,0],wmap[:,1],wmap[:,2],fmt='*')
plt.plot(xdata,ydata,'.',label='Data')
plt.plot(x, A(x, pars_cur), label='Spectrum model')

