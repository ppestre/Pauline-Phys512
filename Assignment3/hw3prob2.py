# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 12:49:41 2020

@author: Pauline
"""

import numpy as np
import camb
from matplotlib import pyplot as plt


def get_spectrum(pars,lmax=2000):
    #print('pars are ',pars)
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


wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')

xdata,ydata,yerror=wmap[:,0],wmap[:,1],wmap[:,2]

pars=np.array([65,0.02,0.1,0.05,2e-9,0.96])

ymodel=get_spectrum(pars)

plt.clf();

#plt.errorbar(wmap[:,0],wmap[:,1],wmap[:,2],fmt='*')
plt.plot(xdata,ydata,'.',label='Data')
plt.plot(ymodel, label='Spectrum model')


chisq=np.sum((ydata-ymodel[2:1201])**2/yerror**2)
