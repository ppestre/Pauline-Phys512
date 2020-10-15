# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 18:24:26 2020

@author: Po
"""

import numpy as np
import camb
from matplotlib import pyplot as plt


def A(x,pars,lmax=2000):
    print(pars)
    # Gets the spectrum for a set of parameters
    H0,ombh2,omch2,As,ns,tau=pars
    pars=camb.CAMBparams()
    pars.set_cosmology(H0=H0,ombh2=ombh2,omch2=omch2,mnu=0.06,omk=0,tau=tau)
    pars.InitPower.set_params(As=As,ns=ns,r=0)
    pars.set_for_lmax(lmax,lens_potential_accuracy=0)
    results=camb.get_results(pars)
    powers=results.get_cmb_power_spectra(pars,CMB_unit='muK')
    cmb=powers['total']
    tt=cmb[:,0]    
    print('done')
    return tt[x.astype(int)]



def get_chisq(data,pars):
    #we need a function that calculates chi^2 for us for the MCMC
    #routine to call
    x=data[0]
    y=data[1]
    yerror=data[2]
    model=A(x,pars)
    chisq=np.sum( (y-model)**2/yerror**2)
    return chisq

def check_pars(pars):
    if pars[-1]<0: print('nooo')
    return


def run_mcmc(pars,data,par_step,chifun,nstep=50000):
    npar=len(pars)
    chain=np.zeros([nstep,npar])
    chivec=np.zeros(nstep)
    print("Starting MCMC")
    
    chi_cur=get_chisq(data,pars)
    for i in range(nstep):
        #print(pars)
        pars_trial=pars+np.random.randn(npar)*par_step
        while pars_trial[-1]<0: print('oop'); pars_trial[-1]=pars[-1]+np.random.randn()*par_step[-1]
        chi_trial=chifun(data,pars_trial)
        #we now have chi^2 at our current location
        #and chi^2 in our trial location. decide if we take the step
        accept_prob=np.exp(-0.5*(chi_trial-chi_cur))
        if np.random.rand()<accept_prob: #accept the step with appropriate probability
            pars=pars_trial
            chi_cur=chi_trial
            print('changed!')
        chain[i,:]=pars
        chivec[i]=chi_cur
    return chain,chivec
 

#%%
'''
[6.84756197e+01, 2.24459632e-02, 1.16026743e-01, 2.05024029e-09, 9.69172419e-01]
[2.37922333e+00, 5.35866653e-04, 5.29149674e-03, 3.91473939e-11, 1.35283042e-02]
'''
       
# load data
wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
xdata,ydata,yerror=wmap[:,0],wmap[:,1],wmap[:,2]
data=[xdata,ydata,yerror]

# initial guess from LM
params=np.array([6.84756197e+01, 2.24459632e-02, 1.16026743e-01, 2.05024029e-09, 9.69172419e-01, 0.05])
errors=np.array([2.37922333, 5.35866653e-04, 5.29149674e-03, 3.91473939e-11, 1.35283042e-02, 0.01])

# first run
chain1,chivec1=run_mcmc(params,data,errorsLM,get_chisq,nstep=10000)
params1=np.mean(chain1,axis=0)
errors1=np.std(chain1,axis=0)

# # improved
# chain2,chivec2=run_mcmc(params1,data,errors1,get_chisq,nstep=20)
# params2=np.mean(chain2,axis=0)
# errors2=np.std(chain2,axis=0)
# print('sigmas after improved chain ', errors2)

# # improved again
# pars_new=params+5*np.random.randn(len(params))*errors2
# chain3,chivec3=run_mcmc(pars_new,data,errors2,get_chisq,nstep=50)
# params3=np.mean(chain3,axis=0)
# errors3=np.std(chain3,axis=0)
# print("sigmas after improved again ", errors3)



#%%

for i in range(6):plt.plot(chain1[:,i])

#%%

for i in range(6):
    ffti=np.abs(np.fft.rfft(chain1[:,i]))
    plt.loglog(ffti/ffti.max(), label='Param %s'%i)
    
plt.legend()

#%%

print('final parameters with errors:')
for par,err in zip(params1,errors1): print('&=%.3g \\pm %.1g \\\\' %(par,err)) 









