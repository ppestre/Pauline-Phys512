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
        
# load data
wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
xdata,ydata,yerror=wmap[:,0],wmap[:,1],wmap[:,2]
data=[xdata,ydata,yerror]

# initial guess
pars=np.array([65,0.02,0.1,2e-9,0.96, 0.05])
pars_sigs0=np.array([1,0.0005,0.005,1e-10,0.01,0.01])

# first run
chain1,chivec1=run_mcmc(pars,data,pars_sigs0,get_chisq,nstep=200x0)
# pars_sigs1=np.std(chain1,axis=0)

# # improved
# chain2,chivec2=run_mcmc(pars,data,pars_sigs1,get_chisq,nstep=20)
# pars_sigs2=np.std(chain2,axis=0)
# print('sigmas after improved chain ', pars_sigs2)

# # improved again
# pars_new=pars+5*np.random.randn(len(pars))*pars_sigs2
# chain3,chivec3=run_mcmc(pars_new,data,pars_sigs2,get_chisq,nstep=50)
# pars_sigs3=np.std(chain2,axis=0)
# print("sigmas after improved again ", pars_sigs3)


plt.plot(chain1[:,0])
plt.plot(chain1[:,1])
plt.plot(chain1[:,2])
plt.plot(chain1[:,3])
plt.plot(chain1[:,4])
plt.plot(chain1[:,5])

#%%




# chains = monte(params, steps)
# params = np.mean(chains,axis=0)
# errors = np.std(chains, axis=0)
# cmb = get_spectrum(ls, params)
# chisq = chisqr(wmap[:,1], cmb, wmap[:,2])
# fig, ax = plt.subplots(2,1, sharex=True)
# res_ax = ax[0]
# dat_ax = ax[1]

# plot_res(res_ax, wmap[:,0], studentized(wmap[:,1], cmb, wmap[:,2]))
# plot_dat(dat_ax, wmap[:,0], wmap[:,1], wmap[:,2])
# plot_func(dat_ax, ls, cmb, label="fit")
# dat_ax.legend()

# print("𝛘² = %.3f" % chisq)
# plt.show(block=True)





