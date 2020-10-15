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
    if len(pars)==5:
        H0,ombh2,omch2,As,ns=pars
        tau=0.05
    elif len(pars)==1:
        H0,ombh2,omch2,As,ns=[68.4756197, .0224459632, .116026743, 2.05024029e-09, .969172419]
        tau=pars[0]
    else:
        H0,ombh2,omch2,As,ns,tau=pars
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
    print('Start derivative...')
    # for ds a list of parameters to vary and their indices
    grad = np.zeros((len(x),len(ds)))
    
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
        
    print('...derivative done\n')
    return grad


#%%

# Load data
wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
xdata,ydata,yerror=wmap[:,0],wmap[:,1],wmap[:,2]

N=np.diag(yerror**2)
Ninv=np.diag(1/yerror**2)

#tau=0.05
pars=np.array([65,0.02,0.1,2e-9,0.96])
ymodel=A(xdata, pars)

resid=ydata-ymodel #data minus current model
chisq=np.sum(resid**2/yerror**2)

# parameters to vary
ds = list(zip([0,1,2,3,4], np.array([65,0.02,0.1,2e-9,0.96])/70.))


print("Initial guess: %s" % pars)
print("Gives 𝛘² = %.3f\n" % chisq)

lm = 0.001
done = False
small = False

for i in range(100):
    
    #This step is slow so I pulled it out of the internal loop, doesn't change
    dAdm=get_grad(xdata,pars, ds)

        
    # Perform next step
    increase = True
    # If lambda was increased
    while increase:
        print("Lambda set to: %e" % lm)
        
        gradchi=dAdm.T@(Ninv@resid)
        curv=dAdm.T@Ninv@dAdm
        step=np.linalg.inv(curv+lm*np.diag(np.diag(curv)))@gradchi
        
        parsnew=pars+step
        ynew=A(xdata, parsnew)
        chisqnew=np.sum((ydata-ynew)**2/yerror**2)
        print("𝛘² = %.3f\n" % chisqnew)
        # Increase lambda is chi squared got worse
        if chisqnew >= chisq: 
            lm *= 10
            print(chisqnew,chisq)
            print('No improvement, increasing lamba')
        else:
            print('Better, decreasing lamba')
            lm *= 0.1
            increase = False
            # If the change in chi was small flag for end
            if chisq - chisqnew < 1e-4:
                # If decrease was small twice, consider converged
                if small:
                    done = True
                    print("Converged!")
                # Otherwise, indicate that it was small once
                else:
                    small = True
                    print("Small change, one more try")
            # Otherwise, reset small flag
            else:
                small = False   
            # Take new values
            ymodel=ynew
            pars=parsnew
            chisq=chisqnew
    # Check end condition
    if done: break
 

print()
print("Final answer: %s" % pars)
print("Gives 𝛘² = %.3f\n" % chisq)

    
#%% 
    
    
pars_errs=np.sqrt(np.diag(np.linalg.inv(curv)))
#print('final parameters are ',pars,' with errors ',pars_errs)

print('final parameters with errors:')
for par,err in zip(pars,pars_errs): print('%.3g +- %.1g' %(par,err)) 


#x=np.linspace(0,2000,2000)
x=xdata

plt.clf();

plt.errorbar(wmap[:,0],wmap[:,1],wmap[:,2],zorder=0)
#plt.plot(xdata,ydata,'.',label='Data')
plt.plot(x, A(x, pars),'r', label='Spectrum model',zorder=10)


#%%


for div in (50.,70.,100.):
    grad=get_grad(xdata,pars, list(zip([0,1,2,3,4], np.array([65,0.02,0.1,2e-9,0.96])/div)))
    plt.plot(xdata,grad[:,0])

