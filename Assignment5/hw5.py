# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 15:04:58 2020

@author: Po
"""
import numpy as np
import matplotlib.pyplot as plt
import h5py
import json
from os.path import join 

from scipy.signal      import tukey
from scipy.signal      import welch
from scipy.interpolate import interp1d
from scipy.integrate   import cumtrapz as integral
from scipy.ndimage     import uniform_filter

# =============================================================================
# Data
# =============================================================================

# Data directory
path  = '.\\LOSC_Event_tutorial'

# Load dictionary of event info
events = json.load(open(join(path,'BBH_events_v3.json')))


#%% ===========================================================================
# Functions
# =============================================================================

#### Reading h5py ####

def read_template(filename):
    '''Read h5py template file'''
    dataFile=h5py.File(filename,'r')
    template=dataFile['template']
    th=template[0]
    tl=template[1]
    return th,tl

def read_file(filename): 
    '''Read h5py data file'''
    dataFile=h5py.File(filename,'r')
    meta=dataFile['meta']

    utc=meta['UTCstart'][()]
    duration=meta['Duration'][()]
    strain=dataFile['strain']['Strain'][()]
    dt=(1.0*duration)/len(strain)

    dataFile.close()
    return strain,dt,utc

#### Helpful ####

def window(n):
    '''Window of n points. Flat center, outter 1/8 drops off like cosine'''
    return tukey(n,0.125)

def psd(signal, dt):
    '''
    Takes psd with overlapping segments. Length divided by 4.
    '''
    freq, psd = welch(signal, 1/dt, window="blackman", nperseg= 4/dt)
    return psd, freq

def smooth(spectrum, freq, n=3):
    '''
    Smooths psd by averaging n points (moving)
    '''
    return interp1d(freq, uniform_filter(spectrum,n))

#### Filtering ####

def matched_filter(strain,template,N_func,dt):
    '''
    d is the windowed strain data
    A is the template
    N is the noise
    in FT space, m=conj(A)*d/N
    '''
    # Window 
    win=window(len(strain))
    
    # Frequencies used
    df = 1/(dt*len(strain))
    freq = np.fft.rfftfreq(len(strain),dt)
    N=N_func(freq)
    
    # FTs with window
    d_ft=np.fft.rfft(win*strain) 
    A_ft=np.fft.rfft(win*template)

    # Matched filter
    mf_ft=np.conj(A_ft)*(d_ft/N)
    mf=np.fft.irfft(mf_ft)
    
    # Integrate to find halfway point
    intg = integral(np.abs(mf_ft), dx=df, initial=0)
    mid_idx = np.argmin(np.abs(intg - max(intg)/2))

    sigsq_template = np.conj(A_ft)*A_ft/N # or is it np.conj(Aft)/(Aft*N)
    sigma = np.sqrt(np.abs(np.sum(sigsq_template)*df) )
    
    sig_data = np.std(np.abs(mf))
    max_data = np.max(np.abs(mf))
    SNR_data = (max_data-np.mean(np.abs(mf)))/sig_data
    
    # Normalize to sigma=1 and shift back to zero
    SNR=np.fft.ifftshift(mf)/sigma
    
    return SNR, SNR_data, freq[mid_idx]

# fig,ax=plt.subplots(figsize=(8,4))
# ax.plot(window(101))
# ax.set_title("Window")
# ax.set_xlabel("Number of Points")
# plt.tight_layout()
# plt.savefig('window')


#%% ===========================================================================
# Noise Model
# =============================================================================

# Empty arrays for noise spectra
Hspectra=[]
Lspectra=[]

# First loop over events to average noise spectra
for name in events:
    
    event=events[name]

    # Load data for Hanford and Livingston
    strainH,dt,utc=read_file(join(path,event['fn_H1']))
    strainL,dt,utc=read_file(join(path,event['fn_L1']))
    
    spec_H, freq = psd(strainH,dt)
    spec_L, freq = psd(strainL,dt)
    
    Hspectra.append(spec_H)
    Lspectra.append(spec_L)

N_H = smooth(np.mean(Hspectra,axis=0),freq,1)
N_L = smooth(np.mean(Lspectra,axis=0),freq,1)

f=np.linspace(10,1000,10000)
fig,ax=plt.subplots(figsize=(6,4))
ax.set_title("Noise Model")
ax.set_xlabel("Frequency")

ax.loglog(f,N_H(f), label='Hanford')
ax.loglog(f,N_L(f), label='Livingston')
plt.legend()
plt.tight_layout()
plt.savefig('noise')

#%% ===========================================================================
# Finding Gravitational Waves!
# =============================================================================

# Loop over events and filter
for name in ['GW150914', 'LVT151012', 'GW151226', 'GW170104']:#['GW150914', 'LVT151012', 'GW151226', 'GW170104']
    event=events[name]

    # Load data for Hanford and Livingston
    strainH,dtH,utc=read_file(join(path,event['fn_H1']))
    strainL,dtL,utc=read_file(join(path,event['fn_L1']))
    
    # Template
    template_name = join(path,event['fn_template'])
    tH,tL=read_template(template_name)

    # Time, same for both
    t=np.arange(0,dt*len(strainH),dt)
    
    # Do the matched filter!
    mf_H, SNR_data_H, mid_freq_H = matched_filter(strainH,tH,N_H,dt)
    mf_L, SNR_data_L, mid_freq_L = matched_filter(strainL,tL,N_L,dt)

    # Plots!
    fig,[ax0,ax1]= plt.subplots(1,2, figsize=(14,4), sharex=True)
    fig.suptitle('Event '+ name)
    
    # Strain
    ax0.set_title("Normalized Strain Data (A.U.)")
    ax0.set_xlabel("Time since %s (s)" %event['utcevent'])
    ax0.plot(t,strainH, label='H')
    ax0.plot(t,strainL, label='L')
    
    # Filter
    ax1.set_title("Matched filter output")
    ax1.set_xlabel("Time shift")
    ax1.plot(t,mf_H, label='H')
    ax1.plot(t,mf_L, label='L')
    
    # Clean up
    plt.legend()
    plt.tight_layout()
    plt.savefig(name)
    
    # print('Event '+ name+'\\\\')
    # print('Hanford expected SNR %.3f \\\\ ' % np.max(np.abs(mf_H)))
    # print('Hanford Data SNR %.3f \\\\' % SNR_data_H)
    # print()
    # print('Livingston expected SNR %.3f \\\\' % np.max(np.abs(mf_L)))
    # print('Livingston Data SNR %.3f \\\\' % SNR_data_L)
    # print() 
    
    #print('Event \quad'+name+'\qquad Hanford:\quad %.3f Hz \qquad Livingston:\quad  %.3f Hz\\\\ ' %(mid_freq_H, mid_freq_L))





