# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 14:30:41 2020

@author: Pauline
"""

import numpy as np
import matplotlib.pyplot as plt

from nbody_class import Space
  
 
nparts = 1000
ncells = 51
# lims = ((0,1),(1,2),(3,4))
lims = ((0,2),(0,2),(0,2))

start= 'many'

if start=='one':
    p  =  Space(pos=np.array([[.8,1.4,3.6]]),
                vel=np.array([[0.,0.,0.]]),
                ncells=ncells, lims=lims)

elif start=='two':
    p  =  Space(pos=np.array([[.8,1,1],[1.2,1,1]]),
                vel=np.array([[0.,.339,0.],[0.,-.339,0.]]),
                ncells=ncells, lims=lims)
    
elif start=='many':
    pos = np.random.uniform(*lims[0],(nparts,3))
    vel = np.zeros((nparts,3))
    p = Space(pos,vel, ncells=ncells, lims=lims)


#I think i'm plotting wrong??

def plotxysum(field,axis=2,lims=lims,**kwargs):
    plt.imshow(field.sum(axis=axis), 
               origin='lower', aspect="equal",  
               cmap='viridis', interpolation='none',
               extent=(*lims[1],*lims[0]), **kwargs)

def plotxyplane(field,z=0,lims=lims,**kwargs):
    plt.imshow(field[:,:,z], 
               origin='lower', aspect="equal", 
               cmap='viridis', interpolation='none',
               extent=(*lims[1],*lims[0]), **kwargs)
    plt.xlabel('y')
    plt.ylabel('x')
    
def plotyzplane(field,x=0,lims=lims,**kwargs):
    plt.imshow(field[x,:,:],
               origin='lower', aspect="equal", 
               cmap='viridis', interpolation='none',
               extent=(*lims[2],*lims[1]), **kwargs)
    plt.xlabel('z')
    plt.ylabel('y')

def plotzxplane(field,y=0,lims=lims,**kwargs): #this one works
    plt.imshow(field[:,y,:], 
               origin='lower', aspect="equal", 
               cmap='viridis', interpolation='none',
               extent=(*lims[2],*lims[0]), **kwargs)
    plt.xlabel('z')
    plt.ylabel('x')


numplots=20
numsteps=5
loop = 1

plt.xlim(.5,1.5)
plt.ylim(.5,1.5)
plt.xlabel('y')
plt.ylabel('x')

if loop:
    for i in range(numplots):
        # print(p.positions)
        plotxyplane(p.density,z=25)
        # plotxysum(p.density)
        
        # plt.plot(p.positions[:,1], p.positions[:,0],'o')
        
        # plt.quiver(x,y,p.forces[2,2,:,:],p.forces[1,2,:,:])#.sum(axis=2)
        
        plt.show()
        
        for j in range(numsteps):
            p.update()
        plt.pause(.01)
        print('step')

else:
    plotzxplane(p.potential,2)
    
    plt.plot(p.positions[:,1], p.positions[:,0],'o')
    
    # plt.quiver(p.x,p.y,p.forces[1,:,:,2],p.forces[0,:,:,2])#.sum(axis=2)
    
    # plt.quiver(p.z,p.x,p.forces[2,:,25,:],p.forces[0,:,25,:])
    
    # plt.quiver(p.z,p.y,p.forces[2,25,:,:],p.forces[1,25,:,:])
