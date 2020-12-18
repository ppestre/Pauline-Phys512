# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 14:41:47 2020

@author: Pauline
"""

import numpy as np
import matplotlib.pyplot as plt

          
class Space:
    def __init__(self, pos, vel, m=1.0, dt=0.01, G=.1,ncells=100, lims=None):
        try: self.m=m.copy()
        except: self.m=m
        
        self.positions  = pos.copy()
        self.velocities = vel.copy()
        self.dt= dt
        self.G = G
        
        if lims is None: print('oopsie')

        self.lims = lims
        self.x = np.linspace(*lims[0],ncells)
        self.y = np.linspace(*lims[1],ncells)
        self.z = np.linspace(*lims[2],ncells)
        self.grid = np.meshgrid(self.x,self.y,self.z)
        
        self.centergrid = np.meshgrid(self.x - np.mean(self.x),
                                      self.y - np.mean(self.y),
                                      self.z - np.mean(self.z))
        
        self.dx = self.x[1] - self.x[0]

        self.ncells = ncells
        self.nparts = len(pos)
        
        self.calc_Greens() # do once

        self.calc_forces() # includes density and potential
        
        print('Initialized!')

    def calc_density(self, pos=None):
        if pos is None: pos=self.positions
        self.density, self.edges = np.histogramdd(pos, bins=self.ncells, range=self.lims)
        return self.density

    def calc_Greens(self, soft=0.1):
        r = np.sqrt(np.sum(np.square(self.centergrid), axis=0))
        r[r<soft]=soft
        self.Greens = -self.G/r 
        return self.Greens
    
    def calc_potential(self,pos=None): # return convolution of density with 1/(r+soft)
        if pos is None: pos=self.positions
        dFT = np.fft.fftn(self.calc_density(pos))
        gFT = np.fft.fftn(np.fft.fftshift(self.Greens))
        self.potential = np.roll(np.real(np.fft.ifftn(dFT * gFT)),1,(0,1,2))
        return self.potential

    def calc_forces(self,pos=None): #take derivative of potential, -gradV=F=ma
        if pos is None: pos=self.positions
        self.forces = -self.grad(self.calc_potential(pos)) # np.array(np.gradient(self.calc_potential(pos)))
        return self.forces
    
    def to_grid(self,r): 

        # r0 =  r - np.mean(lims,axis=1) # this is the center! make the corner
        r0=np.array(self.lims)[:,0]
        span = np.diff(self.lims[0])
        
        i,j,k = ((r-r0)*(ncells)/span).astype(int)
        return i,j,k
    
    def grad(self,pot):
        gx = (np.roll(pot, -1, axis=0) - np.roll(pot, 1, axis=0)) / (2 * self.dx) 
        gy = (np.roll(pot, -1, axis=1) - np.roll(pot, 1, axis=1)) / (2 * self.dx)
        gz = (np.roll(pot, -1, axis=2) - np.roll(pot, 1, axis=2)) / (2 * self.dx)
        return np.array([gx,gy,gz])
    
    def is_in(self, r):
        xmin,xmax = self.lims[0]
        ymin,ymax = self.lims[1]
        zmin,zmax = self.lims[2]

        return ((xmin <= r[0] < xmax) and
                (ymin <= r[1] < ymax) and 
                (zmin <= r[2] < zmax))
    
    def update(self):
        dt= self.dt
        
        rs_mid = self.positions + 0.5*dt*self.velocities
        a_grid = self.calc_forces(rs_mid)/self.m
        
        for n,(r,v) in enumerate(zip(self.positions,self.velocities)):
            if self.is_in(r):
                i,j,k = self.to_grid(r)
                a = a_grid[:,i,j,k]
                v_mid = v + 0.5*a*dt
                self.positions[n]  = r + dt*v_mid
                self.velocities[n] = v + dt*a
                
        self.calc_density()
        # print('updated')
        return 
 
    
 
    
 
nparts = 1000
ncells = 51
# lims = ((0,1),(1,2),(3,4))
lims = ((0,2),(0,2),(0,2))

start= 'two'

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


numplots=200
numsteps=5
loop = 1

plt.xlim(.5,1.5)
plt.ylim(.5,1.5)
plt.xlabel('y')
plt.ylabel('x')

if loop:
    for i in range(numplots):
        # print(p.positions)
        # plotxyplane(p.density,z=25)
        # plotxysum(p.density)
        
        plt.plot(p.positions[:,1], p.positions[:,0],'o')
        
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
