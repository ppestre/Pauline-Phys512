# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 14:41:47 2020

@author: Pauline
"""

import numpy as np
from numba import jit
          
class Space:
    def __init__(self, pos, vel, m=1.0, dt=0.01, G=0.1,ncells=100, lims=None):
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
        
        i,j,k = ((r-r0)*(self.ncells)/span).astype(int)
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
