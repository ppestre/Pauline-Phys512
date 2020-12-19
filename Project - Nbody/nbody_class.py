# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 14:41:47 2020

@author: Pauline
"""

import numpy as np
# from numba import jit
          
class Space:
    def __init__(self, pos, vel, lims, m=1.0, dt=0.01, G=0.1, ncells=100, periodic=True):
        
        try: self.m=m.copy()
        except: self.m=m
        
        self.positions  = pos.copy()
        self.velocities = vel.copy()
        self.dt= dt
        self.G = G
        self.periodic = periodic
        
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

        self.calc_accel() # includes density and potential
        
        print('Initialized!')

    def calc_density(self, pos=None):
        if pos is None: pos=self.positions
        if hasattr(self.m, "__len__"): 
            self.density,  self.edges  = np.histogramdd(pos, bins=self.ncells,
                                                        range=self.lims, 
                                                        weights=self.m)
        else: self.density, self.edges = np.histogramdd(pos, bins=self.ncells,
                                                        range=self.lims)
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
        if self.periodic:
            pot = np.real(np.fft.ifftn(dFT * gFT))
        else:
            dFTpad = np.pad(dFT,((0, self.ncells),))
            gFTpad = np.pad(dFT,((0, self.ncells),))
            potpad = np.real(np.fft.ifftn(dFTpad * gFTpad))
            pot = potpad[:self.ncells,:self.ncells,:self.ncells]
        self.potential = np.roll(pot,1,(0,1,2)) #mysterious off-by-one?
        return self.potential

    def calc_accel(self,pos=None): #take derivative of potential, -gradU=F=ma,-gradV=F=m
        if pos is None: pos=self.positions
        self.accel = -self.grad(self.calc_potential(pos)) 
        # -np.array(np.gradient(self.calc_potential(pos)))
        return self.accel
    
    def to_grid(self,r): 

        r0=np.array(self.lims)[:,0]
        span = np.diff(self.lims[0])
        
        i,j,k = ((r-r0)*(self.ncells)/span).astype(int)
        return i,j,k
    
    def grad(self,pot): #tried using np.gradient but was getting some indexing bugs
        gx = (np.roll(pot, -1, axis=0) - np.roll(pot, 1, axis=0)) / (2 * self.dx) 
        gy = (np.roll(pot, -1, axis=1) - np.roll(pot, 1, axis=1)) / (2 * self.dx)
        gz = (np.roll(pot, -1, axis=2) - np.roll(pot, 1, axis=2)) / (2 * self.dx)
        return np.array([gx,gy,gz])
    
    def is_in(self, i,j,k):
        return ((0 <= i < self.ncells) and
                (0 <= j < self.ncells) and 
                (0 <= k < self.ncells))
    
    def update(self):
        
        dt= self.dt
        n = self.ncells
        
        rs_mid = self.positions + 0.5*dt*self.velocities
        a_grid = self.calc_accel(rs_mid) # returns F/m
        
        for l,(r,v) in enumerate(zip(self.positions,self.velocities)):
            i,j,k = self.to_grid(r)
            if self.periodic:
                i,j,k = i%n,j%n,k%n
                a = a_grid[:,i,j,k]
                v_mid = v + 0.5*a*dt
                self.positions[l]  = r + dt*v_mid
                self.velocities[l] = v + dt*a
            elif self.is_in(i,j,k):
                a = a_grid[:,i,j,k]
                v_mid = v + 0.5*a*dt
                self.positions[l]  = r + dt*v_mid
                self.velocities[l] = v + dt*a
        
        self.Epot = np.sum(self.potential) # if mass=1
        self.Ekin = 0.5*np.sum(np.square(self.potential))         
                
        self.calc_density()
        # print('updated')
        return 
