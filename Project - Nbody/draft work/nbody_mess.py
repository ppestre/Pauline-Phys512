# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 14:41:47 2020

@author: Pauline
"""

import numpy as np
import matplotlib.pyplot as plt

          
class Space:
    def __init__(self,pos=None, vel=None,m=1.0,G=10.0,soft=0.3,dt=0.01,ncells=100,nparts=100000,dim=3,lims=((-2,2),(-2,2),(-2,2))):
        self.m = m
        self.G = G
        self.dt= dt
        self.soft = soft
        self.dim  = dim
        self.lims = lims

        # if nparts is None: nparts = ncells**dim
        self.nparts = nparts
        
        # if ncells is None: ncells = int(nparts**(1/dim))
        self.ncells = ncells
            
        if pos is None: pos = np.random.uniform(*lims[0],(nparts,dim)) #watch lims
        self.positions  = pos
        
        if vel is None: vel = np.zeros((nparts,dim))
        self.velocities = vel
        
        self.x = self.y = self.z = np.linspace(*lims,ncells)
        
        self.grid = np.meshgrid(self.x,self.y,self.z)
        
        
        self.calc_forces()
        
        # self.density    = self.calc_density()
        # self.Greens     = self.calc_Greens()
        # self.potential  = self.calc_potential()
        print('initialized')

    def calc_density(self, pos=None):
        if pos is None: pos=self.positions
        ranges=(self.lims,)*3
        self.density, self.edges = np.histogramdd(pos, bins=self.ncells, range=ranges)
        return self.density

    def calc_Greens(self):
        r = np.sqrt(np.sum(np.square(self.grid), axis=0))
        self.Greens = self.G/(r+self.soft)
        return self.Greens
    
    def calc_potential(self,pos=None): # return convolution of density with 1/(r+soft)
        if pos is None: pos=self.positions
        d = self.calc_density(pos)
        g = self.calc_Greens()
        self.potential = np.abs(np.fft.irfftn(np.fft.rfftn(d) + np.fft.rfftn(g)))
        return self.potential
    
    # def grad(self, v ,bc='wrap'):
    #     h = ( max(self.lims) - min(self.lims) ) / self.ncells
    #     if bc=='wrap':
            
    #         grad = np.zeros((3,*np.shape(v)))
    #         n=len(v)
            
    #         for i in range(n):
    #             for j in range(n):
    #                 for k in range(n):
    #                     grad[0,i,j,k]=(v[i,j,k]-v[i-1,j,k])/(2*h)
    #                     grad[1,i,j,k]=(v[i,j,k]-v[i,j-1,k])/(2*h)
    #                     grad[2,i,j,k]=(v[i,j,k]-v[i,j,k-1])/(2*h)
        
    #     return grad
            
            
    
    def calc_forces(self,pos=None): #take derivative of potential, gradV=F=ma
        if pos is None: pos=self.positions
        self.forces = -np.array(np.gradient(self.calc_potential(pos)))
        return self.forces
    
    def to_grid(self,r, a_grid):
        span= max(self.lims) - min(self.lims)
        i,j,k = (np.round(r*ncells/span)).astype(int)
        return a_grid[:,i,j,k]
    
    def update(self):
        dt= self.dt
        
        rs_mid = self.positions + 0.5*dt*self.velocities
        
        a_grid = self.calc_forces(rs_mid)/self.m
        
        min, max = self.lims
        
        for i,(r,v) in enumerate(zip(self.positions,self.velocities)):
            if r.min()>=min and r.max()<=max:
                a = self.to_grid(r,a_grid)
                v_mid = v + 0.5*a*dt
                self.positions[i]  = r + dt*v_mid
                self.velocities[i] = v + dt*a
                
        self.calc_density()
        # print('updated')
        return 
    
        


    
if __name__=='__main__':
        
    nparts = 1000
    ncells = 50
    lims = (-1.,1.)
    
    # p = Space(nparts=nparts, ncells=ncells, lims=lims)
    
    # p  =  Space(pos=np.array([[0.,0.,0.]]),
    #             vel=np.array([[0.,0.,0.]]),
    #             nparts=1, ncells=ncells, lims=lims)
    
    p  =  Space(pos=np.array([[-.2,0.,0.],[.2,0.,0.]]),
                vel=np.array([[0.,-.2,0.],[0.,0.2,0.]]),
                nparts=2, ncells=ncells, lims=lims)
    
    
    # x=y=z= np.linspace(*p.lims,p.ncells)
    # grid = np.meshgrid(x,y,z)
        

    
    fig,ax = plt.subplots(figsize=(4,4), dpi=100 )
    #c_ax=fig.add_subplot()
    # ax.set_xlabel("y")
    # ax.set_ylabel("x")
    ax.set(xlim=lims, ylim=lims)
    # c=ax.imshow(p.density[2,:,:],origin="lower", aspect="auto", extent=(*lims,*lims))
    # fig.colorbar(c)
    
    plt.tight_layout()
    
    numplots=100
    numsteps=30
    
    for i in range(numplots):
        print(p.positions)
        ax.plot(p.positions[:,0], p.positions[:,1],'o')
        # ax.imshow(p.density.sum(axis=2),origin="lower", aspect="auto", cmap='viridis', interpolation='none',
        #           extent=(*lims,*lims))
        # ax.quiver(x,y,p.forces[2,2,:,:],p.forces[1,2,:,:])#.sum(axis=2)
        plt.show()
        for j in range(numsteps):
            p.update()
        plt.pause(.1)
        print('step')
            
# ncells=5  
# for i in range(ncells):
#     for j in range(ncells):
#         for k in range(ncells):
#             plt.plot(p.forces[0,i,j,k])  
    
