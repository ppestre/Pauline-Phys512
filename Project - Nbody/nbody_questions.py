# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 14:30:41 2020

@author: Pauline
"""

import numpy as np
import matplotlib.pyplot as plt
from nbody_class import Space
from celluloid import Camera


def plot3(field, lims, slice=(0,0,0), sum=False, figaxes=None):
    
    if figaxes is None: fig, axes = plt.subplots(1,3, figsize=(10,3))
    else: fig,axes = figaxes
    axyx, axzy, axzx = axes
    
    axyx.set_xlabel('y')
    axyx.set_ylabel('x')
    axzy.set_xlabel('z')
    axzy.set_ylabel('y')
    axzx.set_xlabel('z')
    axzx.set_ylabel('x')
    
    if sum:
        yx = field.sum(axis=2)
        zx = field.sum(axis=1)
        zy = field.sum(axis=0)

    else:
        x,y,z=slice
        yx = field[:,:,z]
        zx = field[:,y,:]
        zy = field[x,:,:]
    
    axyx.imshow(yx, extent=(*lims[1],*lims[0]),
                origin='lower', aspect="equal", interpolation='none') 
    axzy.imshow(zy, extent=(*lims[2],*lims[1]),
                origin='lower', aspect="equal", interpolation='none')
    axzx.imshow(zx, extent=(*lims[2],*lims[0]),
                origin='lower', aspect="equal", interpolation='none')
    
    return fig, axes


# Change index q to run different parts
q = 3

if q==1:
    
    # Set limits in x,y,z
    lims = ((-1,1),(-1,1),(-1,1))
    # Number of cells on each axis
    ncells = 51
    
    # Single particle at origin, at rest
    rs = [[0.,0.,0.,]]
    vs = [[0.,0.,0.,]]
    m  =  1
    
    # Create an instance of Space
    p = Space(pos=np.array(rs),
              vel=np.array(vs),
              ncells=ncells, 
              lims=lims)
       
    numplots=10
    numsteps=30
    looping = 1
    
    fig = plt.figure()
    camera = Camera(fig)
    plt.xlim(lims[1])
    plt.ylim(lims[0])
    
    
    if looping:
        for i in range(numplots):
            plt.plot(p.positions[:,1], p.positions[:,0],'o')
            plt.pause(.01)
            camera.snap()
            for j in range(numsteps):
                p.update()
    
    # animation = camera.animate()
    # animation.save('single_particle.gif')

    
if q==2:
    
    # Set limits in x,y,z
    lims = ((-1,1),(-1,1),(-1,1))
    # Number of cells on each axis
    ncells = 51
    
    # Two particles moving in opposite directions
    rs = [[-.2,0.,0.],[.2,0.,0.]]
    vs = [[0.,-.35,0.],[0.,.35,0.]]
    
    # Create an instance of Space
    p = Space(pos=np.array(rs),
              vel=np.array(vs),
              ncells=ncells, 
              lims=lims)
    
    # Make some plots!
    
    numplots=25
    numsteps=10
    looping = 1
    
    fig = plt.figure(figsize=(5,5))
    camera = Camera(fig)
    plt.xlim(lims[1])
    plt.ylim(lims[0])
    
    if looping:
        for j in range(numplots):
            plt.plot(p.positions[:,1], p.positions[:,0],'o')
            plt.pause(.01)
            for k in range(numsteps):
                p.update()
            camera.snap()
        # plt.savefig('circular_orbit.png')
        # animation = camera.animate()
        # animation.save('circular_orbit.gif')

if q==3:
    
    # Set limits in x,y,z
    lims = ((-20,20),(-20,20),(-20,20))
    # Number of cells on each axis
    ncells = 51
    # Number of particles
    nparts = 200000
    # Boundary conditions
    periodic = True
    
    # Two particles moving in opposite directions
    rs = np.random.uniform(*lims[0],(nparts,3))
    vs = np.zeros((nparts,3))
    
    # Create an instance of Space
    p = Space(pos=np.array(rs),
              vel=np.array(vs),
              ncells=ncells, 
              lims=lims)

    fig,axes = plot3(p.density, lims, sum=True)
    plt.tight_layout()
    # Make some plots!
    
    numplots=25
    numsteps=30
    looping = 1

    camera = Camera(fig)
    
    Epot=[]
    Ekin=[]
    
    if looping:
        for j in range(numplots):
            plot3(p.density, lims, sum=True, figaxes=(fig,axes))
            #plt.pause(.01)
            for k in range(numsteps):
                p.update()
            camera.snap()
            print(j)
            Epot.append(p.Epot)
            Ekin.append(p.Ekin)

        animation = camera.animate()
        #animation.save('Periodic%s_%iparts_%ix%isteps.gif'%(periodic,nparts,numplots,numsteps))
        
    if looping:
        for j in range(numplots):
            plot3(p.density, lims, sum=True, figaxes=(fig,axes))
            #plt.pause(.01)
            for k in range(numsteps):
                p.update()
            camera.snap()
            print(j)
            Epot.append(p.Epot)
            Ekin.append(p.Ekin)

        animation = camera.animate()
        #animation.save('Periodic%s_%iparts_%ix%isteps.gif'%(periodic,nparts,2*numplots,numsteps))

        print(np.diff(Epot)/np.mean(Epot))
        print(np.diff(Ekin)/np.mean(Ekin))
        

if q==4:
    
    # Set limits in x,y,z
    lims = ((-20,20),(-20,20),(-20,20))
    # Number of cells on each axis
    ncells = 25
    
    # Make a grid in kspace, then calculate 1/k^3 at each point.
    # There are smarter ways to do this but shhhhh...
    x = np.linspace(*lims[0],ncells)
    y = np.linspace(*lims[1],ncells)
    z = np.linspace(*lims[2],ncells)
    kgrid = np.meshgrid(2*np.pi/x,2*np.pi/y,2*np.pi/z)
    gridk = np.sqrt(np.sum(np.square(kgrid), axis=0))
    invk3 = 1/np.power(gridk,3)
    
    print('Please ignore divide by zero warning :P')
    

    
    
    
    noise = np.random.normal(size=invk3.shape)
    
    density = np.abs(np.fft.ifftn(invk3*np.fft.fftn(noise)))
    

    rs = np.array([(xval,yval,zval) for xval in x for yval in y for zval in z])
    ms = density.flatten()
    vs = np.zeros(rs.shape)     
    nparts = len(rs)    
                
    # Create an instance of Space
    p = Space(pos = rs, vel = vs, m = ms,
              ncells=ncells, 
              lims=lims)
    
    fig,axes = plot3(p.density, lims, sum=True)
    plt.tight_layout()
    
    # Make some plots!
    
    numplots=50
    numsteps=15
    looping = 1

    camera = Camera(fig)

    if looping:
        for j in range(numplots):
            plot3(p.density, lims, sum=True, figaxes=(fig,axes))
            plt.pause(.01)
            for k in range(numsteps):
                p.update()
            print(j)
            camera.snap()


        animation = camera.animate()
        animation.save('k3_%iparts_%ix%isteps.gif'%(nparts,numplots,numsteps))

    