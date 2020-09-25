# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:08:34 2020

@author: Pauline
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

#total charge times Coulombs constant
qk=1.0

#Radius of shell
R=1.0

tol=1e-2

def ring_uz(u, z):
    '''
    Electric field at z from a ring of thickness du of a shell of radius R 
    with u=cos(theta) ranging from -1 to  1.
    '''
    return qk * (z-R*u) / (z**2 + R**2-2*R*z*u)**(3/2)




def integrate_step(fun,x1,x2,tol):
    #print('integrating from ',x1,' to ',x2)
    x=np.linspace(x1,x2,5)
    y=fun(x)
    area1=(x2-x1)*(y[0]+4*y[2]+y[4])/6
    area2=(x2-x1)*(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12
    myerr=np.abs(area1-area2)
    if myerr<tol:
        return area2
    else:
        xm=0.5*(x1+x2)
        
        a1=integrate_step(fun,x1,xm,tol/2)
        a2=integrate_step(fun,xm,x2,tol/2)
        return a1+a2

    

z=np.linspace(0,5,1001)

E_mine=[]
for z0 in z:
    try:
        ring_z = lambda u: ring_uz(u,z0)
        E0=integrate_step(ring_z,-1,1,tol)
        E_mine.append(E0) 
    except RecursionError:
        E_mine.append(np.inf)
plt.plot(z,E_mine, '-', label='My integrator')

E_quad=[]    
for z0 in z:
    ring_z = lambda u: ring_uz(u,z0)
    E0=quad(ring_z,-1,1)[0]
    E_quad.append(E0)   
plt.plot(z,E_quad, '--', label='Scipy quad')


plt.xlabel('Distance from center, R=1m')
plt.ylabel('Electric field / $qk_e$')

plt.savefig('Efield.png')
plt.legend()


