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


def ring_du(Q, R, z):
    '''
    Electric field at z from a ring of thickness du of a shell of radius R 
    with u=cos(theta) ranging from -1 to  1.
    '''
    return qk * (z*u-R) / (z**2 + R**2-2*R*z*u) ** (3/2)






def integrate_step(fun,x1,x2,tol):
    print('integrating from ',x1,' to ',x2)
    x=np.linspace(x1,x2,5)
    y=fun(x)
    area1=(x2-x1)*(y[0]+4*y[2]+y[4])/6
    area2=(x2-x1)*( y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12
    myerr=np.abs(area1-area2)
    if myerr<tol:
        return area2
    else:
        xm=0.5*(x1+x2)
        a1=integrate_step(fun,x1,xm,tol/2)
        a2=integrate_step(fun,xm,x2,tol/2)
        return a1+a2

#
        
ans=integrate_step(ring_du,-1,1,0.00001)

print('The field is %s *q*k_e '%ans)
