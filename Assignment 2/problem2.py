# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 16:24:53 2020

@author: Pauline
"""

import numpy as np
from matplotlib import pyplot as plt

def cheb_mat_uniform(nx,ord):
    x=np.linspace(-1,1,nx)
    mat=np.zeros([nx,ord+1])
    mat[:,0]=1.0
    if ord>0:
        mat[:,1]=x
    if ord>1:
        for i in range(1,ord):
            mat[:,i+1]=2*x*mat[:,i]-mat[:,i-1]
    return mat,x


def legmat(n,npt):
    assert(npt>0)
    assert(n>=0)
    
    x=np.linspace(-1,1,npt)
    mat=np.zeros([n+1,npt])
    mat[0,:]=1.0
    if (n>0):
        mat[1,:]=x
    if (n>1):
        for i in range(1,n):
            mat[i+1,:]= ((2*i+1)*x*mat[i,:]-i*mat[i-1,:])/(i+1.0)
    return mat

n=100
ord=10
mat,x=cheb_mat_uniform(n,ord)
xnew=(x+3)/4
y=np.log2(xnew)
# This function over the interval -1 to 1 has the same values as log2(x) from 0.5 to 1




# Chebyshev fit
coeff_c=np.dot(np.linalg.pinv(mat),y)
cheb_fit=np.dot(mat,coeff_c)


coeff_l=np.polynomial.legendre.legfit(x,y,ord)
leg_fit=np.dot(np.transpose(legmat(ord,n)),coeff_l)

plt.clf();plt.plot(xnew,cheb_fit-y);plt.plot(xnew,leg_fit-y);
plt.legend(['Chebyshev Residual','Legendre Residual'])
plt.savefig('cheb_resids.png')

print('RMS error for Chebyshev is ',np.sqrt(np.mean((cheb_fit-y)**2)),' with max error ',np.max(np.abs(cheb_fit-y)))
print('RMS error for Legendre is ',np.sqrt(np.mean((leg_fit-y)**2)),' with max error ',np.max(np.abs(leg_fit-y)))

