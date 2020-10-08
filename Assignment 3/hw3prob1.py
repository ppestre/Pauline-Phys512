# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 18:58:07 2020

@author: Pauline
"""

import numpy as np
import matplotlib.pyplot as plt

data=np.genfromtxt('dish_zenith.txt',unpack=True) # columnwise not rowwise

x,y,z=data[0],data[1],data[2]

A=[np.ones(len(x)), x, y, x**2+y**2]

Ainv=np.linalg.pinv(A)

b=np.dot(np.transpose(Ainv), z)

z_fit=np.dot(np.transpose(A), b)

plt.plot(x,z,'.',label='Data')
plt.plot(x,z_fit,'.', label='Fit')
