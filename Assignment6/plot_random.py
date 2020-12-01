# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 12:20:33 2020

@author: Pauline
"""

import numpy as np
import matplotlib.pyplot as plt


x,y,z = np.genfromtxt('rand_points.txt', unpack=True)

ax = plt.subplot(121,projection='3d')
ax.plot(x,y,z,'.')
ax.set_title('Random points C')


x,y,z = np.random.rand(3,10000)

ax = plt.subplot(122,projection='3d')
ax.plot(x,y,z,'.')
ax.set_title('Random points Python')