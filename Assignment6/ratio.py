# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 11:01:24 2020

@author: Pauline
"""

import numpy as np
import matplotlib.pyplot as plt

n=100000

u=np.random.rand(n)
v=np.random.rand(n)*2/np.e

ratio = v/u
accept= u < np.sqrt(np.exp(-ratio))

# From the ratio of uniforms
x_ratio = ratio[accept]

# From transformation
x_trnsf = -np.log(np.random.rand(len(x_ratio)))


dist,bins=np.histogram(x_ratio,200)

y=dist/max(dist)
x=(bins[1:]+bins[:-1])/2

fig,[ax1,ax2]=plt.subplots(1,2, figsize=(10,4))


ax1.plot(x,y,'.', label='Ratio of Uniforms')
ax1.plot(x,np.exp(-x), label='Analytic')
ax1.set_title('PDF')
ax1.legend()

ax2.plot(np.sort(x_ratio), label='Ratio of Uniforms')
ax2.plot(np.sort(x_trnsf), label='Transformation')
ax2.set_title('CDF')
ax2.legend()

plt.tight_layout()
plt.savefig('ratio.png')


print('Accepted fraction was %.4f'  %(len(x_ratio)/(2*n)) )
