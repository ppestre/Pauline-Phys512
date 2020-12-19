# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 14:49:51 2020

@author: Po
"""

import numpy  as np

def V(x,y): 
    if x==0 and y==0: return 1.
    return -np.log(x**2+y**2)/np.log(16)



#%%