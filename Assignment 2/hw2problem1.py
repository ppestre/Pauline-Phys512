# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 12:34:08 2020

@author: Pauline

Assignment 2

"""

import numpy as np

def integrate_step_old(fun,x1,x2,tol):
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

def integrate_step(fun,x1,x2,tol):
    x=np.linspace(x1,x2,5)
    y=fun(x)
    area1=(x2-x1)*(y[0]+4*y[2]+y[4])/6
    area2=(x2-x1)*(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12
    myerr=np.abs(area1-area2)
    if myerr<tol:
        return area2
    else:
        a1=more_steps(fun,x[:3], y[:3],tol/2)
        a2=more_steps(fun,x[2:], y[2:],tol/2)
        return a1+a2
    
def more_steps(fun,x_old,y_old,tol):
    # this part could definitely be more elegant but oh well
    x0,x2,x4=x_old
    y0,y2,y4=y_old
    
    x1, x3 = (x0+x2)/2., (x2+x4)/2.
    y1, y3 = fun(x1), fun(x3)
    
    x=[x0,x1,x2,x3,x4]
    y=[y0,y1,y2,y3,y4]
    
    #print("saved+=3")
    
    area1=(x[-1]-x[0])*(y[0]+4*y[2]+y[4])/6
    area2=(x[-1]-x[0])*(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12
    myerr=np.abs(area1-area2)
    if myerr<tol:
        return area2
    else:
        a1=more_steps(fun,x[:3], y[:3],tol/2)
        a2=more_steps(fun,x[2:], y[2:],tol/2)
        return a1+a2
    
print('Integrating sin(x) from 0 to pi...') 
integrate_step(np.sin,0,np.pi,1e-5)
print('That saved 78 function calls! \n')



print('Integrating exp(x) from 0 to 1...') 
integrate_step(np.exp,0,1,1e-5)
print('That saved 18 function calls! \n')


print("""
Results:

%timeit(integrate_step(np.sin,0,np.pi,1e-5))
385 µs ± 30.5 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)

%timeit(integrate_step_old(np.sin,0,np.pi,1e-5))
443 µs ± 14.3 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)

%timeit(integrate_step(np.exp,0,1,1e-5))
131 µs ± 5.11 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)

%timeit(integrate_step_old(np.exp,0,1,1e-5))
181 µs ± 7.29 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)
""")
