# Scientific Computing
# Set1, ex. D
# implement a 2-dimensional diffusion equation

# http://www.timteatro.net/2010/10/29/performance-python-solving-the-2d-diffusion-equation-with-numpy/

import random
#import simpy
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, sin, pi
import timeit

nx = 3 # no of interval in x-direction
ny = 3 # # no of interval in x-direction
deltaX = 1/nx # length of an interval
deltaY = deltaX
tstep = 10 # number of time steps
tmax = 10 #max. simulation time
deltaT = tmax/tstep # length of time interval
d = 1 # diffusion coefficient

u = np.zeros((nx, ny)) # 2D array for x-andy-values at current time
ui = np.zeros((nx, ny)) # 2D array for x-andy-values at current time


	
	
# define boundary conditions
for i in range(0, nx):
	for j in range(0, ny):
		ui[i, 0] = 0 # bottom horizontal (x, y=1)
		ui[i, 1] = 1 # top horizontal (x, y=0)
		ui[0, j] = ui[1, j] # periodic bounday conditions along x-axis
		
		
# define initial conditions
# since ui was already initialized with zeros and intial conditions are zero, we don't have to fo anything


# calculate u from ui, calculate Laplacians
for k in range(1, tstep):
	for i in range(1, nx-1):
		for j in range(1, ny-1):
			uxx = (ui[i+1, j] + ui[i-1, j] - 2*ui[i,j]) / (deltaX)**2 # second derivative for x-values
			uyy = (ui[i, j+1] + ui[i, j-1] - 2*ui[i,j]) / (deltaY)**2 # second derivative for y-values
			
			u[i,j] = deltaT * d *(uxx + uyy)
	print (u)
	
