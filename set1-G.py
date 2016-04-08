# Scientific Computing
# Set1, ex. E
# Jacobi iteration, Gauss-seidel iteration, Successive Over-Relaxationanalytic

import random
#import simpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#import matplotlib.image as mpimg
from math import sqrt, sin, pi
#from pylab import *
#import scipy as sp
#from PIL import Image
#import timeit


nx = 50 # no of interval in x-direction
ny = 50 # # no of interval in x-direction
deltaX = 1/nx # length of an interval
deltaX2 = deltaX**2 
deltaY = deltaX
tmax = 10.0 #max. simulation time
d = 1.0 # diffusion coefficient
deltaT = (deltaX2 * 4 * d) / 10  # The max of delta t is deltaX2/(4*d), just to make sure we divide by 10
steps = int(tmax / deltaT) + 1 # calculate the number of time steps
omega = 1.5 #needed for successive over-relaxation method
delta = 0.1
tolerance = 10**(-5)

"""
#----------------START: JACOBI ITERATION-----------------------------------------
u = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at current time
ui = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at previous time
#itListJ = []

for j in range(0, ny+1):
	ui[0, j] = 1

# calculate u from ui, calculate Laplacians
while(tolerance < delta):
	delta = 0
	for j in range(0, ny+1):
		u[0, j] = 1
		for i in range(1, nx):
			u[i,j] = (ui[i-1,j] + ui[i,j-1] + ui[i+1,j] + ui[i,(j+1)%(ny+1)]) / 4
	delta = np.max(np.abs(u - ui))
	ui = np.copy(u)
	print(delta)



img = plt.imshow(u, origin='lower') #RIGHT
plt.colorbar(img) #Right
plt.show() #Right

#---------------------END: JACOBI ITERATION----------------------------------------------------------------

#---------------------START: GAUSS-SEIDEL METHOD----------------------------------------------
g = np.zeros((nx+1, ny+1)) # 2D array for x-and y-values at current time

# calculate u from ui, calculate Laplacians
while(delta > tolerance):
	delta = 0
	for j in range(0, ny+1):
		g[0, j] = 1
		for i in range(1, nx):
			value = np.copy(g[i,j])
			g[i,j] = (g[i-1,j] + g[i,j-1] + g[i+1,j] + g[i,(j+1)%(ny+1)]) / 4
			
			compare = np.abs(g[i,j] - value)
			if (compare > delta):
				delta = compare
	print(delta)


img = plt.imshow(g, origin='lower') #RIGHT
plt.colorbar(img) #Right
plt.show() #Right

#-------------------------------END: GAUSS-SEIDEL METHOD----------------------------------------------

"""

#------------------------------START: SUCCESSIVE OVER-RELAXATION------------------------------------
s = np.zeros((nx+1, ny+1)) # 2D array for x-and y-values at current time


# calculate u from ui, calculate Laplacians
while(delta > tolerance):
	delta = 0
	for j in range(0, ny+1):
		s[0, j] = 1
		for i in range(1, nx):
			value = np.copy(s[i,j])
			s[i,j] = (omega / 4) * (s[i-1,j] + s[i,j-1] + s[i+1,j] + s[i,(j+1)%(ny+1)]) + (1 - omega) * s[i,j]
			
			compare = np.abs(s[i,j] - value)
			if (compare > delta):
				delta = compare


print(countIterations)

#img = plt.imshow(s, origin='lower') #RIGHT
#plt.colorbar(img) #Right
#plt.show() #Right