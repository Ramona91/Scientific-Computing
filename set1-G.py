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


nx = 10 # no of interval in x-direction
ny = 10 # # no of interval in x-direction
deltaX = 1/nx # length of an interval
deltaY = deltaX
nt = 1000 # number of time steps
tmax = 10.0 #max. simulation time
deltaT = tmax/nt # length of time interval
d = 1.0 # diffusion coefficient
deltaX2 = deltaX**2 
deltaTMax = 4 * deltaT * d / deltaX2
omega = 1 #needed for successive over-relaxation method


#----------------START: JACOBI ITERATION-----------------------------------------
u = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at current time
ui = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at previous time
timeListJ = []


# calculate u from ui, calculate Laplacians
for k in range(0, nt):
	for j in range(0, ny+1):
		ui[0, j] = 1
		for i in range(1, nx):
			u[i,j] = (ui[i-1,j] + ui[i,j-1] + ui[i+1,j] + ui[i,(j+1)%(ny+1)]) / 4
	ui = u
	timListJ.append(ui)
	print(ui)


img = plt.imshow(ui, origin='upper') #RIGHT
plt.colorbar(img) #Right
plt.show() #Right

#---------------------END: JACOBI ITERATION----------------------------------------------------------------

#---------------------START: GAUSS-SEIDEL METHOD----------------------------------------------
g = np.zeros((nx+1, ny+1)) # 2D array for x-and y-values at current time
gi = np.zeros((nx+1, ny+1)) # 2D array for x-and y-values at previous time
timeListG = []

# calculate u from ui, calculate Laplacians
for k in range(0, nt):
	for j in range(0, ny+1):
		gi[0, j] = 1
		for i in range(1, nx):
			g[i,j] = (g[i-1,j] + g[i,j-1] + gi[i+1,j] + gi[i,(j+1)%(ny+1)]) / 4
	gi = g
	timListG.append(gi)
	print(ui)


img = plt.imshow(gi, origin='upper') #RIGHT
plt.colorbar(img) #Right
plt.show() #Right

#-------------------------------END: GAUSS-SEIDEL METHOD----------------------------------------------


#------------------------------START: SUCCESSIVE OVER-RELAXATION------------------------------------
s = np.zeros((nx+1, ny+1)) # 2D array for x-and y-values at current time
si = np.zeros((nx+1, ny+1)) # 2D array for x-and y-values at previous time
timeListS = []

# calculate u from ui, calculate Laplacians
for k in range(0, nt):
	for j in range(0, ny+1):
		si[0, j] = 1
		for i in range(1, nx):
			s[i,j] = (s[i-1,j] + s[i,j-1] + si[i+1,j] + si[i,(j+1)%(ny+1)]) / 4
			s[i,j] = (1 - omega) * si[i,j] + omega * s[i,j]
	si = s
	timListG.append(gi)
	print(ui)


img = plt.imshow(si, origin='upper') #RIGHT
plt.colorbar(img) #Right
plt.show() #Right






