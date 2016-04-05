# Scientific Computing
# Set1, ex. D
# implement a 2-dimensional diffusion equation  

# http://www.timteatro.net/2010/10/29/performance-python-solving-the-2d-diffusion-equation-with-numpy/
# plot colormaps: http://stackoverflow.com/questions/7229971/2d-grid-data-visualization-in-python
# color map documentation: http://matplotlib.org/examples/api/colorbar_only.html

import random
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import sqrt, sin, pi

nx = 30 # no of interval in x-direction
ny = 30 # # no of interval in y-direction
deltaX = 1/nx # length of an interval
deltaY = deltaX
tmax = 1.0 #max. simulation time
d = 1.0 # diffusion coefficient
deltaX2 = deltaX**2
deltaY2 = deltaY**2
deltaT = (deltaX2/(4*d)) / 10
#deltaTMax = (4 * deltaT * d) / deltaX2
timeList = []
nt = int(tmax / deltaT)

print(deltaT, nt)
u = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at current time
ui = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at previous time


# define length of time interval
#if deltaTMax <= 1:
#	deltaT = tmax/nt
#else:
#	deltaT = deltaX2/(4*d)


a1 = int (0.001 / deltaT)
a2 = int (0.01 / deltaT)
a3 = int (0.1 / deltaT)
a4 = int (1 / deltaT)

print(a1, a2, a3, a4)

# calculate u from ui, calculate Laplacians
for k in range(0, nt + 1):
	for j in range(0, ny+1):
		ui[0, j] = 1
		for i in range(1, nx):
			u[i,j] = ui[i,j] + (deltaT *d/deltaX2) * (ui[i+1,j] + ui[i-1,j] + ui[i,(j+1)%(ny+1)] + ui[i,j-1] - 4*ui[i,j]) #for j+1 take j=0 value to take periodic boundary into account
	ui = u
	#if(k == 0) or (k == a1) or (k == a2) or (k == a3) or (k == a4):
	timeList.append(ui)
	#print(ui)



#for n in range(1, nt):
#	u[1:-1,1:-1] = ui[1:-1, 1:-1] + deltaT * d *  ((ui[2:,1:-1] - 2 * ui[1:-1, 1:-1] + ui[:-2,1:-1])/(deltaX2) + (ui[1:-1,:-2] - 2 * ui[1:-1,1:-1] + ui[1:-1,2:])/(deltaY2))
#	ui = sp.copy(u)
	#print(u)

print(timeList[0],timeList[50])

fig = plt.figure()
plt.subplot(231)
plt.imshow(timeList[0], origin='lower') #
plt.subplot(232)
plt.imshow(timeList[1], origin='lower')
plt.subplot(233)
plt.imshow(timeList[2], origin='lower') 
plt.subplot(234)
plt.imshow(timeList[3], origin='lower')
plt.subplot(235)
plt.imshow(timeList[4], origin='lower')
plt.colorbar() 
plt.show() 

