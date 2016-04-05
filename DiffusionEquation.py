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

nx = 50 # no of interval in x-direction
ny = 50 # # no of interval in y-direction
deltaX = 1/nx # length of an interval
deltaY = deltaX
nt = 1000 # number of time steps
tmax = 10.0 #max. simulation time
deltaT = 0.0001 # length of time interval
d = 1.0 # diffusion coefficient
deltaX2 = deltaX**2
deltaY2 = deltaY**2
deltaTMax = 4 * deltaT * d / deltaX2

u = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at current time
ui = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at previous time


# calculate u from ui, calculate Laplacians
for k in range(0, nt):
	for j in range(0, ny+1):
		ui[0, j] = 1
		for i in range(1, nx):
			u[i,j] = ui[i,j] + (deltaT *d/deltaX2) * (ui[i+1,j] + ui[i-1,j] + ui[i,(j+1)%(ny+1)] + ui[i,j-1] - 4*ui[i,j]) #for j+1 take j=0 value to take periodic boundary into account
	ui = u
	print(ui)


#for n in range(1, nt):
#	u[1:-1,1:-1] = ui[1:-1, 1:-1] + deltaT * d *  ((ui[2:,1:-1] - 2 * ui[1:-1, 1:-1] + ui[:-2,1:-1])/(deltaX2) + (ui[1:-1,:-2] - 2 * ui[1:-1,1:-1] + ui[1:-1,2:])/(deltaY2))
#	ui = sp.copy(u)
	#print(u)


img = plt.imshow(ui, origin='upper') #RIGHT
plt.colorbar(img) #Right
plt.show() #Right

