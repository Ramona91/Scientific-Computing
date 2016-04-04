# Scientific Computing
# Set1, ex. D
# implement a 2-dimensional diffusion equation

# http://www.timteatro.net/2010/10/29/performance-python-solving-the-2d-diffusion-equation-with-numpy/
# plot colormaps: http://stackoverflow.com/questions/7229971/2d-grid-data-visualization-in-python
# color map documentation: http://matplotlib.org/examples/api/colorbar_only.html

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

nx = 100 # no of interval in x-direction
ny = 100 # # no of interval in x-direction
deltaX = 1/nx # length of an interval
deltaY = deltaX
nt = 10 # number of time steps
tmax = 10.0 #max. simulation time
deltaT = tmax/nt # length of time interval
d = 1.0 # diffusion coefficient
deltaX2 = deltaX**2 
deltaTMax = 4 * deltaT * d / deltaX2

u = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at current time
ui = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at previous time




# define intital conditions
for i in range(0, nx+1):
	for j in range(0, ny+1):
		ui[i,j] = 0
	
	
# define boundary conditions
for i in range(0, nx+1):
	for j in range(0, ny+1):
		ui[nx, j] = 0 # bottom horizontal (x, y=0)
		ui[0, j] = 1 # top horizontal (x, y=1)
		ui[i, 0] = ui[i, ny] # periodic boundary conditions along x-axis
		
		
# define initial conditions
# since ui was already initialized with zeros and intial conditions are zero, we don't have to fo anything


# calculate u from ui, calculate Laplacians
for k in range(1, nt):
	for i in range(1, nx):
		for j in range(1, ny):
			u[i,j] = ui[i,j] + (deltaT *d/deltaX2) * (ui[i+1,j] + ui[i-1,j] + ui[i,j+1] + ui[i,j-1] - 4*ui[i,j])
			ui[i,j] = u[i,j]
	#print(u)

	
#ListedColormap(colors, name='from_list', N=None) http://matplotlib.org/api/colors_api.html#matplotlib.colors.LinearSegmentedColormap.from_list

#----------------------VERSION 1: DOESN'T WORK YET -------------------------------------------------------------------------------------------------------
# make a color map of continuous colors
'''
cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', ['blue','black','red'],256)

# tell imshow about color map so that only set colors are used
img2 = plt.imshow(ui,interpolation='nearest', cmap = cmap2, origin='lower')

# make a color bar
plt.colorbar(img2,cmap=cmap2)

plt.show()
'''
#----------------------------------------------------------------------------------------------------------------------------------------------------------------


#----------------------version 2: works-----------------------------------------------------------
# make a color map of fixed colors
cmap = mpl.colors.ListedColormap(['blue','black','red'])
bounds=[-6,-2,2,6] # numbers have to increase monotonically, amount of numbers has to be one larger than the number of colours
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# tell imshow about color map so that only set colors are used
img = plt.imshow(ui,interpolation='nearest', cmap = cmap,norm=norm, origin='lower')

# make a color bar
plt.colorbar(img,cmap=cmap, norm=norm,boundaries=bounds,ticks=[-5,0,5])


plt.show()
#------------------------------------------------------------------------------------------------------------------





	#print (u)
	# u[i,j] = ui[i,j] + deltaT *d/deltaX2 * (ui[i+1,j] + ui[i-1,j], ui[i,j+1], ui[i,j-1] - 4*ui[i,j])
	
'''
			uxx = (ui[i+1, j] + ui[i-1, j] - 2*ui[i,j]) / (deltaX)**2 # second derivative for x-values
			uyy = (ui[i, j+1] + ui[i, j-1] - 2*ui[i,j]) / (deltaY)**2 # second derivative for y-values
			
			u[i,j] = deltaT * d *(uxx + uyy)
'''
