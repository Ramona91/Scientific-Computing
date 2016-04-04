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

nx = 10 # no of interval in x-direction
ny = 10 # # no of interval in x-direction
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



'''
# define intital conditions
for i in range(0, nx+1):
	for j in range(0, ny+1):
		ui[i,j] = 0
'''
	
# define boundary conditions
for i in range(0, nx+1):
	for j in range(0, ny+1):
		ui[nx, j] = 0 # bottom horizontal (x, y=0)
		ui[0, j] = 1 # top horizontal (x, y=1)
		#ui[i, 0] = ui[i, ny] # periodic boundary conditions along x-axis



# calculate u from ui, calculate Laplacians
for k in range(1, nt):
	for i in range(1, nx):
		for j in range(0, ny+1):
			u[i,j] = ui[i,j] + (deltaT *d/deltaX2) * (ui[i+1,j] + ui[i-1,j] + ui[i,(j+1)%(ny+1)] + ui[i,j-1] - 4*ui[i,j]) #for j+1 take j=0 value to take periodic boundary into account
	ui = u
	#print(u)
	
img = plt.imshow(ui, origin='lower') #RIGHT
#img.set_cmap('autumn') #optional

plt.colorbar(img) #Right
plt.show() #Right

