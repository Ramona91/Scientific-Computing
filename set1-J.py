# Scientific Computing
# Set1, ex. J
# implement a 2-dimensional time-independent diffusion equation with object

import random
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, sin, pi

nx = 50 # no of interval in x-direction
ny = 50 # # no of interval in x-direction
deltaX = 1/nx # length of an interval
deltaY = deltaX
tmax = 10.0 #max. simulation time
d = 1.0 # diffusion coefficient
deltaX2 = deltaX**2
deltaT = (4*d*deltaX2)/ 10 # length of time interval
steps = int(tmax / deltaT) # number of iteration steps
delta = 0.1
tolerance = 10**(-5)
upperLeftX = 20
upperRightX = upperLeftX
upperLeftY = 20
lowerLeftY = upperLeftY
lowerLeftX = 30
lowerRightX = lowerLeftX
upperRightY = 30
lowerRightY = upperRightY
omega = 1.5 #needed for successive over-relaxation method
#lastColumnJ = []



#------------------------------START: SUCCESSIVE OVER-RELAXATION------------------------------------
s = np.zeros((nx+1, ny+1)) # 2D array for x-and y-values at current time
partObjectY = []
partObjectX = []

for j in range(upperLeftY, upperRightY+1):
	for i in range(upperLeftX, lowerLeftX+1):
		partObjectY.append(j)
		partObjectX.append(i)

# calculate u from ui, calculate Laplacians
while (delta > tolerance):
	delta = 0
	for j in range(0, ny+1):
		s[0, j] = 1
		for i in range(1, nx):
			oldS = np.copy(s[i,j])
			
			if (j in partObjectY) and (i in partObjectX):
				s[i,j] = 0
			else:
				s[i,j] = (s[i-1,j] + s[i,j-1] + s[i+1,j] + s[i,(j+1)%(ny+1)]) / 4
				s[i,j] = (1 - omega) * s[i,j] + omega * s[i,j]
				deltaNew = np.abs(s[i,j] - oldS)
	
			if deltaNew > delta:
				delta = deltaNew
#------------------------------END: SUCCESSIVE OVER-RELAXATION------------------------------------

img = plt.imshow(s, origin='lower') #RIGHT
plt.colorbar(img) #Right
plt.show() #Right
