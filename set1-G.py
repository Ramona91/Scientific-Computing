# Scientific Computing
# Set1, ex. E
# Jacobi iteration, Gauss-seidel iteration, Successive Over-Relaxationanalytic

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import sqrt, sin, pi



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

tolerance = 10**(-5) # do main loop until tolerance level is small enough

#----------------START: JACOBI ITERATION-----------------------------------------
u = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at current time
ui = np.zeros((nx+1, ny+1)) # 2D array for x-andy-values at previous time


for j in range(0, ny+1):
	ui[0, j] = 1

deltau = 1
deltauList = [] # for the graph

# calculate u from ui, calculate Laplacians
while(tolerance < deltau):
	deltau = 0
	for j in range(0, ny+1):
		u[0, j] = 1
		for i in range(1, nx):
			u[i,j] = (ui[i-1,j] + ui[i,j-1] + ui[i+1,j] + ui[i,(j+1)%(ny+1)]) / 4
	deltau = np.max(np.abs(u - ui))

	deltauList.append(deltau)
	ui = np.copy(u)



#img = plt.imshow(u, origin='lower') #RIGHT
#plt.colorbar(img) #Right
#plt.show() #Right

#---------------------END: JACOBI ITERATION----------------------------------------------------------------


#---------------------START: GAUSS-SEIDEL METHOD----------------------------------------------
g = np.zeros((nx+1, ny+1)) # 2D array for x-and y-values at current time

deltag = 1
deltagList = [] # for the graph

# calculate u from ui, calculate Laplacians
while(deltag > tolerance):
	deltag = 0
	for j in range(0, ny+1):
		g[0, j] = 1
		for i in range(1, nx):
			value = np.copy(g[i,j])
			g[i,j] = (g[i-1,j] + g[i,j-1] + g[i+1,j] + g[i,(j+1)%(ny+1)]) / 4
			
			compare = np.abs(g[i,j] - value)
			if (compare > deltag):
				deltag = compare

	deltagList.append(deltag)


#img = plt.imshow(g, origin='lower') #RIGHT
#plt.colorbar(img) #Right
#plt.show() #Right

#-------------------------------END: GAUSS-SEIDEL METHOD----------------------------------------------



#------------------------------START: SUCCESSIVE OVER-RELAXATION------------------------------------
s = np.zeros((nx+1, ny+1)) # 2D array for x-and y-values at current time

deltas = 1
deltasList = []

# calculate u from ui, calculate Laplacians
while(deltas > tolerance):
	deltas = 0
	for j in range(0, ny+1):
		s[0, j] = 1
		for i in range(1, nx):
			value = np.copy(s[i,j])
			s[i,j] = (omega / 4) * (s[i-1,j] + s[i,j-1] + s[i+1,j] + s[i,(j+1)%(ny+1)]) + (1 - omega) * s[i,j]
			
			compare = np.abs(s[i,j] - value)
			if (compare > deltas):
				deltas = compare

	deltasList.append(deltas)


#img = plt.imshow(s, origin='lower') #RIGHT
#plt.colorbar(img) #Right
#plt.show() #Right


plt.hold(True)
jacobi, = plt.plot(deltauList, 'b.')
gauss, = plt.plot(deltagList, 'g.')
sor, = plt.plot(deltasList, 'r.')
plt.legend([jacobi, gauss, sor], ['Jacobi', 'Gauss-Seidel', 'SOR'])
plt.yscale('log', nonposy='clip')
plt.xlabel('Number of iterations')
plt.ylabel('Precision of delta')
plt.title('The convergence of Jacobi, Gauss-Seidel & Successive Over-Relaxation compared.')
plt.show()