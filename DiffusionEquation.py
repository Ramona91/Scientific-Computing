# Scientific Computing
# Set1, ex. D
# implement a 2-dimensional diffusion equation

# http://www.timteatro.net/2010/10/29/performance-python-solving-the-2d-diffusion-equation-with-numpy/
# plot colormaps: http://stackoverflow.com/questions/7229971/2d-grid-data-visualization-in-python
# color map documentation: http://matplotlib.org/examples/api/colorbar_only.html

import random
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, sin, pi

nx = 10 # no of interval in x-direction
ny = 10 # # no of interval in y-direction
deltaX = 1/nx # length of an interval
tmax = 1.0 #max. simulation time
d = 1.0 # diffusion coefficient
deltaX2 = deltaX**2
deltaT = (deltaX2/(4*d)) / 10  # The max of delta t is deltaX2/(4*d), just to make sure we divide by 10
timeList = [] # Make an empty list for the plots later
nt = int(tmax / deltaT) + 1 # calculate the number of time steps

# Start with an nx + 1 by ny + 1 matrix with zero's
u = np.zeros((nx+1, ny+1)) # 2D array for x-and-y-values at current time
ui = np.zeros((nx+1, ny+1)) # 2D array for x-and-y-values at previous time


# We want the figures at t is exactly 0.001, 0.01, 0.1 and 1, calculating here the related k value
a1 = int (0.001 / deltaT)
a2 = int (0.01 / deltaT)
a3 = int (0.1 / deltaT)
a4 = int (1 / deltaT)

for j in range(0, ny+1):
	ui[0, j] = 1

# calculate u from ui, calculate Laplacians
for k in range(0, nt + 1):
	for j in range(0, ny+1):
		u[0, j] = 1
		for i in range(1, nx):
			u[i,j] = ui[i,j] + (deltaT *d/deltaX2) * (ui[i+1,j] + ui[i-1,j] + ui[i,(j+1)%(ny+1)] + ui[i,j-1] - 4*ui[i,j]) #for j+1 take j=0 value to take periodic boundary into account
	ui = np.copy(u)
	# Add the values for the plots in the list
	if(k == 0) or (k == a1) or (k == a2) or (k == a3) or (k == a4):
		timeList.append(np.copy(u))


# Prepare a single column of each saved time step for the figure where the analytic solution is compared
time1 = [row[3] for row in timeList[0]]
time2 = [row[3] for row in timeList[1]]
time3 = [row[3] for row in timeList[2]]
time4 = [row[3] for row in timeList[3]]
time5 = [row[3] for row in timeList[4]]

# List have to be flipped for a nice graph later
time1 = list(reversed(time1))
time2 = list(reversed(time2))
time3 = list(reversed(time3))
time4 = list(reversed(time4))
time5 = list(reversed(time5))

# For the percentage form 0 to 1 on the x-axis
perc = np.linspace(0,1,len(time1))

# Make figure at different time steps in 1 dimension
plt.hold(True)
first, = plt.plot(perc, time1, 'b')
second, = plt.plot(perc, time2, 'g')
third, = plt.plot(perc, time3, 'r')
fourth, = plt.plot(perc, time4, 'c')
fifth, = plt.plot(perc, time5, 'm')
plt.legend([first, second, third, fourth, fifth], ['t = 0', 't = 0.001', 't = 0.01', 't = 0.1', 't = 1'], loc = 2)
plt.xlabel('y - coordinate of the area')
plt.ylabel('Concentration over the y-coordinate')
plt.title('The discrete solution of the concentration, as a function of the y-coordinate, for different t-values.')
plt.show()


# Make colour figure of the entire window (2D) with colorbar
fig = plt.figure()
plt.subplot(231)
plt.title('Diffusion at t = 0')
plt.imshow(timeList[0], origin='lower') #
plt.subplot(232)
plt.title('Diffusion at t = 0.001')
plt.imshow(timeList[1], origin='lower')
plt.subplot(233)
plt.title('Diffusion at t = 0.01')
plt.imshow(timeList[2], origin='lower') 
plt.subplot(234)
plt.title('Diffusion at t = 0.1')
plt.imshow(timeList[3], origin='lower')
plt.subplot(235)
plt.title('Diffusion at t = 1')
plt.imshow(timeList[4], origin='lower')
plt.colorbar() 
plt.show() 

