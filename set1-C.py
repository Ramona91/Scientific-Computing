# Scientific Computing
# Set1, ex. Computing
# implement a 1-dimensional wave equation

import random
#import simpy
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, sin, pi
import timeit

# http://stackoverflow.com/questions/26393545/python-graphing-the-1d-wave-equation-beginner

# define wave equation parameters
L = 1.0 # length of string
c = 1.0 # speed of wave with which it traverses through space and time
xmin = 0 # spatial starting position of wave
xmax = 3 # spatial end poition of wave
tmin = 0 # temporal starting point
tmax = 3 # temporal end point
nx = 100 # number of space intervals # (xmax - xmin) / deltaX
nt = 100 # number of time intervals # int((tmax - tmin) / deltaT) 
deltaX = L / nx # size of space interval
deltaT = 0.01 # size of time interval 

u = np.zeros((nx,nt)) #solution to WE, create matrix u with nx rows and nt columns, space: row; time: column

#print deltaX

# define initial condition
for i in range(0, nx):
	u[i, 0] = sin(2 * pi * i)
	u[i, 1] = u[i, 0]

for j in range(0, nt):
	u[0, j] = 0 # define left boundary condition
	u[nx-1, j] = 0 # define right boundary condition, do I have to replace nx by L????

#populate the solution matrix u, i.e. run the simultion, exclude initial values (already calculated) and boundary values (also already calculated)
for j in range(1, nt-1):
	for i in range(1, nx-1):
		u[i, j+1] = (c * deltaT / deltaX)**2 * (u[i+1, j] + u[i-1, j] - 2 * u[i,j]) / (deltaX)**2 - u[i, j-1] + 2 * u[i,j]
		print (u[i, j+1])

		
#plot doesn't work
plt.plot(u[i, 0], u[0, j])
plt.show()