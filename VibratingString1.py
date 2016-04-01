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

#check the whole program by setting nx=3, nt=4 (or other smaller values), print u after each for-loop

# define wave equation parameters
L = 1.0 # length of string
c = 1.0 # speed of wave with which it traverses through space and time
nx = 100 # number of space intervals in x and y direction
nt = 1000 # number of time intervals
deltaX = L / nx # size of space interval
deltaT = 0.00001 # size of time interval 


u = np.zeros((nx + 1, nt + 1)) # solution to Wave Equation, create matrix u with nx+1 rows and nt+1 columns, space: row; time: column


# define initial condition
for i in range(0, nx + 1):  #for nx+1 rows
	u[i, 0] = sin(2 * pi * i * deltaX)
	u[i, 1] = u[i, 0] #  1 need 3 time values: previous, current and future (calculated below)
	
# Define boundary conditons
for j in range(0, nt + 1):
	u[0, j] = 0 # upper boundary condition
	u[nx, j] = 0 # lower boundary condition, do I have to replace nx by L????


# populate the solution matrix u, i.e. run the simultion, exclude initial values (already calculated) and boundary values 
# (also already calculated)
for j in range(1, nt):  # We wanna iterate through to the end of the array, -> nt, but since we calculate j+1, we need to subtract 1 on the iteration scheme
	for i in range(1, nx):
		u[i, j + 1] = ((c * deltaT) / deltaX)**2.0 * ((u[i+1, j] + u[i-1, j] - 2.0 * u[i,j]) / ((deltaX)**2.0)) - u[i, j-1] + 2.0 * u[i,j]

		#print (u[i,j])
		
#print (u)


# extract data from different columns in order to plot wave at different time points
t0 = u[:,[0]]		# t = 0
t25 = u[:,[24]]	 	# t = 25
t50 = u[:,[49]]		# t = 50
t75 = u[:,[74]]		# t = 75
t100 = u[:,[99]]	# t = 100
t250 = u[:,[249]]
t400 = u[:,[399]]
t500 = u[:,[499]]
t800 = u[:,[799]]
#print (t75)



 #plot the wave amplitudes over spatiatl locations (x) for different time points
plt.plot(t0, label='t = 0')
plt.plot(t25, label='t = 25')
plt.plot(t50, label='t = 50')
plt.plot(t75, label='t = 75')
plt.plot(t100, label='t = 100')
plt.plot(t250, label='t = 250')
plt.plot(t400, label='t = 400')
plt.plot(t500, label='t = 500')
plt.plot(t800, label='t = 800')
plt.legend(loc='lower left')

plt.xlabel('x')
plt.ylabel('Amplitude')

plt.tight_layout() #get all plots into the same one

plt.show()