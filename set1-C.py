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
nx = 100 # number of space intervals # (xmax - xmin) / deltaX
nt = 100 # number of time intervals # int((tmax - tmin) / deltaT) 
deltaX = L / nx # size of space interval
deltaT = 0.01 # size of time interval 

u = np.zeros((nx,nt)) #solution to WE, create matrix u with nx rows and nt columns, space: row; time: column

#print (u)



# define initial condition
for i in range(0, nx):  #for nx+1 rows
	u[i, 0] = sin(2 * pi * i)
	u[i, 1] = u[i, 0] #why do I have to define stuff for t=1??? (because 1 need 3 time values: previous, current and future (calculted below)

#print (u)
	

for j in range(0, nt):
	u[0, j] = 0 # define left boundary condition
	u[nx-1, j] = 0 # define right boundary condition, do I have to replace nx by L????

#print(u)


#populate the solution matrix u, i.e. run the simultion, exclude initial values (already calculated) and boundary values (also already calculated)
for j in range(1, nt-1):  #actually we wanna iterate through to the end of the arrau, i.e. nt, but since we calculate j+1, we need to subtract 1 on the iteration scheme
	for i in range(1, nx-1):
		u[i, j+1] = (c * deltaT / deltaX)**2 * (u[i+1, j] + u[i-1, j] - 2 * u[i,j]) / (deltaX)**2 - u[i, j-1] + 2 * u[i,j]
		#print (u[i, j+1])
		
#print (u)


# extract data from different columns in order to plot wave at different time points
t0 = u[:,[0]]# t = 0
t25 = u[:,[19]] # t = 25
t50 = u[:,[49]]# t = 50
t75 = u[:,[74]]# t = 75 #doesn't work
t100 = u[:,[99]]# t = 100, # doesn't work
#print (t50)




# plot the wave amplitudes over spatiatl locations (x) for different time points
plt.plot(t0, label='t = 0')
plt.plot(t25, label='t = 25')
#plt.plot(t50, label='t = 50')
plt.legend(loc='upper right')

plt.xlabel('x')
plt.ylabel('Amplitude')

plt.tight_layout() #get all plots into the same one

plt.show()
