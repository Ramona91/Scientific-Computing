# Scientific Computing
# Set1, ex. E
# error analysis: Compare analytic with numerical solution for time-dependent diffusion


import random
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, sin, pi
from scipy.special import erfc

nx = 10 # no of interval in x-direction
ny = 10 # # no of interval in y-direction
deltaX = 1/nx # length of an interval
tmax = 1.0 #max. simulation time
D = 1.0 # diffusion coefficient
deltaX2 = deltaX**2
deltaT = (deltaX2/(4*D)) / 10  # The max of delta t is deltaX2/(4*d), just to make sure we divide by 10
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
			u[i,j] = ui[i,j] + (deltaT *D/deltaX2) * (ui[i+1,j] + ui[i-1,j] + ui[i,(j+1)%(ny+1)] + ui[i,j-1] - 4*ui[i,j]) #for j+1 take j=0 value to take periodic boundary into account
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


#--------------Calculate analytic solution to time-dependent diffusion equation------------------------------------------------------------
#code taken from blackboard
# analytc expression for the concentration c(x, t)
# from the lecture slides
def c(x,t):
	N = 30 # upper limit for the sum
	c = 0
	d = 2*np.sqrt(D*t)
	
	if t == 0: #avoid division by zero
		d = 2*np.sqrt(D*10**(-100))
	
	for i in range (0,N):
		c += erfc((1-x+2*i) / d) - erfc((1+x+2*i) / d)

	return c

# apply c on a vector
def C(X, t):
	C = [c(x,t) for x in X]
	return C

# get y-values of analytic solution
def yvalues(t):
	X = np.linspace (0, 1, 11)
	Y = C(X, t)
	return Y
	

for t in [0, .001, 0.01, 0.1, 1]:
	if t == 0:
		y1 = yvalues(t)
	if t == 0.001:
		y2 = yvalues(t)
	if t == 0.01:
		y3 = yvalues(t)
	if t == 0.1:
		y4 = yvalues(t)
	if t == 1:
		y5 = yvalues(t)

#-------------------------------------------------------------------------------


#[i - j for i, j in zip(a, b)] #Numerical - Analytic
errorTime1 = [i - j for i, j in zip(time1, y1)]
errorTime2 = [i - j for i, j in zip(time2, y2)]
errorTime3 = [i - j for i, j in zip(time3, y3)]
errorTime4 = [i - j for i, j in zip(time4, y4)]
errorTime5 = [i - j for i, j in zip(time5, y5)]

#Squared error for numerical solutions at different time points
squaredError1 = [i ** 2 for i in errorTime1]
squaredError2 = [i ** 2 for i in errorTime2]
squaredError3 = [i ** 2 for i in errorTime3]
squaredError4 = [i ** 2 for i in errorTime4]
squaredError5 = [i ** 2 for i in errorTime5]

#Mean-Squared error for numerical solutions at different time points
mseT1 = (sum(squaredError1)/ len(squaredError1))
mseT2 = (sum(squaredError2)/ len(squaredError2))
mseT3 = (sum(squaredError3)/ len(squaredError3))
mseT4 = (sum(squaredError4)/ len(squaredError4))
mseT5 = (sum(squaredError5)/ len(squaredError5))

print(mseT1, mseT2, mseT3, mseT4, mseT5)

# Plot MSEs
plt.hold(True)
p1 = plt.plot([1],[mseT1], 'ob')
p2 = plt.plot([2],[mseT2], 'ob')
p3 = plt.plot([3],[mseT3], 'ob')
p4 = plt.plot([4],[mseT4], 'ob')
p5 = plt.plot([5],[mseT5], 'ob')
plt.xlim(0, 6)
plt.show()

# For the percentage form 0 to 1 on the x-axis
percAnalytic = np.linspace(0,1,len(time1))





# Make figure for different methods in 1 dimension squared error
plt.hold(True)
first, = plt.plot(percAnalytic, squaredError1 , 'b')
second, = plt.plot(percAnalytic, squaredError2 , 'g')
third, = plt.plot(percAnalytic, squaredError3, 'r')
fourth, = plt.plot(percAnalytic, squaredError4, 'c')
fith, = plt.plot(percAnalytic, squaredError5, 'm')
plt.legend([first, second, third, fourth, fith], ['t=0', 't=0.001', 't=0.01', 't=0.1', 't=1'], loc = 2)
plt.xlabel('y - coordinate of the area')
plt.ylabel('Squared error')
#plt.title('The discrete solution of the concentration, as a function of the y-coordinate, for different methods.')
plt.show()

