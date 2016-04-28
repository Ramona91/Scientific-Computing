# Scientific Computing
# Set2, ex. N
# Gray-Scott model



import numpy as np
from numpy import random
import numpy.ma as ma
import matplotlib.pyplot as plt
from math import sqrt, sin, pi
from scipy.special import erfc

''' 
nx = 10 # no of interval in x-direction
ny = 10 # # no of interval in y-direction
deltaX = 1/nx # length of an interval
tmax = 1.0 #max. simulation time
Du = 0.16 # diffusion coefficient for u
Dv = 0.08 # # diffusion coefficient for v
deltaX2 = deltaX**2
deltaT = (deltaX2/(4*D)) / 10  # The max of delta t is deltaX2/(4*d), just to make sure we divide by 10
timeList = [] # Make an empty list for the plots later
nt = int(tmax / deltaT) + 1 # calculate the number of time steps
''' 

nx = 100 # no of interval in x-direction
ny = 100 # # no of interval in y-direction
deltaT = 1
deltaX = 1 # length of an interval

Du = 0.16 # diffusion coefficient for u
Dv = 0.08 # # diffusion coefficient for v
f = 0.035
k = 0.06

timeListu = [] # Make an empty list for the plots later u
timeListv = [] # Make an empty list for the plots later v
tmax = 1000.0 #max. simulation time
nt = int(tmax / deltaT) + 1 # calculate the number of time steps


# Start with an nx + 1 by ny + 1 matrix with zero's
u = np.zeros((nx+1, ny+1)) # 2D array for x-and-y-values at current time
ui = np.zeros((nx+1, ny+1)) # 2D array for x-and-y-values at previous time
v = np.zeros((nx+1, ny+1)) # 2D array for x-and-y-values at current time
vi = np.zeros((nx+1, ny+1)) # 2D array for x-and-y-values at previous time

neighbours = ()

# squre in centre for intial condition
uLx = int(((nx+1)/2)-2) # upper left x
uLy = int(((ny+1)/2)-2) # upper left y
uRx = int(((nx+1)/2)-2)
uRy = int(((ny+1)/2)+2)
lLx = int(((nx+1)/2)+2)
lLy = int(((ny+1)/2)-2)
lRx = int(((nx+1)/2)+2)
lRy = int(((ny+1)/2)+2)

# We want the figures at t is exactly 0.001, 0.01, 0.1 and 1, calculating here the related k value
a1 = int (1 / deltaT)
a2 = int (5 / deltaT)
a3 = int (10 / deltaT)
a4 = int (15 / deltaT)
a5 = int (20 / deltaT)
a6 = int (25 / deltaT)
a7 = int (50 / deltaT)
a8 = int (75 / deltaT)
a9 = int (100 / deltaT)
a10 = int (250 / deltaT)
a11 = int (500 / deltaT)
a12 = int (750 / deltaT)
a13 = int (1000 / deltaT)
# a9 = int (10000 / deltaT)
# a6 = int (100000 / deltaT)

times = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13]


# adding tuples: http://stackoverflow.com/questions/1380860/add-variables-to-tuple
# hierarchical list of nearest neighbours (in which order should nearest neighbours be accessed 
# in order to provide replacement values of NaNs)
for i in range(1, nx+1): 
	newNeighbour = ((0,i), (0,-i), (i,0), (-i,0), (i,i), (-i,i), (i,-i), (-i,-i))
	neighbours = neighbours + newNeighbour

for j in range(uLy, uRy+1): # intial conditions
	for i in range(uLx, lLx+1):
		ui[i,j] = 0.5
		vi[i,j] = 0.25



# calculate u from ui, calculate Laplacians
for t in range(0, nt + 1): # keeping track of time, range(0, nt + 1)
	for j in range(0, ny+1): 
		for i in range(1, nx): # boundary conditions in y-direction: i=0, i=ny+1 --> always 0		
			u[i,j] = (ui[i,j] + (deltaT *Du/deltaX**2) * (ui[i+1,j] + ui[i-1,j] + ui[i,(j+1)%(ny+1)] + ui[i,j-1] - 4*ui[i,j]) 
			- ui[i,j]*(vi[i,j])**2 + f*(1-ui[i,j])) + round(random.uniform()) #for j+1 take j=0 value to take periodic boundary into account
			v[i,j] = (vi[i,j] + (deltaT *Dv/deltaX**2) * (vi[i+1,j] + vi[i-1,j] + vi[i,(j+1)%(ny+1)] + vi[i,j-1] - 4*vi[i,j]) 
			+ ui[i,j]*(vi[i,j]**2) - (f+k)*vi[i,j]) + round(random.uniform())

	u_masked = ma.masked_invalid(u)
	v_masked = ma.masked_invalid(v)


# http://stackoverflow.com/questions/3662361/fill-in-missing-values-with-nearest-neighbour-in-python-numpy-masked-arrays
# replace NaNs in masked array with nearest neighbour value
	for hor_shift, vert_shift in neighbours: 
		if not np.any(u_masked.mask) or not np.any(v_masked.mask):
			break
		u_shifted = np.roll(u_masked, shift = hor_shift, axis=1)
		u_shifted = np.roll(u_shifted, shift = vert_shift, axis=0)
		idx =~ u_shifted.mask * u_masked.mask
		u_masked[idx] = u_shifted[idx]
		v_shifted = np.roll(v_masked, shift = hor_shift, axis=1)
		v_shifted = np.roll(v_shifted, shift = vert_shift, axis=0)
		idx =~ v_shifted.mask * v_masked.mask
		v_masked[idx] = v_shifted[idx]


	#print(v_masked)

	ui = np.copy(u_masked)
	vi = np.copy(v_masked)

	
	# Add the values for the plots in the list
	if (t == 0) or (t == a1) or (t == a2) or (t == a3) or (t == a4) or (t == a5) or (t == a6) or (t == a7) or (t == a8) or (t == a9) or (t == a10) or (t == a11) or (t == a12) or (t == a13):
		timeListu.append(np.copy(u_masked))
		timeListv.append(np.copy(v_masked))

	if t%100 == 0:
		print(t)


# make plots
for (u,t) in zip(timeListu, times):
	img = plt.imshow(u, origin='lower')
	plt.colorbar(img) 
	plt.xlabel('x')
	plt.ylabel('y')
	plt.savefig('GrayScott_u_atTime' + str(t) + '.png')  # now it's printing 4 figure!!! 
	plt.show()

for (v,t) in zip(timeListv, times):
	img = plt.imshow(v, origin='lower')
	plt.colorbar(img) 
	plt.xlabel('x')
	plt.ylabel('y')
	plt.savefig('GrayScott_v_atTime' + str(t) + '.png')  # now it's printing 4 figure!!! 
	plt.show()
