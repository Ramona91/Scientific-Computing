# Scientific Computing
# Set2, ex. N
# Gray-Scott model



import numpy as np
from numpy import random
import matplotlib.pyplot as plt


nx = 200 # no of interval in x-direction
ny = 200 # # no of interval in y-direction
deltaT = 1
deltaX = 1 # length of an interval

Du = 0.16 # diffusion coefficient for u
Dv = 0.08 # # diffusion coefficient for v
f = 0.035
k = 0.06

timeListu = [] # Make an empty list for the plots later u
timeListv = [] # Make an empty list for the plots later v
tmax = 10000 #max. simulation time
nt = int((tmax / deltaT)+1)# calculate the number of time steps


# Start with an nx + 1 by ny + 1 matrix with zero's
u = np.zeros((nx+1, ny+1)) # 2D array for x-and-y-values at current time
ui = np.zeros((nx+1, ny+1)) # 2D array for x-and-y-values at previous time
v = np.zeros((nx+1, ny+1)) # 2D array for x-and-y-values at current time
vi = np.zeros((nx+1, ny+1)) # 2D array for x-and-y-values at previous time




a1 = int (1 / deltaT)
a10 = int (10 / deltaT)
a50 = int (50 / deltaT)
a100 = int (100 / deltaT)
a500 = int (500 / deltaT)
a1000 = int (1000 / deltaT)
a1500 = int (1500 / deltaT)
a2000 = int (2000 / deltaT)
a2500 = int (2500 / deltaT)
a3000 = int (3000 / deltaT)
a3500 = int (3500 / deltaT)
a4000 = int (4000 / deltaT)
a4500 = int (4500 / deltaT)
a5000 = int (5000 / deltaT)
a5500 = int (5500 / deltaT)
a6000 = int (6000 / deltaT)
a6500 = int (6500 / deltaT)
a7000 = int (7000 / deltaT)
a7500 = int (7500 / deltaT)
a8000 = int (8000 / deltaT)
a8500 = int (8500 / deltaT)
a9000 = int (9000 / deltaT)
a9500 = int (9500 / deltaT)
a10000 = int (10000 / deltaT)

times = [a1, a10, a50, a100, a500, a1000, a1500, a2000, a2500, a3000, a3500, a4000, a4500, a5000, a5500, a6000, a6500, a7000, a7500, a8000, a8500, a9000, a9500, a10000]#a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13]


r = (nx+1)/10

ui += 0.50
vi[(nx+1)/2-r:(nx+1)/2+r,(nx+1)/2-r:(nx+1)/2+r] = 0.25
ui += 0.05*random.random(((nx+1),(nx+1)))
vi += 0.05*random.random(((nx+1),(nx+1)))



# calculate u from ui, calculate Laplacians
for t in range(0, nt+1): # keeping track of time, range(0, nt + 1)
	for j in range(0, ny+1): 
		for i in range(1, nx): # boundary conditions in y-direction: i=0, i=ny+1 --> always 0		
			u[i,j] = deltaT * ((Du/deltaX**2) * (ui[i+1,j] + ui[i-1,j] + ui[i,(j+1)%(ny+1)] + ui[i,j-1] - 4*ui[i,j])
				- ui[i,j] * (v[i,j])**2 + f*(1-ui[i,j])) + ui[i,j]
			v[i,j] = deltaT * ((Dv/deltaX**2) * (vi[i+1,j] + vi[i-1,j] + vi[i,(j+1)%(ny+1)] + vi[i,j-1] - 4*vi[i,j])
				+ ui[i,j] * (v[i,j])**2 - (f+k) * vi[i,j]) + vi[i,j]

	
	ui = np.copy(u)
	vi = np.copy(v)

	
	# Add the values for the plots in the line
	if (t == 0) or (t == a1) or (t == a10) or (t == a50) or (t == a100) or (t == a500) or (t == a1000) or (t == a1500) or (t == a2000) or (t == a2500) or (t == a3000) or (t == a3500) or (t == a4000) or (t == a4500) or (t == a5000) or (t == a5500) or (t == a6000) or (t == a6500) or (t == a7000) or (t == a7500) or (t == a8000) or (t == a8500) or (t == a9000) or (t == a9500) or (t == a10000):
		timeListu.append(np.copy(u))
		timeListv.append(np.copy(v))

	if t%1000 == 0:
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
