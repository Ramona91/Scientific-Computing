# Scientific Computing
# Set1, ex. E
# Jacobi iteration, Gauss-seidel iteration, Successive Over-Relaxationanalytic

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import sqrt, sin, pi



nx = 20 # no of interval in x-direction
ny = 20 # # no of interval in x-direction
deltaX = 1/nx # length of an interval
deltaX2 = deltaX**2 
deltaY = deltaX
tmax = 10.0 #max. simulation time
d = 1.0 # diffusion coefficient
deltaT = (deltaX2 * 4 * d) / 10  # The max of delta t is deltaX2/(4*d), just to make sure we divide by 10
omega = 1 #needed for successive over-relaxation method
length_previous_deltaList = 10000000 # some big number  to calculate the optimal omega

tolerance = 10**(-5) # do main loop until tolerance level is small enough


#------------------------------START: SUCCESSIVE OVER-RELAXATION------------------------------------

differentOmegas = []

while(omega < 2):
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

		deltasList.append(np.copy(deltas))

	
	#if we want to check if we found the optimal omega
	#if(length_previous_deltaList < len(deltasList)):
	#	omega -= 0.01
	#	break
	#length_previous_deltaList = len(deltasList)

	differentOmegas.append(np.copy(deltasList))

	omega += 0.1

print(omega)

# Print the convergence for all the omega's
plt.hold(True)
omega0, = plt.plot(differentOmegas[0], 'b.')
omega1, = plt.plot(differentOmegas[1], 'g.')
omega2, = plt.plot(differentOmegas[2], 'r.')
omega3, = plt.plot(differentOmegas[3], 'c.')
omega4, = plt.plot(differentOmegas[4], 'm.')
omega5, = plt.plot(differentOmegas[5], 'k.')
omega6, = plt.plot(differentOmegas[6], 'y.')
omega7, = plt.plot(differentOmegas[7], 'b*')
omega8, = plt.plot(differentOmegas[8], 'g*')
omega9, = plt.plot(differentOmegas[9], 'r*')

plt.legend([omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7, omega8, omega9], ['Omega = 1', 'Omega = 1.1', 'Omega = 1.2', 'Omega = 1.3', 'Omega = 1.4', 'Omega = 1.5', 'Omega = 1.6', 'Omega = 1.7', 'Omega = 1.8', 'Omega = 1.9'])
plt.yscale('log', nonposy='clip')
plt.xlabel('Number of iterations')
plt.ylabel('Precision of delta')
plt.title('The convergence of Successive Over-Relaxation for different omega"s.')
plt.show()
