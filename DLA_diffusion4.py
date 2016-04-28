# Scientific Computing, assignment 2
# Set1, ex. J
# Object in diffusion space, implemented with SOR

import numpy as np
from numpy import random
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import sqrt, sin, pi
#from scipy.special import erfc



nx = 200 # no of interval in x-direction
ny = 200 # # no of interval in x-direction
d = 1.0 # diffusion coefficient
omega = 1.6 #needed for successive over-relaxation method
eta = 0 # determines shape of object
tolerance = 10**(-5) # do main loop until tolerance level is small enough



# results for different etas (cluster characteristics)
concentrations = []
widths = []
heights = []
surfaces = []
etas = []







while eta <= 3:
	etas.append(eta)

	object = np.zeros((nx+1, ny+1)) # presence of object encoded as 1 in object array, object in centre cell of matrix
	for j in [(ny+1)*0.5]: # at the centre column
		for i in [1]: # at the centre row, (nx+1)*0.5
			object[i,j] = 1

	s = np.zeros((nx+1, ny+1)) # concentrations array

	xIndex = [] # keep track of indeces where objects are located (important to calculate cluster properties)
	yIndex = []




	for k in range(0, 5): # number of growth steps
		concentrationCandidates = []
		#candidates = []
		concentrationSum = 0
		delta = 1
		#z = 0
		#pgs = []


		# access all growth candidates
		# save their corresponding concentrations in list
		for j in range(0, ny+1):
			for i in range(1, nx):
				if (object[i,j] != 1) and ((object[i-1,j] == 1) or (object[i+1, j] == 1) or (object[i, j-1] == 1) or (object[i, (j+1)%(ny+1)] == 1)):
					if s[i,j] == 0:
						c = (10**(-10))**eta #c = erfc((1-x+2*i) / div) - erfc((1+x+2*i) / div) # what is i?
					else:
						c = (s[i,j])**eta # growth candidate to the power of eta (part of probability formula)
					concentrationCandidates.append(c)	
		concentrationSum = sum(concentrationCandidates)



		
		# select growth candidates based on probability	
		for j in range(0, ny+1):
			for i in range(1, nx):
				if ((object[i,j] != 1) and ((object[i-1,j] == 1) or (object[i+1, j] == 1) or (object[i, j-1] == 1) or (object[i, (j+1)%(ny+1)] == 1))):
					if s[i,j] == 0:
						c = (10**(-10))**eta #c = erfc((1-x+2*i) / div) - erfc((1+x+2*i) / div) # what is i?
					else:
						c = (s[i,j])**eta
					if (c/concentrationSum) > random.uniform():
						object[i, j] = 1
						s[i,j] = 0



		# calculate u from ui, calculate Laplacians
		while(delta > tolerance):
			delta = 0
			for j in range(0, ny+1):
				s[0, j] = 1
				for i in range(1, nx):
					value = np.copy(s[i,j])
					
					if object[i,j] == 1:
						s[i,j] = 0
					else:
						s[i,j] = (omega / 4) * (s[i-1,j] + s[i,j-1] + s[i+1,j] + s[i,(j+1)%(ny+1)]) + (1 - omega) * s[i,j]
						compare = np.abs(s[i,j] - value)

					if (compare > delta):
						delta = compare




	# calculate width, height and surface of cluster
	for j in range(0, ny+1):
		for i in range(0, nx+1):
			if object[i,j] == 1:
				xIndex.append(i)
				yIndex.append(j)

	#print(xIndex)
	width = max(yIndex) - min(yIndex)
	height = max(xIndex) - min(xIndex)
	surface = len(xIndex)

	widths.append(width)
	heights.append(height)
	surfaces.append(surface)
	concentrations.append(s)

	print(eta)
	eta += 0.5







for (r,e) in zip(concentrations, etas):
	img = plt.imshow(r, origin='lower')
	plt.colorbar(img) 
	plt.xlabel('x')
	plt.ylabel('y')
	plt.savefig('cluster_diffusion_eta' + str(e) + '.png')  # now it's printing 4 figure!!! 
	plt.show()

plt.hold(True)
w = plt.plot(etas, widths, 'b')
h = plt.plot(etas, heights, 'g')
s = plt.plot(etas, surfaces, 'r')
plt.xlabel('Eta')
plt.ylabel('Cluster property')
plt.legend(["Width", "Height", "Surface"])
plt.savefig('cluster_diffusion_comparison_clusterProperties.png')
plt.show()

		

