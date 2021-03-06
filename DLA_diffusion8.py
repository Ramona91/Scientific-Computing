# Scientific Computing, assignment 2
# Set1, ex. J
# Object in diffusion space, implemented with SOR

import numpy as np
from numpy import random
import numpy.ma as ma
import matplotlib.pyplot as plt
#import matplotlib as mpl
from math import sqrt
from scipy.special import erfc



nx = 200 # no of interval in x-direction
ny = 200 # # no of interval in x-direction
d = 1.0 # diffusion coefficient
omega = 1.6 #needed for successive over-relaxation method
etha = 2 # determines shape of object




# results for different etas (cluster characteristics)
ethas = [etha]
growthSteps = 0
concentrations = []



cluster = np.zeros((nx+1, ny+1)) # presence of object encoded as 1 in object array, object in centre cell of matrix
cluster[nx, (ny+1)*0.5] = 1 # object starts to grow at the bottom centre

s = np.zeros((nx+1, ny+1)) # concentrations array
s[0,:] = 1 # boundary condition: source
s[nx, :] = 0 # boundary condition: sink






def find_growthCandidates(cluster):
	concentrationCandidates = []
	concentrationSum = 0
	xCandidates = [] # keep track of indeces where objects are located (important to calculate cluster properties)
	yCandidates = []
	for j in range(0, ny+1):
		for i in range(1, nx):
			if (cluster[i,j] != 1) and ((cluster[i, j-1] == 1) or (cluster[i, (j+1)%(ny+1)] == 1) or (cluster[i-1,j] == 1) or (cluster[i+1, j] == 1)):
				c = max((s[i,j])**etha, 0) # growth candidate to the power of eta (part of probability formula)
				concentrationCandidates.append(c)
				yCandidates.append(j)
				xCandidates.append(i)
	concentrationSum = sum(concentrationCandidates)
	if concentrationSum == 0:
		concentrationSum = 10**(-3)
	#print(concentrationCandidates)
	# print(concentrationSum)
	return(concentrationCandidates, concentrationSum, yCandidates, xCandidates)


def calculate_p(concentrationCandidates, concentrationSum):
	probabilities = []
	for c in concentrationCandidates:
		p = (c/concentrationSum)
		probabilities.append(p)
	#print(probabilities)
	return(probabilities)


def add_to_cluster(probabilities, yCandidates, xCandidates, s, cluster):
	for (p, i, j) in zip(probabilities, xCandidates, yCandidates):
		if p > random.uniform():
			cluster[i,j] = 1
			s[i,j] = 0
	#print(cluster)
	return (s, cluster)



def do_analytic(t, s):
	for i in range(0, nx+1):
		s[i,:] = (nx-i)/nx
	return s



def do_SOR(s, cluster):
	compare = 0
	value = 0
	delta = 1
	tolerance = 10**(-5)
	while (delta > tolerance):
		delta = 0
		for j in range(0, ny+1):
			s[0, j] = 1
			for i in range(1, nx):
				value = np.copy(s[i,j])
					
				if cluster[i,j] == 1: # put beginning seed in cluster matrix
					s[i,j] = 0
					#s[0, (ny+1)*0.5] = 0 # not needed if specified above
				else:
					s[i,j] = (omega / 4) * (s[i-1,j] + s[i,j-1] + s[i+1,j] + s[i,(j+1)%(ny+1)]) + (1 - omega) * s[i,j]
					compare = np.abs(s[i,j] - value)

				if (compare > delta):
					delta = compare

	return(s)




def calculate_clusterProperties(cluster):
	xIndex = []
	yIndex = []
	for j in range(0, ny+1):
		for i in range(0, nx+1):
			if cluster[i,j] == 1:
				xIndex.append(i+1)
				yIndex.append(j+1)

	#print(xIndex)
	width = max(yIndex) - min(yIndex)
	height = max(xIndex)
	surface = len(xIndex)

	print("width:", width)
	print("height:", height)
	print("surface:", surface)



def plot_dla(concentrations, ethas):
	for (c,e) in zip(concentrations, ethas):
		cmap = plt.cm.autumn
		cmap.set_bad(color='black')
		im = plt.imshow(c, origin = 'upper', interpolation='none', cmap=cmap) #origin='lower'
		plt.colorbar(im)
		plt.xlabel('x')
		plt.ylabel('y')
		#ax.invert_yaxis()
		plt.savefig('cluster_diffusion_eta' + str(e) + '.png')  # now it's printing 4 figure!!! 
		plt.show()





while growthSteps > (-1):
	if growthSteps == 0:
		s = do_analytic(1, s) # t = 1
		#s = do_SOR(s, cluster)
		print(s)
	else:
		(concentrationCandidates, concentrationSum, yCandidates, xCandidates) = find_growthCandidates(cluster)
		probabilities = calculate_p(concentrationCandidates, concentrationSum)
		(s, cluster) = add_to_cluster(probabilities, yCandidates, xCandidates, s, cluster)
		s = do_SOR(s, cluster)

	if growthSteps%100 == 0:
		print(growthSteps)



	if np.any(cluster[1,:] == 1) or np.any(cluster[1:nx-1, 0] == 1) or np.any(cluster[1:nx-1, ny] == 1): # stop condition
		s = ma.masked_where(s == 0, s)
		concentrations.append(s)
		print(growthSteps)
		break
	else:
		growthSteps += 1

calculate_clusterProperties(cluster)
plot_dla(concentrations, ethas)












		

