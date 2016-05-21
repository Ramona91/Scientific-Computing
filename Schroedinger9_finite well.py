
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, exp
from math import factorial, pi, sin
from scipy.linalg import eig, eigh
from numpy.polynomial.hermite import hermval


# define grid
N = 144 # number of grid points
a0, a1 = -10, 10 # grid boundaries
deltaX = abs((a0 - a1))/ (N-1) # grid spacing (distance between gridpoints)
x = np.linspace(a0, a1, N) # N grid points evenly spaced between a0 and a1
H_kin = np.zeros([N,N]) # initialize Hamiltonian matrix as 2D array of zeros
quantum_numbers = list(range(10))
potentials = ['finite']





def set_V():
	V = np.zeros([N,N])
	for j in range(0, int(N/4)):
		V[j,j] = 1
	for j in range(int(N/4), int(N - N/4)):
		V[j,j] = 0
	for j in range(int(N - N/4), N):
		V[j,j] = 1

	return V





def setup_hamiltonian(p): # setup Hamiltonian
	for k in range(0, N):
		if k > 0:
			H_kin[k-1,k] = -0.5/deltaX**2 # diagonal above
		H_kin[k,k] = 1/deltaX**2 # main diagonal
		if k+1 < N:
			H_kin[k+1,k] = -0.5/deltaX**2 # diagonal below

		H = H_kin + V # add finite difference part to potential
		#print(H)

	E, Psi = eigh(H, eigvals_only=False) # E is eigenvalue, Psi is eigenvector

	return(E, Psi)







def plot_eigenvectors(Psi, p): # plot some eigenvectors
	plt.hold(True)
	for n in [0, 1, 2]: # columns in eigenvector matrix
		c=np.zeros(n+1)
		c[n]=1
		plt.plot(x, Psi[:, n], label='$n=%i$' % n) #Psi[:, n]/np.sqrt(deltaX)
	plt.legend(ncol=2, loc = 3)
	plt.xlabel(r'$x$')
	plt.ylabel(r'$\Psi_n(x)$')
	#plt.title('Eigenvectors: %s' % p)
	plt.savefig('eigenvectors_finite.png')
	plt.show()


# plot probability distribution of Psi for numerical and analytic solutions, as well as potential
def plot_probabilityDistribution(Psi, p):
	plt.hold(True)
	plt.axvline(x=a0/2, ymin=0, ymax=1, hold=None, linestyle = '--', color = 'k', label='$V$')
	plt.axvline(x=a1/2, ymin=0, ymax=1, hold=None, linestyle = '--', color = 'k')
	for n in [0, 1, 2]: # columns in eigenvector matrix
		c=np.zeros(n+1)
		c[n]=1
		plt.plot(x, (Psi[:, n])**2, label='$n=%i$' % n) # numerical solution of Psi, each column represent an eigenvector corresponding to a quantum number

	plt.xlabel(r'$x$')
	plt.ylabel(r'$|\Psi_n(x)|^2$')
	#plt.title('Probability: %s' % p)
	plt.legend(loc = 1)
	plt.savefig('probabilityDistributionPsi_finite.png')
	plt.show()



		



for p in potentials:
	V = set_V()
	E, Psi = setup_hamiltonian(p)
	plot_eigenvectors(Psi, p)
	plot_probabilityDistribution(Psi, p)


