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
potentials = ['infinite', 'parabolic']





def V(x, p): # caluclate potential
	if p == 'infinite':
		return 0*x 

	if p == 'parabolic':
		return 0.5*x**2




def setup_hamiltonian(p): # setup Hamiltonian
	for k in range(0,N):
		if k > 0:
			H_kin[k-1,k] = -0.5/deltaX**2 # diagonal above
		H_kin[k,k] = 1/deltaX**2 # main diagonal
		if k+1 < N:
			H_kin[k+1,k] = -0.5/deltaX**2 # diagonal below

		H_potential = np.diag(V(x, p))
		H = H_kin + H_potential # add finite difference part to potential
		#print(H)

	E, Psi = eigh(H, eigvals_only=False) # E is eigenvalue, Psi is eigenvector

	return(E, Psi)


def analytic(p):
	if p == 'infinite':
		psi_analytic = []
		psi2_analytic = []
		for n in quantum_numbers:
			psi_one_n = []
			for i in x:
				Psi = sqrt(2/(a1*2)) * sin(n*pi*i/(a1*2))
				psi_one_n.append(Psi)

				if i == x[-1]:
					psi_analytic.append(psi_one_n)
		psi2_analytic = [[e*e for e in m] for m in psi_analytic]
		return(psi_analytic, psi2_analytic)


	if p == 'parabolic':
		psi_analytic = []
		psi2_analytic = []
		for n in quantum_numbers:
			psi_one_n = []
			for i in x:
				Psi = - sqrt(8/(a1*2)) * sin(n*pi*i/(a1*2))
				psi_one_n.append(Psi)

				if i == x[-1]:
					psi_analytic.append(psi_one_n)
		psi2_analytic = [[e*e for e in m] for m in psi_analytic]
		return(psi_analytic, psi2_analytic)






def plot_eigenvalues(E): #plot first eigenvalues
	n=10
	plt. plot(np.array(range(0, n)), np.array(range(0, n))+0.5, 'r', label='analytic')
	plt. plot(np.array(range(0, n)), np.sort(E)[0:n], 'b' ,label='numerical')
	plt.xlabel('index')
	plt.ylabel('energy')
	plt.legend() #bbox_to_anchor=(0., 1.02, 1., .102), loc=3, borderaxespad=0.
	plt.savefig('eigenvalues.png')
	plt.show()


def plot_eigenvectors(Psi, psi_analytic, p): # plot some eigenvectors
	plt.hold(True)
	for n in [0, 1, 2]: # columns in eigenvector matrix
		c=np.zeros(n+1)
		c[n]=1
		plt.plot(x, Psi[:, n], label='$n_{num}=%i$' % n) #Psi[:, n]/np.sqrt(deltaX)
		plt.plot(x, psi_analytic[n], label='$n_{anal}=%i$' % n)
	plt.legend(ncol=2, loc = 4)
	plt.xlabel(r'$x$')
	plt.ylabel(r'$\Psi_n(x), \Psi_{exact}(x)$')
	#plt.title('Eigenvectors: %s' % p)
	plt.savefig('eigenvectors_%s.png' % p)
	plt.show()


# plot probability distribution of Psi for numerical and analytic solutions, as well as potential
def plot_probabilityDistribution(Psi, psi2_analytic, p):
	plt.hold(True)
	plt.plot(x, V(x, p)*0.001, linestyle = '--', color = 'k', label='$V$') # potential
	for n in [0, 1, 2]: # columns in eigenvector matrix
		c=np.zeros(n+1)
		c[n]=1
		plt.plot(x, (Psi[:, n])**2, label='$n=%i$' % n) # numerical solution of Psi, each column represent an eigenvector corresponding to a quantum number
		#plt.plot(x, (psi2_analytic[n]), label='$n_{anal}=%i$' % n)
		
		#plt.plot(hermval(x, c)*exp(-0.5*x**2)/sqrt(sqrt(pi)*2**n*factorial(n))**2) # analytic solution
		#plt.plot(x, )
	plt.xlabel(r'$x$')
	plt.ylabel(r'$|\Psi_n(x)|^2$')
	#plt.title('Probability: %s' % p)
	plt.legend()
	plt.savefig('probabilityDistributionPsi_%s.png' % p)
	plt.show()


def plot_error(error, p): # plot error of some eigenvectors
	for n in [0, 1, 2]:
		plt.plot(x, error[n], label='$n=%i}$' % n)
		# c=np.zeros(n+1)
		# c[n]=1
		#plt.plot(x, (Psi[:, n]/sqrt(deltaX)**2) -(hermval(x, c)*exp(-0.5*x**2)/sqrt(sqrt(pi)*2**n*factorial(n))**2), label='$n=%i$' % n)
		# plt.plot(x, (Psi[:, n]**2) - (psi2_analytic[n]), label='$n=%i$' % n)
	plt.legend()
	plt.xlabel(r'$x$')
	plt.ylabel(r'$|\Psi_n(x)|^2-|\Psi_{n,\mathrm{exact}}(x)|^2$')
	#plt.title('Error: %s' % p)
	plt.savefig('probability_error_%s.png' % p)
	plt.show()



def error(Psi, psi2_analytic, psi_analytic):
	error_prob = []
	error_psi = []
	for n in [0, 1, 2]:
		Psi_list = Psi[:,n].tolist()
		Psi2_list = [i ** 2 for i in Psi_list]
		error_prob_n = [a - b for a, b in zip(psi2_analytic[n], Psi2_list)]
		error_psi_n = [a - b for a, b in zip(psi_analytic[n], Psi_list)]
		error_prob.append(error_prob_n)
		error_psi.append(error_psi_n)
	#print(len(error_all))
	return error_prob, error_psi
		



for p in potentials:
	psi_analytic, psi2_analytic = analytic(p)
	E, Psi = setup_hamiltonian(p)
	error_prob, error_psi = error(Psi, psi2_analytic, psi_analytic)
	plot_eigenvectors(Psi, psi_analytic, p)
	plot_probabilityDistribution(Psi, psi2_analytic, p)
	plot_error(error_prob, p)

