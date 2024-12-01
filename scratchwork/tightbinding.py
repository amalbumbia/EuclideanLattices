import numpy as np
from matplotlib import pyplot as plt
# Tight Binding Chain

L = 100 # n. of sites
H = np.zeros((L,L), dtype = complex) # Hamiltonian
for i in range(L):
  j = (i+1) % L # neighbor site, with periodic boundary conditions
  H[i,j] = H[j,i] = 1 # hopping amplitude
evals, evecs = np.linalg.eigh(H) # diagonalize a Hermitian matrix

plt.plot(evals, '.')
plt.ylabel(r'$E_i$')
plt.xlabel(r'$i$')
plt.show()
