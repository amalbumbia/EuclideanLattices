import numpy as np
import matplotlib.pyplot as plt
from math import gcd


class Triangular_Hamiltonian:
    """ Triangular lattice simulation with Anderson localization and magnetic field """

    def __init__(self, length: int, t: float, W: float, phi: float, q: int):
        """
          Initialize Triangular_Hamiltonian class.

          Parameters:
              length (int): Lattice size (L x L).
              t (float): Hopping parameter.
              W (float): Disorder parameter.
              phi (float): Magnetic flux per plaquette (in units of flux quantum).
              q (int): Maximum denominator for phi values in Hofstadter butterfly.
        """

        self.L = length
        self.N = self.L * self.L  # Total number of sites
        self.t = t  # Hopping parameter
        self.disorder = W  # Disorder strength
        self.phi = phi  # Magnetic flux per plaquette
        self.max_q = q  # Maximum denominator for phi values
        self.matrix = np.zeros((self.N, self.N), dtype=complex)
        self.on_site_potential = np.zeros(self.N)

    """ Defining and diagonalizing the Hamiltonian for the system """

    def disorder_setter(self):
        # Incorporate the disorder parameter into matrix elements as an on-site disorder potential
        self.on_site_potential = self.disorder * (2 * np.random.rand(self.N) - 1)

    def peierls_phase(self, i, j, delta_i, delta_j):
        """Calculate Peierls phase factor for hopping between sites.

        Parameters:
        i (int): x-index of starting site
        j (int): y-index of starting site
        delta_i (int): Change in x-coordinate between sites.
        delta_j (int): Change in y-coordinate between sites.

        Returns:
        complex: Phase factor for the hopping term
        """

        # Using the Landau gauge
        # Average x position during hop
        i_avg = i + delta_i / 2

        # Phase accumulated during hop
        phase = np.exp(2j * np.pi * self.phi * i_avg * delta_j)

        return phase

    def construct_hamiltonian(self):
        """Construct the Hamiltonian matrix with hopping and disorder terms."""
        self.disorder_setter()
        self.matrix = np.zeros((self.N, self.N), dtype=complex)

        for i in range(self.L):
            for j in range(self.L):
                n = i * self.L + j

                # On-site potential
                self.matrix[n, n] = self.on_site_potential[n]

                # List of the six nearest neighbors
                neighbors = [
                    (1, 0),    # Right (+i)
                    (0, 1),    # Up (+j)
                    (-1, 1),   # Up-Left (-i, +j)
                    (-1, 0),   # Left (-i)
                    (0, -1),   # Down (-j)
                    (1, -1)    # Down-Right (+i, -j)
                ]

                for delta_i, delta_j in neighbors:
                    i_neighbor = (i + delta_i) % self.L
                    j_neighbor = (j + delta_j) % self.L
                    n_neighbor = i_neighbor * self.L + j_neighbor

                    # Calculate Peierls phase
                    phase = self.peierls_phase(i, j, delta_i, delta_j)

                    # Add hopping term
                    self.matrix[n, n_neighbor] += -self.t * phase

        # Ensure Hamiltonian is Hermitian
        self.H = (self.matrix + self.matrix.conj().T) / 2

        # Compute eigenvalues and eigenvectors
        self.evals, self.evecs = np.linalg.eigh(self.H)
        
    def plot_hofstadter_butterfly(self):
        # Plot Hofstadter butterfly.
        plt.figure(figsize=(10, 8))
        phis = []
        energies = []

        original_phi = self.phi  # Save original phi value

        for q in range(1, self.max_q + 1):
            for p in range(q):
                if gcd(p, q) == 1:
                    self.phi = p / q
                    self.construct_hamiltonian()
                    phis.extend([self.phi] * self.N)
                    energies.extend(self.evals.tolist())

        self.phi = original_phi  # Restore original phi value

        plt.scatter(phis, energies, s=0.1, color='black')
        plt.xlabel('Flux per Plaquette $\\phi$')
        plt.ylabel('Energy $E$')
        plt.title(f'Hofstadter Butterfly for $\\phi = p / {self.max_q}$ and $W = {self.disorder}$')
        plt.grid(True)
        plt.show()
    def prepare_outputs(self):
        """
        Package all relevant parameters and diagonalization 
        outputs in a tuple to pass onto independent plotting functions.

        Returns:
            tuple: Parameter inputs for plotting functions.
        """
        self.evals, self.evecs = self.construct_hamiltonian()
        
        outputs = (self.L, self.t, self.disorder, self.phi, 
                   self.max_q, self.evals, self.evecs)
        
        return outputs