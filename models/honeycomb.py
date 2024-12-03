import numpy as np
import matplotlib.pyplot as plt
from math import gcd
import os

class Honeycomb_Hamiltonian:
    """Honeycomb lattice simulation with Anderson localization and a magnetic field."""

    def __init__(self, length: int, t: float, W: float, phi: float, q: int, save = False):
        """
        Initialize the Honeycomb_Hamiltonian class.

        Parameters:
            length (int): Lattice size.
            t (float): Hopping parameter.
            W (float): Disorder strength.
            phi (float): Magnetic flux per plaquette (in units of flux quantum).
            q (int): Maximum denominator for phi values in Hofstadter butterfly.
        """
        self.L = length  # Lattice dimension
        self.N = 2 * self.L * self.L  # Total number of sites (factor of 2 for two sublattices)
        self.t = t  # Hopping parameter
        self.disorder = W  # Disorder strength
        self.phi = phi  # Magnetic flux per plaquette
        self.max_q = q  # Maximum denominator for phi values
        self.matrix = np.zeros((self.N, self.N), dtype=complex)  # Hamiltonian matrix
        self.lattice_type = 'Honeycomb'

        # Initialize on-site disorder potentials
        self.on_site_potential = np.zeros(self.N)
        self.save = save
        if self.save:
            # Base directory for saving plots
            base_dir = os.path.join('plots', self.lattice_type)
            
            # Determine subdirectory based on disorder state
            if self.disorder == 0:
                sub_dir = os.path.join(base_dir, 'No_Disorder', f'L{self.L}_t{self.t}_phi{self.phi}_q{self.max_q}')
            else:
                sub_dir = os.path.join(base_dir, 'Disorder', f'L{self.L}_t{self.t}_phi{self.phi}_q{self.max_q}_dis{self.disorder}')
            
            # Set the path and ensure the directory exists
            self.path = sub_dir
            
            # Create the directory if it doesn't already exist
            os.makedirs(self.path, exist_ok=True)
    
    
    def saving(self, title, save = False):
        
        if self.save == True and save == True and title != None:
            plt.savefig(os.path.join(self.path, title+'.png'))
            plt.show()
        else:
            plt.show()

    #  Defining and diagonalizing the Hamiltonian for the system

    def disorder_setter(self):
        # Apply on-site disorder potentials.
        self.on_site_potential = self.disorder * (2 * np.random.rand(self.N) - 1)

    def peierls_phase(self, delta_x, delta_y, x, y):
        """
        Calculate the Peierls phase.

        Parameters:
            delta_x (int): Change in x-coordinate between sites.
            delta_y (int): Change in y-coordinate between sites.
            x (int): x-coordinate of the starting site.
            y (int): y-coordinate of the starting site.

        Returns:
            complex: Phase factor to be applied to the hopping term.
        """
        # Using Landau gauge
        # Phase accumulated is phi * x * delta_y
        phase = 2 * np.pi * self.phi * (x * delta_y)
        return np.exp(1j * phase)

    def construct_hamiltonian(self):
        """
        Construct the Hamiltonian matrix with hopping,
        Peierls phases, and disorder.

        Returns:
            list: Eigenvalues and eigenvectors of Hamiltonian matrix.
        """
        self.disorder_setter()
        self.matrix = np.zeros((self.N, self.N), dtype=complex)

        for i, j in np.ndindex((self.L, self.L)):
            n = i * self.L + j
            A = 2 * n    # Sublattice A index
            B = A + 1    # Sublattice B index

            # On-site potentials
            self.matrix[A, A] = self.on_site_potential[A]
            self.matrix[B, B] = self.on_site_potential[B]

            # Hopping from A to B in the same unit cell
            phase = self.peierls_phase(0, 0, i, j)
            self.matrix[A, B] = -self.t * phase
            self.matrix[B, A] = -self.t * np.conj(phase)

            # Horizontal hopping from A to B (delta_x = 1, delta_y = 0)
            i_x = (i + 1) % self.L
            n_x = i_x * self.L + j
            B_x = 2 * n_x + 1
            phase = self.peierls_phase(1, 0, i, j)
            self.matrix[A, B_x] = -self.t * phase
            self.matrix[B_x, A] = -self.t * np.conj(phase)

            # Vertical hopping from A to B (delta_x = 0, delta_y = 1)
            j_y = (j + 1) % self.L
            n_y = i * self.L + j_y
            B_y = 2 * n_y + 1
            phase = self.peierls_phase(0, 1, i, j)
            self.matrix[A, B_y] = -self.t * phase
            self.matrix[B_y, A] = -self.t * np.conj(phase)

        # Constructed Hamiltonian
        self.H = self.matrix

        # Compute eigenvalues and eigenvectors
        self.evals, self.evecs = np.linalg.eigh(self.H)

        return self.evals, self.evecs
    
    def plot_hofstadter_butterfly(self, title = None, save = False):

        if title == None:
            title = 'Hofstadter Butterfly for $\phi = p / '+ str(self.max_q) + '$ and $W = '+ str(self.disorder) + '$'
            path = 'Hofstadter Butterfly'

        # Plot the Hofstadter butterfly
        plt.figure(figsize=(10, 8))
        phis = []
        energies = []

        for q in range(1, self.max_q + 1):
            for p in range(q + 1):
                if gcd(p, q) == 1:
                    self.phi = p / q
                    self.construct_hamiltonian()
                    phis.extend([self.phi] * self.N)
                    energies.extend(self.evals.tolist())

        plt.scatter(phis, energies, s=0.1, color='black')
        plt.xlabel('Magnetic Flux per Plaquette $\\phi$')
        plt.ylabel('Energy $E$')
        plt.title(title)
        plt.grid(True)
        self.saving(path, save)

    def prepare_outputs(self):
        """
        Package all relevant parameters and diagonalization 
        outputs in a tuple to pass onto independent plotting functions.

        Returns:
            tuple: Parameter inputs for plotting functions.
        """
        self.evals, self.evecs = self.construct_hamiltonian()
        
        outputs = (self.L, self.t, self.disorder, self.phi, 
                   self.max_q, self.evals, self.evecs, self.lattice_type)
        
        return outputs

