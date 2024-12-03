import numpy as np
import matplotlib.pyplot as plt
from math import gcd
import os

class Square_Hamiltonian:
    """ Square lattice simulation with Anderson localization and a magnetic field"""

    def __init__(self, length: int, t: float, W: float, phi: float, q: int, save = False):
        """
          Initialize Square_Hamiltonian class.

          Parameters:
              length (int): Lattice size (L x L).
              t (float): Hopping parameter.
              W (float): Disorder parameter.
              phi (float): Magnetic flux per plaquette (in units of flux quantum).
              q (int): Maximum denominator for phi values in Hofstadter butterfly.
          """

        self.L = length  # Lattice dimension
        self.N = self.L * self.L  # Total number of sites
        self.t = t  # Hopping parameter
        self.disorder = W  # Disorder strength
        self.phi = phi  # Magnetic flux per plaquette
        self.max_q = q  # Maximum denominator for phi values
        self.matrix = np.zeros((self.N, self.N), dtype=complex)  # Hamiltonian matrix
        self.lattice_type = "Square"
        # Initialize on-site disorder potentials
        self.on_site_potential = np.zeros(self.N)

        # Initialize boundary fluxes - these do not contribute to the flux plaquett
        self.phi_x = 0.0  # Flux through x direction
        self.phi_y = 0.0  # Flux through y direction
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


    """ Defining and diagonalizing the Hamiltonian for the system """

    def disorder_setter(self):
        # Incorporate the disorder parameter into matrix elements as an on-site disorder potential
        self.on_site_potential = self.disorder * (2 * np.random.rand(self.N) - 1)

    def peierls_phase(self, i, j, direction):
        """
        Calculate the Peierls phase for hopping between sites.

        Parameters:
            i (int): x-index of the starting site.
            j (int): y-index of the starting site.
            direction (str): 'x' for horizontal hopping, 'y' for vertical hopping.

        Returns:
            float: Phase factor for hopping term.
        """
        # Using Landau gauge
        if direction == 'x':
            # Hopping in the x-direction
            phase = 0.0
            if (i + 1) >= self.L:
                # Boundary hopping in x-direction
                phase += 2 * np.pi * self.phi_x
            return np.exp(1j * phase)
        elif direction == 'y':
            # Hopping in the y-direction
            phase = 2 * np.pi * self.phi * i
            if (j + 1) >= self.L:
                # Boundary hopping in y-direction
                phase += 2 * np.pi * self.phi_y
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
            n = i * self.L + j  # Current site index

            # On-site potential
            self.matrix[n, n] = self.on_site_potential[n]

            # Hopping in x-direction (to site (i+1, j))
            m_x = ((i + 1) % self.L) * self.L + j
            phase_x = self.peierls_phase(i, j, 'x')
            self.matrix[n, m_x] = -self.t * phase_x

            # Hopping in y-direction (to site (i, j+1))
            m_y = i * self.L + (j + 1) % self.L
            phase_y = self.peierls_phase(i, j, 'y')
            self.matrix[n, m_y] = -self.t * phase_y

        # Ensure the Hamiltonian is Hermitian
        self.H = self.matrix + self.matrix.conj().T

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
                    self.construct_hamiltonian() # Reconstruct hamiltonian for each allowed phi
                    phis.extend([self.phi] * self.N)
                    energies.extend(self.evals.tolist())

        plt.scatter(phis, energies, s=0.1, color='black')
        plt.xlabel('Flux per Plaquette $\phi$')
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

