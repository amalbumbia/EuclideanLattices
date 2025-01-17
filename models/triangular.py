import numpy as np
import matplotlib.pyplot as plt
from math import gcd
import os

class Triangular_Hamiltonian:
    """ Triangular lattice simulation with Anderson localization and magnetic field """

    def __init__(self, t: float, W: float, phi: float, max_q: int, save = False):
        """
        Initialize Triangular lattice. 
        
        Parameters:
            t (float): Hopping parameter
            W (float): Disorder strength
            phi (float): Magnetic flux per plaquette (in units of flux quantum).
            q (int): Maximum denominator for phi values in Hofstadter butterfly.
            save (bool): If True, save plots to disk
        """
        self.t = t  # Hopping parameter
        self.disorder = W  # Disorder strength
        self.phi = phi  # Magnetic flux per plaquette
        self.max_q = max_q  # Maximum denominator for phi values
        self.save = save

        # Directory handling if saving
        self.lattice_type = 'Triangular'
        if self.save:
            # Base directory for saving plots
            base_dir = os.path.join('plots', self.lattice_type)
            
            # Determine subdirectory based on disorder state
            if self.disorder == 0:
                sub_dir = os.path.join(base_dir, 'No_Disorder', f't{self.t}_phi{self.phi}_q{self.max_q}')
            else:
                sub_dir = os.path.join(base_dir, 'Disorder', f't{self.t}_phi{self.phi}_q{self.max_q}_dis{self.disorder}')
            
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

    def disorder_setter(self, N):
        # Incorporate the disorder parameter into matrix elements as an on-site disorder potential
        return self.disorder * (2 * np.random.rand(N) - 1)

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

    def construct_hamiltonian(self, p, q):
        """
        Construct the Hamiltonian matrix with hopping,
        Peierls phases, and disorder.

        Returns:
            evals (1D np.array): Eigenvalues of the Hamiltonian
            evecs (2D np.array): Corresponding eigenvectors
        """
        L = q # Set the matrix dimension equal to q (denominator of phi)
        N = L * L # Dimension of adjacency matrix
        self.phi = p / q  # Flux per plaquette

        # Initialize Hamiltonian
        on_site_disorder = self.disorder_setter(N)
        matrix = np.zeros((N, N), dtype=complex)

        for i in range(L):
            for j in range(L):
                n = i * L + j

                # On-site potential
                matrix[n, n] = on_site_disorder[n]

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
                    i_neighbor = (i + delta_i) % L
                    j_neighbor = (j + delta_j) % L
                    n_neighbor = i_neighbor * L + j_neighbor

                    # Calculate Peierls phase
                    phase = self.peierls_phase(i, j, delta_i, delta_j)

                    # Add hopping term
                    matrix[n, n_neighbor] += -self.t * phase

        # Ensure Hamiltonian is Hermitian
        H = (matrix + matrix.conj().T) / 2

        # Compute eigenvalues and eigenvectors
        self.evals, self.evecs = np.linalg.eigh(H)

        return self.evals, self.evecs

    """Defining plotting functions dependent on matrix construction"""

    def plot_hofstadter_butterfly(self, title=None, save=False):
        # Plotting the Hofstadter butterfly
        if title is None:
            title = rf'Hofstadter Butterfly for $\phi = p / '+ str(self.max_q) + '$ and $W = '+ str(self.disorder) + '$'
            path = "Hofstadter_Butterfly"

        phis = []
        energies = []

        for q in range(1, self.max_q + 1):
            for p in range(q + 1):
                
                if gcd(p, q) == 1:
                    evals, _ = self.construct_hamiltonian(p, q) # Reconstruct hamiltonian for each allowed phi
                    
                    phis.extend([p / q] * len(evals))
                    energies.extend(evals.tolist())

        plt.figure(figsize=(10, 8))
        plt.scatter(phis, energies, s=0.1, color='black')
        plt.xlabel(r'Flux per plaquette $\phi = \frac{p}{q}$')
        plt.ylabel(r'Energy $E$')
        plt.title(title)
        plt.grid(True)

    def prepare_outputs(self):
        """
        Package all relevant parameters and diagonalization 
        outputs in a tuple to pass onto independent plotting functions.

        Returns:
            tuple: Parameter inputs for plotting functions.
        """
        self.evals, self.evecs = self.construct_hamiltonian()
        
        outputs = (self.t, self.disorder, self.phi, 
                   self.max_q, self.evals, self.evecs, self.lattice_type)
        
        return outputs