import numpy as np
import matplotlib.pyplot as plt
from math import gcd
import os

class Triangular_Hamiltonian:
    """ Triangular lattice simulation with Anderson localization and magnetic field """

    def __init__(self, mode: str, t: float, W: float, max_q: int = None, L: int = None, p: int = 0, q: int = 0, save=False):
        """
        Initialize Triangular lattice. 
        
        Parameters:
            mode (str): 'butterfly' for Hofstadter butterfly, 'specific' for modeling a specific system.
            t (float): Hopping parameter.
            W (float): Disorder strength.
            max_q (int, optional): Maximum denominator for phi values in Hofstadter butterfly.
            L (int, optional): Lattice size (L x L for the honeycomb lattice).
            p (int, optional): Numerator of the flux fraction (specific mode).
            q (int, optional): Denominator of the flux fraction (specific mode).
            save (bool, optional): If True, save plots to disk.
        """
        self.t = t  # Hopping parameter
        self.disorder = W  # Disorder strength
        self.save = save
        self.mode = mode.lower()

        if self.mode == 'butterfly':
            # Hofstadter butterfly mode
            if max_q is None:
                raise ValueError("max_q must be specified in 'butterfly' mode.")
            self.max_q = max_q
            self.L = 0  # Lattice size is determined by q for each flux calculation
            self.p = 0
            self.q = 0
            self.phi = 0
            self.N = 0

        elif self.mode == 'specific':
            # Specific system mode
            if L is None:
                raise ValueError("Lattice size (L) must be specified in 'specific' mode.")
            self.L = L
            self.N = self.L*self.L
            self.max_q = None  # Not needed for a specific system
            self.p = p
            self.q = q
            if self.p == 0 and self.q == 0:
                self.phi = 0
            else:
                self.phi = self.p / self.q

        else:
            raise ValueError("Invalid mode. Choose 'butterfly' or 'specific'.")


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

    def construct_hamiltonian(self, p_idx, q_idx):
        """
        Construct the Hamiltonian matrix with hopping,
        Peierls phases, and disorder.

        Returns:
            evals (1D np.array): Eigenvalues of the Hamiltonian
            evecs (2D np.array): Corresponding eigenvectors
        """
        if self.L == 0 :
            L_dim = q_idx # Set the matrix dimension equal to q (denominator of phi)
        else:
            L_dim = self.L

        N_dim = L_dim * L_dim # Dimension of adjacency matrix

        if p_idx == 0 and q_idx == 0:
            self.phi = 0
        else:
            self.phi = p_idx / q_idx

        # Initialize Hamiltonian
        on_site_disorder = self.disorder_setter(N_dim)
        matrix = np.zeros((N_dim, N_dim), dtype=complex)

        for i in range(L_dim):
            for j in range(L_dim):
                n = i * L_dim + j

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
                    i_neighbor = (i + delta_i) % L_dim
                    j_neighbor = (j + delta_j) % L_dim
                    n_neighbor = i_neighbor * L_dim + j_neighbor

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
            title = rf'Hofstadter Butterfly for Triangular Lattice $\phi = p / '+ str(self.max_q) + '$ and $W = '+ str(self.disorder) + '$'
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
        plt.grid(False)

        self.saving(path, save)

    def prepare_outputs(self):
        """
        Package all relevant parameters and diagonalization 
        outputs in a tuple to pass onto independent plotting functions.

        Returns:
            tuple: Parameter inputs for plotting functions.
        """
        self.evals, self.evecs = self.construct_hamiltonian(self.p, self.q)
        
        outputs = (self.t, self.disorder, self.p, self.q,
                   self.max_q, self.L, self.N, self.evals, self.evecs, self.lattice_type)
        
        return outputs