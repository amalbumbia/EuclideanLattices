import numpy as np
import matplotlib.pyplot as plt
from math import gcd
import os

class Honeycomb_Hamiltonian:
    """Honeycomb lattice simulation with Anderson localization and a magnetic field."""
    def __init__(self, mode: str, t: float, W: float, max_q: int = None, L: int = None, p: int = 0, q: int = 0, save=False):
        """
        Initialize Honeycomb lattice.

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
            self.N = 2 * (self.L*self.L)
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
        self.lattice_type = 'Honeycomb'
        if self.save:
            base_dir = os.path.join("plots", self.lattice_type) # Base directory for saving plots

             # Determine subdirectory based on disorder state
             
            if self.disorder == 0:
                sub_dir = os.path.join(base_dir, "No_Disorder",
                                       f"t{self.t}_phi{self.phi}_maxq{self.max_q}")
            else:
                sub_dir = os.path.join(base_dir, "Disorder",
                                       f"t{self.t}_phi{self.phi}_maxq{self.max_q}_W{self.disorder}")
                
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
        # Phase accumulated is phi * x * delta_y
        phase = 2 * np.pi * self.phi * (x * delta_y)
        return np.exp(1j * phase)

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

        N_dim = 2 * (L_dim * L_dim) # Dimension of adjacency matrix

        if p_idx == 0 and q_idx == 0:
            self.phi = 0
        else:
            self.phi = p_idx / q_idx

        # Initialize Hamiltonian
        on_site_disorder = self.disorder_setter(N_dim)
        matrix = np.zeros((N_dim, N_dim), dtype=complex)

        for i, j in np.ndindex((L_dim, L_dim)):
            n = i * L_dim + j
            A = 2 * n    # Sublattice A index
            B = A + 1    # Sublattice B index

            # On-site potentials
            matrix[A, A] = on_site_disorder[A]
            matrix[B, B] = on_site_disorder[B]

            # Hopping from A to B in the same unit cell
            phase = self.peierls_phase(0, 0, i, j)
            matrix[A, B] = -self.t * phase
            matrix[B, A] = -self.t * np.conj(phase)

            # Horizontal hopping from A to B (delta_x = 1, delta_y = 0)
            i_x = (i + 1) % L_dim
            n_x = i_x * L_dim + j
            B_x = 2 * n_x + 1
            phase = self.peierls_phase(1, 0, i, j)
            matrix[A, B_x] = -self.t * phase
            matrix[B_x, A] = -self.t * np.conj(phase)

            # Vertical hopping from A to B (delta_x = 0, delta_y = 1)
            j_y = (j + 1) % L_dim
            n_y = i * L_dim + j_y
            B_y = 2 * n_y + 1
            phase = self.peierls_phase(0, 1, i, j)
            matrix[A, B_y] = -self.t * phase
            matrix[B_y, A] = -self.t * np.conj(phase)

        # Constructed Hamiltonian
        H = matrix

        # Compute eigenvalues and eigenvectors
        self.evals, self.evecs = np.linalg.eigh(H)

        return self.evals, self.evecs
    
    """Defining plotting functions dependent on matrix construction"""

    def plot_hofstadter_butterfly(self, title=None, save=False):
        # Plotting Hofstadter butterfly
        if title is None:
            title = rf'Hofstadter Butterfly for Honeycomb Lattice $\phi = p / '+ str(self.max_q) + '$ and $W = '+ str(self.disorder) + '$'
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