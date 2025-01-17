import numpy as np
import matplotlib.pyplot as plt
from math import gcd
import os


class Square_Hamiltonian:
    """ Square lattice simulation with Anderson localization and a magnetic field"""
    def __init__(self, t: float, W: float, phi: float, max_q: int, save=False):
        """
        Initialize Square Lattice.
        
        Parameters:
            t (float): Hopping parameter
            W (float): Disorder strength
            phi (float): Magnetic flux per plaquette (in units of flux quantum).
            q (int): Maximum denominator for phi values in Hofstadter butterfly.
            save (bool): If True, save plots to disk
        """
        self.t = t
        self.disorder = W
        self.phi = phi
        self.max_q = max_q
        self.save = save

         # Initialize boundary fluxes - these do not contribute to the flux plaquett
        self.phi_x = 0.0  # Flux through x direction
        self.phi_y = 0.0  # Flux through y direction
        
        # Directory handling if saving
        self.lattice_type = "Square"
        if self.save:
            base_dir = os.path.join("plots", self.lattice_type) # Base directory for saving plots

             # Determine subdirectory based on disorder state
             
            if self.disorder == 0:
                sub_dir = os.path.join(base_dir, "No_Disorder",
                                       f"t{self.t}_maxq{self.max_q}")
            else:
                sub_dir = os.path.join(base_dir, "Disorder",
                                       f"t{self.t}_maxq{self.max_q}_W{self.disorder}")
                
            # Set the path and ensure the directory exists
            self.path = sub_dir

            # Create the directory if it doesn't already exist
            os.makedirs(self.path, exist_ok=True)


    def saving(self, title, save=False):
       
        if self.save and save and title is not None:
            plt.savefig(os.path.join(self.path, title + ".png"))
            plt.show()
        else:
            plt.show()

    """ Defining and diagonalizing the Hamiltonian for the system """

    def disorder_setter(self, N):
        # Incorporate the disorder parameter into matrix elements as an on-site disorder potential
        return self.disorder * (2 * np.random.rand(N) - 1)

    def peierls_phase(self, i, j, L, direction):
        """
        Calculate the Peierls phase for hopping between sites.

        Parameters:
            i (int): x-index of the starting site.
            j (int): y-index of the starting site.
            direction (str): 'x' for horizontal hopping, 'y' for vertical hopping.

        Returns:
            complex: exp(i * phase)
        """
        if direction == 'x':
            # Hopping in the x-direction
            phase = 0.0
            if (i + 1) >= L:
                # Boundary hopping in x-direction  
                phase += 2 * np.pi * self.phi_x
            return np.exp(1j * phase)

        elif direction == 'y':
            # Hopping in the y-direction
            phase = 2 * np.pi * self.phi * i
            if (j + 1) >= L:
                # Boundary hopping in y-direction  
                phase += 2 * np.pi * self.phi_y
            return np.exp(1j * phase)


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

        for i, j in np.ndindex((L, L)):
            n = i * L + j  # Current site index

            # On-site potential
            matrix[n, n] = on_site_disorder[n]

            # Hopping in x-direction (to site (i+1, j))
            m_x = ((i + 1) % L) * L + j
            phase_x = self.peierls_phase(i, j, L, 'x')
            matrix[n, m_x] = -self.t * phase_x

            # Hopping in y-direction (to site (i, j+1))
            m_y = i * L + ((j + 1) % L)
            phase_y = self.peierls_phase(i, j, L, 'y')
            matrix[n, m_y] = -self.t * phase_y

        # Ensure the Hamiltonian is Hermitian
        H = matrix + matrix.conj().T

        # Compute eigenvalues and eigenvectors
        evals, evecs = np.linalg.eigh(H)
        return evals, evecs
    
    """Defining plotting functions dependent on matrix construction"""

    def plot_hofstadter_butterfly(self, title=None, save=False):
        # Plotting the Hoftsadter butterfly
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

    def prepare_outputs(self, p, q):
        """
        Package all relevant parameters and diagonalization 
        outputs in a tuple to pass onto independent plotting functions.

        Returns:
            tuple: Parameter inputs for plotting functions.
        """
        self.evals, self.evecs = self.construct_hamiltonian(p, q)
        
        outputs = (self.t, self.disorder, self.phi, 
                   self.max_q, self.evals, self.evecs, self.lattice_type)
        
        return outputs