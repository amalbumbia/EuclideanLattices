import numpy as np
import matplotlib.pyplot as plt
from math import gcd


class Triangular_Hamiltonian:
    """ Triangular lattice simulation with Anderson localization and magnetic field."""

    def __init__(self, length: int, t: float, W: float, phi: float, q: int):
        """
        Initialize the Triangular_Hamiltonian class.

        Parameters:
            length (int): Lattice size.
            t (float): Hopping parameter.
            W (float): Disorder strength.
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

    def disorder_setter(self):
        """Set random on-site potentials for Anderson localization."""
        self.on_site_potential = self.disorder * (2 * np.random.rand(self.N) - 1)

    def peierls_phase(self, x, y, delta_x, delta_y):
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
        # Magnetic flux per plaquette
        phi = self.phi
        x_i = x
        x_f = (x + delta_x)
        y_i = y
        y_f = (y + delta_y)

        # Average x position during hopping
        x_avg = x_i + delta_x / 2

        # Phase accumulated during hopping
        phase = np.exp(2j * np.pi * phi * x_avg * delta_y)

        return phase

    def construct_hamiltonian(self):
        """
        Construct the Hamiltonian matrix with hopping,
        Peierls phases, and disorder.

        """
        self.disorder_setter()
        self.matrix = np.zeros((self.N, self.N), dtype=complex)

        for x in range(self.L):
            for y in range(self.L):
                n = x * self.L + y

                # On-site potential
                self.matrix[n, n] = self.on_site_potential[n]

                # List of neighbor displacements (all six nearest neighbors)
                neighbor_displacements = [
                    (1, 0),    # Right (+x)
                    (0, 1),    # Up (+y)
                    (-1, 1),   # Up-Left (-x, +y)
                    (-1, 0),   # Left (-x)
                    (0, -1),   # Down (-y)
                    (1, -1)    # Down-Right (+x, -y)
                ]

                for delta_x, delta_y in neighbor_displacements:
                    x_neighbor = (x + delta_x) % self.L
                    y_neighbor = (y + delta_y) % self.L
                    n_neighbor = x_neighbor * self.L + y_neighbor

                    # Calculate Peierls phase
                    phase = self.peierls_phase(x, y, delta_x, delta_y)

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