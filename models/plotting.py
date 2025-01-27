import numpy as np
import matplotlib.pyplot as plt
import os

class Plotting_Functions:
    
    def __init__(self, outputs, save = False):
        self.t = outputs[0]  # Hopping parameter
        self.disorder = outputs[1]  # Disorder strength
        self.p = outputs[2]
        self.q = outputs[3]
        
        if self.p == 0 and self.q == 0:
            self.phi = 0
        else:
            self.phi = self.p / self.q

        self.max_q = outputs[4]  # Maximum denominator for phi values
        self.L = outputs[5]
        self.N = outputs[6]
        self.evals = outputs[7]
        self.evecs = outputs[8]
        self.lattice_type = outputs[9]
        self.save = save


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
        
        if save == True and title != None:
            plt.savefig(os.path.join(self.path, title))
            plt.show()
        else:
            plt.show()

    """ Basic plotting functions """

    def plot_evals(self, title = None, save = False):
        
        if title == None:
            title = 'Eigenvalues of the ' + self.lattice_type + ' Lattice Hamiltonian'

        # Plot eigenvalues of hamiltonian matrix
        legend = f'L={self.L}, t={self.t}, W={self.disorder}, $\phi$={self.phi}'
        plt.plot(self.evals, '.')
        plt.ylabel(r'Eigenvalues $E_i$')
        plt.xlabel('Index $i$')
        plt.title(title)
        plt.legend([legend])
        plt.grid(True)

        self.saving(title, save)


    def plot_evec(self, title = None, save = False):

        if title == None:
            title = 'Arbitrary Eigenvector'

        # Plot some eigenvector in the middle of the spectrum in the presence of disorder
        #self.psi = self.evecs[:,self.L//2] # Some eigenvector in the middle of the spectrum
        #self.psi = self.evecs[:, min(self.L//2, self.evecs.shape[1]-1)]
        middle_index = min(len(self.evals) // 2, self.evecs.shape[1] - 1)
        self.psi = self.evecs[:, middle_index]
        fig, ax = plt.subplots(2,1,sharex=True)
        ax[0].plot(np.abs(self.psi)**2)
        ax[1].semilogy(np.abs(self.psi)**2)
        ax[1].set_xlabel('x')
        ax[0].set_ylabel(r'$ |\psi(x)|^2$')
        ax[1].set_ylabel(r'$ |\psi(x)|^2$')

        legend = f'L={self.L}, t={self.t}, W={self.disorder}, $\phi$={self.phi}'
        plt.title('Arbitrary Eigenvector')
        plt.legend([legend])
        plt.grid(True)

        self.saving(title, save)

    def plot_pr(self, title = None, save = False):

        if title == None:
            title = 'Localization Properties of the '+self.lattice_type+' Lattice'

        # Plot Participation Ratio
        self.PR = 1./np.sum(np.abs(self.evecs)**4, axis=0) # 'evecs' is a matrix of $\psi_i(x)$ amplitudes, 1st axis is x. This does the sum over x.

        legend = f'L={self.L}, t={self.t}, W={self.disorder}, $\phi$={self.phi}'
        plt.plot(self.evals, self.PR, 'o')
        plt.xlabel('Energy $E$')
        plt.ylabel('Inverse Participation Ratio (IPR)')
        plt.title(title)
        plt.legend([legend])
        plt.grid(True)

        self.saving(title, save)

    def plot_density_of_states(self, sigma=0.1, num_points=1000, title = None, save = False):
        """
        Plot the density of states.

        Parameters:
            sigma (float): Standard deviation for Gaussian broadening.
            num_points (int): Number of points in the energy grid.
        """
        if title == None:
            title = 'Density of States vs Energy for '+ self.lattice_type + ' Lattice'

        energy_min = np.min(self.evals) - 1
        energy_max = np.max(self.evals) + 1
        E_vals = np.linspace(energy_min, energy_max, num_points) # Artificial energy space seperate to eigenvalues
        dos = np.zeros_like(E_vals)

        for E_n in self.evals:
            dos += np.exp(-((E_vals - E_n) ** 2) / (2 * sigma ** 2)) / (np.sqrt(2 * np.pi) * sigma) # Using gaussian broadening

        legend = f'L={self.L}, t={self.t}, W={self.disorder}, $\phi$={self.phi}'
        plt.figure(figsize=(8, 6))
        plt.plot(E_vals, dos)
        plt.xlabel('Energy $E$')
        plt.ylabel('Density of States $g(E)$')
        plt.title(title)
        plt.legend([legend])
        plt.grid(True)
        
        self.saving(title, save)


    def plot_hall_conductance(self, num_points=100, title = None, save = False):
        """
        Compute and plot the Hall conductance.

        Parameters:
            num_points (int): Number of points in the energy grid.
        """
        if title == None:
            title = 'Hall Conductance vs Energy for '+ self.lattice_type + ' Lattice'

        # Exploiting density of states N = \int g(E) ; then let hall conductance = N e^2/h
        h = 6.62607015e-34  # Planck's constant (Js)
        e = 1.602176634e-19  # Elementary charge (C)

        energy_min = np.min(self.evals) - 1
        energy_max = np.max(self.evals) + 1
        energies = np.linspace(energy_min, energy_max, num_points) # Artificial linspace 
        hall_conductances = np.zeros_like(energies)

        # Cumulative density of states
        cumulative_dos = np.array([np.sum(self.evals < E) for E in energies]) / self.N

        # Hall conductance (in units of e^2/h)
        hall_conductances = cumulative_dos * (e ** 2 / h)

        legend = f'L={self.L}, t={self.t}, W={self.disorder}, $\phi$={self.phi}'
        plt.figure(figsize=(8, 6))
        plt.plot(energies, hall_conductances)
        plt.xlabel('Energy $E$')
        plt.ylabel('Hall Conductance $\sigma_{xy}$ (S)')
        plt.title(title)
        plt.legend([legend])
        plt.grid(True)

        self.saving(title, save)    

    def plot_hall_conductance_test(self, num_points=100, title=None, save=False, eta=None):
        """
        Compute and plot the Hall conductance using the Kubo formula, with improved stability.
        """
        if title == None:
            title = 'Hall Conductance vs Energy for '+ self.lattice_type + ' Lattice'

        H = self.H
        L = self.L
        N = L * L

        evals, evecs = np.linalg.eigh(H)

        if eta is None:
            energy_range = evals.max() - evals.min()
            eta = 0.05 * energy_range / N if N > 1 else 1e-3

        x_coords = np.array([k // L for k in range(N)], dtype=float)
        y_coords = np.array([k % L for k in range(N)], dtype=float)
        X = np.diag(x_coords)
        Y = np.diag(y_coords)

        vx = 1j * (H @ X - X @ H)
        vy = 1j * (H @ Y - Y @ H)

        for i in range(N):
            idx = np.argmax(np.abs(evecs[:, i]))
            phase = np.exp(-1j * np.angle(evecs[idx, i]))
            evecs[:, i] *= phase

        vx_eig = evecs.conj().T @ vx @ evecs
        vy_eig = evecs.conj().T @ vy @ evecs

        dE = evals[None, :] - evals[:, None]
        np.fill_diagonal(dE, np.inf)
        denominator = dE**2 + eta**2
        cross_term = np.imag(vx_eig * vy_eig.conj())

        E_min, E_max = evals.min(), evals.max()
        energies = np.linspace(E_min - 0.5, E_max + 0.5, num_points)
        hall_conductances = np.zeros_like(energies)
        prefactor = 1.0 / (2.0 * np.pi)

        for i, E in enumerate(energies):
            f = (evals < E).astype(float)
            f_diff = f[:, None] - f[None, :]
            integrand = cross_term * f_diff / denominator
            hall_conductances[i] = prefactor * np.sum(integrand)

        plt.figure(figsize=(8, 6))
        plt.plot(energies, hall_conductances, lw=1.5)
        plt.xlabel('Energy $E$')
        plt.ylabel('$\\sigma_{xy}$ ($e^2/h$)')
        plt.title(title)
        plt.grid(True)

        self.saving(title, save)   