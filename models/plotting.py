import numpy as np
import matplotlib.pyplot as plt
import os

class Plotting_Functions:
    
    def __init__(self, outputs, save = False):
        
        
        self.L = outputs[0]  # Lattice dimension
        self.N = self.L * self.L  # Total number of sites
        self.t = outputs[1]  # Hopping parameter
        self.disorder = outputs[2]  # Disorder strength
        self.phi = outputs[3]  # Magnetic flux per plaquette
        self.max_q = outputs[4]  # Maximum denominator for phi values
        self.evals = outputs[5]
        self.evecs = outputs[6]
        self.lattice_type = outputs[7]
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
        self.psi = self.evecs[:,self.L//2] # Some eigenvector in the middle of the spectrum

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
