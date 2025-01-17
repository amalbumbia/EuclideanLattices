def compute_thouless_conductance(self):
        """
        Compute and plot the Thouless conductance.
        """
        # Shift in eigenvalues due to a small flux perturbation
        original_evals = self.evals.copy()

        # Apply a small flux perturbation
        delta_phi = 1e-5
        self.phi_x += delta_phi  # Perturb flux in x-direction
        self.construct_hamiltonian()
        perturbed_evals = self.evals.copy()

        # Reset phi_x
        self.phi_x -= delta_phi
        self.construct_hamiltonian()

        # Compute energy differences
        delta_E = np.abs(perturbed_evals - original_evals)

        # Level spacings
        level_spacings = np.diff(original_evals)
        mean_spacing = np.mean(level_spacings)

        # Thouless conductance
        g_thouless = delta_E / mean_spacing

        # Plotting
        plt.figure(figsize=(8, 6))
        plt.plot(original_evals, g_thouless, '.', markersize=5)
        plt.xlabel('Energy $E$')
        plt.ylabel('Thouless Conductance $g_{xx}$')
        plt.title('Thouless Conductance vs Energy')
        plt.grid(True)
        plt.show()