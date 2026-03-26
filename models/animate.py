import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import gcd
class DisorderEvolution:
    """
    Animate the evolution of Hofstadter butterflies with increasing disorder
    for a given lattice class (e.g. Square_Hamiltonian in 'butterfly' mode).
    """

    def __init__(self, outputs, lattice_class, max_W,
                 save=False, frames=51, interval=100, style='default'):
        """
        Parameters
        ----------
        outputs : tuple
            (t, disorder, p, q, max_q, L, N, evals, evecs, lattice_type)
            Typically from model.prepare_outputs(), or constructed manually.
        lattice_class : class
            The Hamiltonian class to use (e.g. Square_Hamiltonian).
        max_W : float
            Maximum disorder strength for the animation.
        save : bool
            Whether to use the model's saving logic (not for the animation file).
        frames : int
            Number of frames in the animation.
        interval : int
            Delay between frames in ms.
        style : str
            Matplotlib style to use (e.g. 'default', 'dark_background').
        """
        # Unpack outputs
        self.t, self.disorder, self.p, self.q, self.max_q, \
            self.L, self.N, self.evals, self.evecs, self.lattice_type = outputs

        self.lattice_class = lattice_class
        self.frames = frames
        self.interval = interval
        self.save = save
        self.max_W = max_W

        # Matplotlib figure/artists
        self.fig = None
        self.ax = None
        self.scatter = None
        self.title = None

        # Set style (use a valid style name)
        plt.style.use(style)

    def _collect_butterfly_data(self, W):
        """
        Compute Hofstadter butterfly data for a given disorder W.
        """
        # Instantiate the model in butterfly mode
        model = self.lattice_class(
            mode='butterfly',
            t=self.t,
            W=W,
            max_q=self.max_q,
            L=None,
            save=False
        )

        phis = []
        energies = []

        # Sweep over p/q just like in plot_hofstadter_butterfly
        for q_val in range(1, model.max_q + 1):
            for p_val in range(q_val + 1):
                if gcd(p_val, q_val) == 1:
                    evals, _ = model.construct_hamiltonian(p_val, q_val)
                    phi_val = p_val / q_val
                    phis.extend([phi_val] * len(evals))
                    energies.extend(evals.tolist())

        return np.array(phis), np.array(energies)

    def _init_plot(self):
        """Initialize the animation figure and first frame."""
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.ax.set_xlabel('Flux per plaquette φ = p/q')
        self.ax.set_ylabel('Energy E')

        # initial W = 0
        phis, energies = self._collect_butterfly_data(W=0.0)

        self.scatter = self.ax.scatter(phis, energies, s=0.2, alpha=0.6, c="pink")
        self.title = self.ax.set_title(
            f'Hofstadter Butterfly for {self.lattice_type} Lattice (W = 0.00)'
        )

        self.ax.set_xlim(0, 1)
        if len(energies) > 0:
            self.ax.set_ylim(energies.min()*1.1, energies.max()*1.1)
        self.ax.grid(False)

        return self.scatter, self.title

    def _update(self, frame):
        """Update function for each animation frame."""
        W = (frame / (self.frames - 1)) * self.max_W
        phis, energies = self._collect_butterfly_data(W)
        offsets = np.column_stack((phis, energies))
        self.scatter.set_offsets(offsets)
        self.title.set_text(
            f'Hofstadter Butterfly for {self.lattice_type} Lattice (W = {W:.2f})'
        )
        return self.scatter, self.title

    def animate(self, save_path=None):
        """
        Run the animation. If save_path is given, save as a GIF using pillow.
        """
        self._init_plot()

        anim = FuncAnimation(
            self.fig,
            self._update,
            frames=self.frames,
            interval=self.interval,
            blit=False
        )

        if save_path:
            anim.save(save_path, writer='pillow', fps=10)
        else:
            plt.show()

        plt.close(self.fig)
