# Euclidean Lattice Simulations

This repository is a work in progress. My aim is to simulate Anderson localization and the quantum hall effect on a few common euclidean lattices (square, honeycomb, triangular, kagome). A primary result is the Hofstadter butterfly for each of these lattice types given different input parameters such as the presence of an on-site disorder potential. 
 
As it stands, the Thouless conductance and chern number additions are still in progress. The kagome lattice model is nearly complete. 

There is an odd quirk in the butterflies due to an additional periodic smattering of points overlayed on top of the correct butterfly structure for each model. I'm trying to figure out exactly why that's happening so I can fix it. 

# Getting Started

To use the notebooks, ensure you have a `conda` environment with:
- `numpy`
- `matplotlib`
- `jupyter`

# Introduction
Condensed matter physics realizes the importance of 2D materials when it comes to the observation of interesting phenomena (ie. optical, magnetic, etc) and topological properties as well as technological applications. 
Different 2D materials have different underlying lattice structures. Some common ones are the square lattice, the honeycomb lattice, triangular/hexagonal lattice, and the kagome lattice. 

These lattices can be described via Hamiltonians encoding how a species could "hop" between adjacent sites (nearest neighbors). A useful way to understand these systems is by diagonalizing the associated Hamiltonians to obtain allowed energies and states, and use these eigenvalues and eigenvectors to compute various parameters. In general, it is useful to plot the eigenvalue spectrum, different eigenvectors, density of states, and the inverse participation ratio. In the case of a magnetic field application perpendicular to the lattice, plotting the Hofstader butterfly is a primary result. Generally, the band structure of these systems is useful to compute, but energy band theory is dependent on the concept of crytal momentum --- a relevant parameter when it comes to describing translational symmetry in the lattice via Bloch's theorem. However, Bloch's theorem does not hold in cases like hyperbolic lattice connectivity. There is also little support for computing Hofstadter butterflies, while there exist multiple Python packages referring to band theory calculations. This brings us to what is interesting about the implementation of lattice simulations in this repository:

- We create base cases for the most common, physically realizable lattice types without reliance on crystal momentum to allow for extending this code to cases where Bloch's theorem no longer holds.

- We account for the presence of additional on-site disorder as well as a constant, perpendicular magnetic field. This allows us to simulate Anderson localization and plot Hofstader butterflies with and without disorder. The magnetic field support also allows for further invesitgation into the quantum hall effect.

- The classes creating the hamiltonians for each lattice type can be modified or refactored to support additional parameters or more accurately reflect a particular type of material such as graphene.

## Tight-Binding Hamiltonian
Suppose we take an atom. There are orbitals associated with it that describe regions where electrons could exist -- these are eigenfunctions of the Hamiltonian describing our atom. But what if we configure our atom with other atoms in a crystalline configuration? 
How do we describe the energy of a system where particles are placed in a crystalline configuration?
The following tight-binding Hamiltonian describes lattices in the absence of external interactions while accounting for hopping:

$$H = \omega_0 \sum_i a_i^{\dagger} a_i - t \sum_{\langle i,j\rangle} (a_i^{\dagger} a_j + a_j^{\dagger} a_i)$$

where $\omega_0$ is the on-site energy. The first sum runs over all lattice sites. The second sum describes hopping between nearest neighbors with a hopping amplitude $t$. It also encodes lattice geometry. 

When we're simulating specific lattice types, it's actually easier to construct a Hamiltonian in terms of a matrix populated by considering the lattice geometry. By this, we mean populating an $N$ x $N$ matrix based on the nearest neighbor hopping as it would occur on a particular lattice.

## Anderson Localization
The key idea behind Anderson Localization is that certain materials can undergo a phase transition from conductor to insulator if the system passes a disorder threshold. Thus, in systems of sufficiently large disorder (such as defected semiconductors) the electronic wavefunction associated with a now spatially localized state becomes localized. This localization influences the aforementioned phase transition. \textbf{In other words, spatial localization of the electronic wavefunction causes a change in the conductance of a highly disordered material.}

The Anderson tight-binding model allows us to describe this phenomenon more effectively. Here, electrons can tunnel between neighboring lattice sites until high disorder in the lattice causes quantum amplitudes associated with the tunneling paths to cancel amongst each other. A localized wavefunction follows. Equivalently, we can say that the incoming wave is scattered by potentials arising from the disorder, and the scattered wavelets destructively interfere as they propel forward at high disorder. This causes an exponential decay of the wavefunction. In this way, \textbf{Anderson localization} can be thought of as an interference phenomenon. 

Experimentally, electron localization has mostly been observed in a 1D case. 

The Anderson Hamiltonian can help us describe the localization in more technical terms. We write it as such

$$H = W \sum_n (\epsilon_n c_n^\dagger c_n) + t \sum_{<n,m>} (c_n^\dagger c_m + h.c)$$

where $t$ is the parameter describing the nearest hopping neighbor, $W$ is the disorder parameter, and $\epsilon_n$ is the random on-site energy in the range $[-\frac{1}{2},\frac{1}{2}]$.

## Inverse Participation Ratio (IPR)

The "participation ratio" gives an estimation of the localization length:
$$IPR = \frac{(\sum_x |\psi(x)|^2)^2}{ \sum_x |\psi(x)|^4}$$
(the numerator is not necessary if wavefunctions are normalized).

## Hofstadter Butterflies
What happens when we apply a perpendicular, uniform magnetic field onto a lattice? The general tight-binding hamiltonian will now involve a "Peierls phase" accounting for the magnetic flux through each plaquette as well as relevant changes in the boundary conditions.
An interesting result is that if we plot the energies as a function of magnetic flux ratios ($\phi = p/q$) such that $p$ and $q$ are coprime integers, we obtain a fractal pattern. It is a recursive structure. The way we constructed the butterfly involved choosing a maximum value for $q$, iterating through all the coprime $p$, $q$ pairs leading up to that point, and then reconstructing the hamiltonian for each consequent $\phi = p/q$. 

## Thouless Conductance

## Hall Conductance

## Chern Number

# Basic Usage

# Plot Examples

# Credits and Acknowledgements

Special thanks to Dr. Matteo Ippoliti for being my advisor and encouraging me to put everything on Github! A prior version of this repository is availible here as my project for PHY 329 (Computational Physics) at UT Austin ([class website](https://www.wgilpin.com/cphy/)). Parts of the code and README are the same because I wrote both of them. Consider this repository an extension of my old work. I'd like to acknowledge my friends Rida Siddiqi and Maddox Wroblewski for their assistance and input towards the initial repository (especially with organization and plotting logistics), as well as my friends Bishoy Kousa and Aravind Karthigeyan for their help and advice with troubleshooting.

I will compile a bibliography of all the sources I learned from and place it here when complete.
