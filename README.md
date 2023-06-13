# Classical Simulation of Dipolar Atoms in a Trap
This program classically simulates N dipolar atomic atoms in a harmonic trapping potential using the Metropolis-Hastings Monte Carlo algorithm for sampling the Boltzmann distribution.

I implemented the algorithm in C as well as calculations for the number density, pair density and reblocking. I used Python and various scientific libraries for data analysis and the production of plots, examples of which are shown below.

![Low temperature number density contour plot for 15 atoms in a pancake-shaped trap.](images/density_contour.png)

![Plot of the pair correlation function at various different temperature for a pancake-shaped trap showing the phase transition of the system.](images/interparticle_distance.png)


## Dependencies
- `gsl`
- `matplotlib`
- `numpy`
- `scipy`
- `seaborn`
- `tomlc99`
- `tomli`
