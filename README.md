# Classical Simulation of Dipolar Atoms in a Trap
A Metropolis-Hastings classical Monte Carlo algorithm for sampling the Boltzmann distribution for dipolar atoms in a trapping potential.

## Outline
- `main`: File containing the main Metropolis-Hastings algorithm as well as some functions to process energies (reblocking) and to calculate error.
- `analysis`: Script that reads data produced from main to produce a number of plots.
- `config`: Input file to determine various properties of the system being simulated in main.

## Libraries Used
- C libraries:
    - [tomlc99](https://github.com/cktan/tomlc99)
    - [gsl](https://www.gnu.org/software/gsl/doc/html/index.html)
- Python libraries:
    - tomli
    - NumPy
    - matplotlib
