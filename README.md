[![Documentation](https://readthedocs.org/projects/simos/badge/?version=latest)](https://simos.readthedocs.io/en/latest/?badge=latest)
[![CI](https://github.com/spinsimulation/simos/actions/workflows/main.yml/badge.svg)](https://github.com/spinsimulation/simos/actions/workflows/main.yml)
[![PyPI](https://img.shields.io/pypi/v/simos?logo=PyPI&color=0F81C1)](https://pypi.org/project/simos/)
[![Website](https://img.shields.io/badge/website-website)](https://simos.kherb.io)
[![Virtual Lab](https://img.shields.io/badge/Virtual%20Lab-8A2BE2)](https://simos.kherb.io/virtual/lab/index.html?path=Welcome.ipynb)

# SimOS Library

SimOS (**SIM**ulation of **O**ptically-adressable **S**pins) is a library for 
the simulation of open quantum systems consisting of spins, electronic levels 
and combinations thereoff. It can simulate the spin dynamics of conventional 
magnetic or electron paramagnetic resonance, but is further capable to simulate 
optically adressable spins or purely optical systems. ODMR combines electron 
spin resonance with optical measurements and is a highly sensitive technique 
to study systems possessing spin-dependent fluorescence. In our examples 
section, we feature two prototypical examples of optically adressable 
spins - the nitrogen vacancy center in diamond and photogenerated 
spin-correlated radical pairs.

Modelling the dynamics of open quantum system can be a non-trivial task. Their 
incoherent interaction with the system environment interferes with their coherent 
time evolution and an accurate simulation of their dynamics requires a quantum 
master equation. SimOS provides an interface for the construction of operators and 
super-operators of arbitrarily complex systems of spins and electronic levels and 
facilitates simulations of their dynamics on various levels of theory. Pythons 
popular numpy, scipy, qutip and sympy libraries are readily integrated and can be 
chosen as backends. We try to minimize high-level functions and keep the style of 
the simulations as close to a pen-and paper style as possible.

Our main focus lies on the QME in Lindblad form, for which we provide various engines 
for computationally efficient time propagation. In addition spatial dynamics such 
as rotational diffusion, linear flow or magic angle spinning can be simulated with 
our Fokker-Planck framework.


## How to Cite
If you find this project useful, then please cite:

Laura A. Völker, John M. Abendroth, Christian L. Degen and Konstantin Herb: <br>
*SimOS: A Python Framework for Simulations of Optically Addressable Spins* <br>
arXiv:2501.05922 <br>
https://doi.org/10.48550/arXiv.2501.05922


## Contribution guidelines
Documentation is important! 

Please use expressive and meaningful variable names. This makes it easier for 
all who are using the code and try to understand what's happening. In addition,
please comment you implementations. When you use constants, approximations or 
a special implementation scheme, please consider including a DOI of a
publication from which you retrieved your information.

In general, please follow the recommendations given in the Python Enhancement 
Proposals, especially [PEP-8](https://www.python.org/dev/peps/pep-0008/?) and 
[PEP-257](https://www.python.org/dev/peps/pep-0257/).

## Installation

### As end-user
The SimOS library is available at [PyPI](https://pypi.org/project/simos/). It can be 
installed via pip, the Python package installer. As the backend structure is modular, 
packages for the backends have to be installed separately. To install the recommended 
set of dependencies, run the following command:

```bash
pip install simos[recommended]
```

To install the minimal set of dependencies, run the following command:

```bash
pip install simos
```
Note, that if you install the minimal set of dependencies, executing the functions 
with the default backend method “qutip” will raise an error and you have to specify 
the backend method explicitly.


### As developer
To install in developer mode clone the repository and set your working directory to 
the root of the repository. Then install the package with the following command:
```bash
pip install numpy scipy sympy qutip pytest
pip install -e .
```
