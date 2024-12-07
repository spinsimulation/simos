import numpy as _np
from .gyromagnetic_ratios import *

golden_ratio = (1+_np.sqrt(5))/2


###########################################################
# PHYSICAL CONSTANTS
###########################################################

#: Planck's constant (CODATA 2018)
hbar = 1.054571817e-34 #[Js]
h = 6.62607015e-34 #[J/Hz]

#: Speed of light (CODATA 2018)
c = 299792458  #[m/s]

#: Elementary charge (CODATA 2018)
elementary_charge = 1.602176634e-19 #[C]

#: electron mass (CODATA 2018)
me = 9.1093837015e-31 #[kg]

#: Bohr magneton (CODATA 2018)
mub = 927.40100783e-26 #[J/T]

#: g factor electron (CODATA 2018)
ge = 2.00231930436256 #[1]

#: Atomic mass u (CODATA 2018)
u = 1.66053906660e-27 #[kg]

#: Nuclear magneton (CODATA 2018)
munuc = 5.0507837461e-27 #[J/T]

#: Boltzman constant (CODATA 2018)
kB = 1.380649e-23 #[J/K]

#: Magic angle in radians
magic_angle = _np.arccos(1/_np.sqrt(3))

# MAS axis
mas_axis = _np.array([1/_np.sqrt(3), 1/_np.sqrt(3), 1/_np.sqrt(3)])

# Vacuum Permitivity (CODATA 2018)
eps_0 = 8.8541878128e-12 #[F/m]

# Vacuum Magnetic Permeability (CODATA 2018)
mu_0 = 1.25663706212e-6 #[N/A**2]

# About the sign of the gyromagn. ratio: 
# Electron has a negative gyromagnetic ratio, H1 positive
#  
# The nuclear magnetic moments are taken from the sources given in the IAEA table
# International Atomic Energy Agency, INDC(NDS)-0658, February 2014
# (https://www-nds.iaea.org/publications/indc/indc-nds-0658.pdf)

#: Electron gyromagn. ratio (CODATA 2018)
ye = -1.76085963023e11 #[1/(sT)]
#: Gyromagn. ratio 1H in H20 (CODATA 2018)
yH1_H20 = 2.675153151e8 #[1/(sT)]





