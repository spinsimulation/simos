import numpy as _np
from .gyromagnetic_ratios import *
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
#: NV gyromagn. ratio (DOI: 10.1103/PhysRevB.79.075203)
yNV = ye * (1+357e-6) #[1/(sT)] 
#: Gyromagn. ratio 1H in H20 (CODATA 2018)
yH1_H20 = 2.675153151e8 #[1/(sT)]


#: Zero field splitting of NV (no reference, know-how of group)
D = 2.871e+9*2*_np.pi #[1/s]

#: Lattice constant Diamond at 300 K.
#: Source: http://7id.xray.aps.anl.gov/calculators/crystal_lattice_parameters.html
a = 3.56683e-10 #[m]

# Electric constant of diamond (DOI: 10.1063/1.3253121)
eps_r = 5.66 
golden_ratio = (1+_np.sqrt(5))/2

# Axial component of ground state permanent electric dipole moment of NV (DOI : 10.1016/0009-2614(90)85665-Y)
eps_z = 3.5e-3  # Hz/(V/m)
# Non-axial component of ground state permanent electric dipole moment of NV (DOI : 10.1016/0009-2614(90)85665-Y)
eps_xy = 170e-3  # Hz/(V/m)




