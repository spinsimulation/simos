import numpy as _np
import scipy as _sc
from ..constants import *
from ..trivial import *
from ..trivial import spher2cart, cart2spher
from ..propagation import evol, rot
from ..qmatrixmethods import applySuperoperator, ptrace, expect, tensor, data, tidyup, identity
from .. import backends
from ..core import System, subsystem, reverse_subsystem
from ..states import pol_spin
from itertools import combinations
from ..coherent import auto_zeeman_interaction, dipolar_coupling, dipolar_spatial
from ..incoherent import tidyup_ratedict, transition_operators
import copy



######################################################################
# Please note  - this is an NV center specific submodule. 
# Functions are not necessarily applicable to other spin/level systems.
######################################################################




###########################################################
# NV Specific Constants 
###########################################################

# NV gyromagnetic ratio (DOI: 10.1103/PhysRevB.79.075203).
yNV = ye * (1+357e-6) #[1/(sT)] 


# Zero field splitting fo the NV Center 
D = 2.871e+9*2*_np.pi #[1/s]; ZFS in the electronic ground state 
DESpara =  1.44e9*2*_np.pi # [1/s] ; electronic excited state parallel component 
DESperp =  1.541e9/2*2*_np.pi # [1/s] ; electronic excited state perpendicular component

# Spin-orbit interaction in the electronic excited state.
lambdaESpara = 5.33e9*2*_np.pi # [1/s]
lambdaESperp = 0.154e9*2*_np.pi # [1/s] 
gorb = 0.1 # a.u. , g value of the excited state orbit. 

# Permanent electric dipole moment of the NV. These are only known for gs, es values have
# to the best of our knowledge not been reported so far. 
eps_z  = 3.5e-3  # Hz/(V/m) ; ground state parallel component
eps_xz = 170e-3  # Hz/(V/m) ; ground state perpendicular component 

#: Lattice constant Diamond at 300 K.
a = 3.56683e-10 # [m] Lattice constant Diamond at 300 K. Source: http://7id.xray.aps.anl.gov/calculators/crystal_lattice_parameters.html
eps_r =  5.66  # Electric constant of diamond (DOI: 10.1063/1.3253121)


#########################
# NV  System Child Class  
#########################

class NVSystem(System):
    """ A class for initialising quantum systems of NV centers.
    """
    # Constructor. 
    def __init__(self, optics = True, orbital = False, nitrogen= True, natural = False, further_spins = [], method = 'qutip'):
        # NV center electron spin.
        S = {'val': 1, 'name':'S', 'type': 'NV-'} 
        # Electronic states.
        GS = {'val': 0 , 'name':'GS'}
        if orbital == False:
            ES = {'val': 0 , 'name':'ES'}
        else:
            ES = {'val': 1/2 , 'name':'ES'}
        SS = {'val': 0 , 'name':'SS'}
        if optics == True:
            system_array = ([(GS, ES), S], SS)
        else:
            system_array = [S]
        # Further Spins.
        if nitrogen == True:
            if natural == True:
                N = {'val': 1, 'name':'N', 'type': 'N'}
            else:
                N = {'val': 1/2, 'name':'N', 'type': 'N'}            
            further_spins = [N] + further_spins
        if len(further_spins) > 0:
            names = [i["name"] for i in further_spins]
            if "S" in names:
                raise ValueError("Name S is reserved for the NV Center Electron Spin.")
            system_array = [system_array] + further_spins
        self.photooptions = {"optics": optics, "orbital": orbital}
        self.spinoptions =  {"nitrogen": nitrogen, "furtherspins": len(further_spins) > 0}
        super().__init__(system_array, method = method)

    # Function to get Hamiltonian that includes ZFS 
    # as well as interactions with magnetic and electric fields. 
    def field_hamiltonian(self, Bvec = _np.array([0,0,0]), EGS_vec = _np.array([0,0,0]), EES_vec = _np.array([0,0,0])):
        """ Hamiltonian of the NV center electron spin in the presence of magnetic and electric/strain fields. 

            :Keyword Arguments:
            * *B_vec*  --  Cartesian magnetic field vector, in Tesla.
            * *EGS_vec* -- Cartesian vector representing products of electrin/strain fields and NV electric dipole moment in the GS, in 1/s.
            * *EES_vec* -- Cartesian vector representing products of electrin/strain fields and NV electric dipole moment in the ES, in 1/s.
            
            :returns: Field Hamiltonian of the NV center.
        """
        # No differentiation btw. GS, ES in this case. Just call the NV 
        if self.photooptions["optics"] is False:
            HGS =  D * ((self.Sz)**2 - 2/3*self.Sid) # Zero field splitting
            HGS +=  -yNV*(self.Sx*Bvec[0] + self.Sy*Bvec[1] + self.Sz*Bvec[2])  # Zeeman interaction 
            HGS +=  EGS_vec[0]*(self.Sy**2 - self.Sx**2) + EGS_vec[1]*(self.Sx*self.Sy - self.Sy-self.Sx) + EGS_vec[2]*((self.Sz)**2 - 2/3*self.Sid) # Electric field interaction
            return HGS
        else:
            # Ground state Hamiltonian 
            HGS =  D * ((self.Sz*self.GSid)**2 - 2/3*self.Sid*self.GSid)  # ZFS GS
            HGS +=  -yNV * (self.Sx*self.GSid*Bvec[0] + self.Sy*self.GSid*Bvec[1] + self.Sz*self.GSid*Bvec[2])  # Zeeman GS
            HGS +=  EGS_vec[0]*self.GSid*(self.Sy**2 - self.Sx**2)  + EGS_vec[1]*self.GSid*(self.Sx*self.Sy - self.Sy*self.Sx) + EGS_vec[2]*self.GSid*((self.Sz)**2 - 2/3*self.Sid) # Electric field interaction
            # Room temperature (no orbital substructure in ES)
            if self.photooptions["orbital"] is False:
                # Excited state Hamiltonian
                HES =  DESpara * ((self.Sz*self.ESid)**2 - 2/3*self.Sid*self.ESid)  # ZFS ES 
                HES +=   -yNV * (self.Sx*self.ESid*Bvec[0] + self.Sy*self.ESid*Bvec[1] + self.Sz*self.ESid*Bvec[2]) # Zeeman 
                HES +=   EES_vec[2]*2*self.ESid*self.Sid + _np.sqrt(EES_vec[0]**2 + _np.sqrt(EES_vec[1]**2))*(self.ESid*(self.Sy**2 - self.Sx**2))
                return HGS, HES    
            # Variable Temperature (+ orbital in ES)         
            else:
                # Excited state Hamiltonian
                HES =  DESpara * ((self.Sz*self.ESid)**2 - 2/3*self.Sid*self.ESid) + DESperp*(2*self.ESz*(self.Sy**2 - self.Sx**2) - 2*self.ESx*(self.Sy*self.Sx + self.Sx*self.Sy))   # ZFS
                HES += - lambdaESpara*2*self.ESy*self.Sz  + lambdaESperp*(2*self.ESz*(self.Sx*self.Sz + self.Sz*self.Sx) - 2*self.ESx*(self.Sy*self.Sz + self.Sz*self.Sy))  # Spin-Orbit 
                HES +=   -yNV * (self.Sx*self.ESid*Bvec[0] + self.Sy*self.ESid*Bvec[1] + self.Sz*self.ESid*Bvec[2]) # Zeeman 
                HES +=   ((mub*gorb)/hbar)*Bvec[2]*2*self.ESy*self.Sid # Zeeman orbital 
                HES +=   EES_vec[0]*2*self.ESz*self.Sid  +  EES_vec[1]*2*self.ESx*self.Sid  +  EES_vec[2]*self.ESid*self.Sid # Strain / Electric Field 
                return HGS, HES

    def transition_operators(self, T=298, beta = 0.2, Bvec = _np.array([0,0,0]), Evec = _np.array([0,0,0])):
        """ Collapse operators for laser- and phonon-induced incoherent transitions between electronic levels of the NV center.
        For a classical rate model of the NV center, only beta must be specified and all other parameters can be omitted. 

        :Keyword Arguments:
            * *T*  --  Temperature (in Kelvin), defaults to 298 K.
            * *beta* -- A parameter from 0-1 to indicate laser power relative to saturation power (saturation for beta = 1).
            * *Bvec* -- Cartesian magnetic field vector, in Tesla.
            * *Evec* -- Cartesian vector representing products of electrin/strain fields and NV electric dipole moment in the ES, in Hz.

        :returns: Two lists of collapse operators, first list does include laser excitation, second list does not include laser excitation. 
        """
        # Test whether NV was build with electronic levels.
        # Otherwise no collapse operators can be returned. 
        if self.photooptions["optics"] is False:
            raise ValueError("Rate dictionary requires included electronic levels.")
        rates = decay_rates(T)
        NV_rates = {}
        laser_rates = {}
        # Classical Rate Model if orbital branches are not being considered. 
        if self.photooptions["orbital"] == False:
            if T < 200:
                raise Warning("Classical Model may be inaccurate for low temperatures. Consider using a full quantum mechanical model.")
            NV_rates = {}
            #  GS <-> ES transitions
            NV_rates["ES,S[1]->GS,S[1]"]  = rates["optical_emission"]
            NV_rates["ES,S[-1]->GS,S[-1]"] = rates["optical_emission"]
            NV_rates["ES,S[0]->GS,S[0]"] = rates["optical_emission"]
            # Spin dependent rates for the ISC (Intersystem crossing) process
            NV_rates["ES,S[1]->SS"] = rates["ISC_ms1"]
            NV_rates["ES,S[-1]->SS"] =rates["ISC_ms1"]
            NV_rates["ES,S[0]->SS"] = rates["ISC_ms0"]
            # Decay from shelving to ground state
            NV_rates["SS->GS,S[0]"] = rates["ssdecay_ms0"]
            NV_rates["SS->GS,S[1]"] = rates["ssdecay_ms1"]
            NV_rates["SS->GS,S[-1]"] = rates["ssdecay_ms1"]
            # Tidyup the rate dictionary
            NV_rates = tidyup_ratedict(self, NV_rates)
            # Laser excitation
            laser_rates["GS,S[0]->ES,S[0]"] =  beta*rates["optical_emission"]
            laser_rates["GS,S[1]->ES,S[1]"] =  beta*rates["optical_emission"]
            laser_rates["GS,S[-1]->ES,S[-1]"] = beta*rates["optical_emission"]
            laser_rates = tidyup_ratedict(self, laser_rates)
            return transition_operators(self, fuse_dictionaries(NV_rates, laser_rates)), transition_operators(self, NV_rates)
        # QM Rate Model 
        elif self.photooptions["orbital"] == True:
            # Calculate the perpendicular strain.
            E_perp = _np.sqrt(Evec[0]**2 + Evec[1]**2)
            # Get phonon rates.  
            kup, kdown = phonon_rates(T, E_perp)
            # If the system contains additional spins, project onto the subspace.
            if self.spinoptions["nitrogen"] or self.spinoptions["furtherspins"]:
                _ , rop , reverse = subsystem(self, self.id, ["S", "GS", "ES", "SS"])
                smallNV = NVSystem(optics = True, orbital = True, nitrogen = False, method = self.method)
            else:
                smallNV  = self
            # Set up the alternate bases.  
            # ZF basis (independent of Hamiltonian).   
            T_zftoez =  _np.array([[ 1,    0,    0,     0,    0,    0,    0,    0,    0 , 0 ],
                                [ 0,    1,    0,     0,    0,    0,    0,    0,    0 , 0 ],
                                [ 0,    0,    1,     0,    0,    0,    0,    0,    0 , 0 ],
                                [ 0,    0,    0, -1j/2,  1/2,    0,    0,-1j/2, -1/2  , 0],
                                [ 0,    0,    0,     0,    0,    0,    1,    0,    0 ,0 ],
                                [ 0,    0,    0, -1j/2, -1/2,    0,    0,-1j/2,  1/2 ,0 ],
                                [ 0,    0,    0,   1/2, 1j/2,    0,    0, -1/2, 1j/2  ,0],
                                [ 0,    0,    0,     0,    0,    1,    0,    0,    0 ,0 ],
                                [ 0,    0,    0,  -1/2, 1j/2,    0,    0,  1/2, 1j/2 ,0 ],
                                [ 0,    0,    0,     0,    0,    0,    0,    0,    0 , 1 ]])
            T_eztozf = _np.transpose(T_zftoez)
            smallNV.add_basis(T_eztozf, "zf", ["GS1", "GS0", "GSm1", "E1", "E2", "Ey0", "Ex0", "A1", "A2", "SS"])
            # HF basis (dependent on ES Hamiltonian)
            HGS, HES = smallNV.field_hamiltonian(Bvec, Evec)
            HES_sub_orbit, HES_remaining, reverse_HES_sub_orbit = subsystem(smallNV, HES, "ES")
            HES_remainingu  = [identity(r.shape[0], dims = r.dims) for r in HES_remaining]
            vals, vecs = HES_sub_orbit.eigenstates()
            T_sub = _np.zeros(HES_sub_orbit.shape, dtype = complex) 
            for ind_s, s in enumerate(vecs):
                T_sub[ind_s,:] = data(s)[:,0]
            T_sub = tidyup(T_sub, method = self.method, dims = HES_sub_orbit.dims)
            T_hftoez = data(reverse_subsystem(T_sub, HES_remainingu, reverse_HES_sub_orbit))
            T_eztohf = _np.transpose(T_hftoez) 
            smallNV.add_basis(T_eztohf, "hf", ["GS1", "GS0", "GSm1", "Ex1", "Ex0", "Exm1", "Ey1", "Ey0", "Eym1", "SS"])              
            # Define the rates.
            NV_rates = {}
            laser_rates = {}
            # Laser excitation, here we assume symmetric excitation of both orbital branches. 
            #  GS <-> ES transitions (laser excitation and fluorescent decay)
            laser_rates["hf_GS0 -> hf_Ex0"] =  0.5*beta*rates["optical_emission"]
            laser_rates["hf_GS0 -> hf_Ey0"] =  0.5*beta*rates["optical_emission"]
            laser_rates["hf_GS1 -> hf_Ex1"] = 0.5*beta*rates["optical_emission"]
            laser_rates["hf_GS1 -> hf_Ey1"] = 0.5*beta*rates["optical_emission"]
            laser_rates["hf_GSm1 -> hf_Exm1"]  = 0.5*beta*rates["optical_emission"]
            laser_rates["hf_GSm1 -> hf_Eym1"]  = 0.5*beta*rates["optical_emission"]
            # All decays: 
            NV_rates["GS,S[0] <- ES[-0.5],S[0]"] =  rates["optical_emission"]
            NV_rates["GS,S[0] <- ES[0.5],S[0]"] =  rates["optical_emission"]
            NV_rates["GS,S[1] <- ES[-0.5],S[1]"] = rates["optical_emission"]
            NV_rates["GS,S[1] <- ES[0.5],S[1]"] =   rates["optical_emission"]
            NV_rates["GS,S[-1] <- ES[-0.5],S[-1]"]  =rates["optical_emission"]
            NV_rates["GS,S[-1] <- ES[0.5],S[-1]"]  = rates["optical_emission"]
            # Spin dependent rates for the ISC (Intersystem crossing) process
            NV_rates["zf_A1->zf_SS"] = rates["ISC_ms1"]/0.52
            NV_rates["zf_E1->zf_SS"] = rates["ISC_ms1"]
            NV_rates["zf_E2->zf_SS"] = rates["ISC_ms1"]
            NV_rates["zf_Ex0->zf_SS"] = rates["ISC_ms0"]
            NV_rates["zf_Ey0->zf_SS"] = rates["ISC_ms0"]
            # Spin dependent decay from shelving to ground state
            NV_rates["SS->GS,S[0]"] = rates["ssdecay_ms0"]
            NV_rates["SS->GS,S[1]"] = rates["ssdecay_ms1"]
            NV_rates["SS->GS,S[-1]"] = rates["ssdecay_ms1"]
            # Phonon induced transition rates
            NV_rates["hf_Ey1 -> hf_Ex1"] = kup
            NV_rates["hf_Ey0 -> hf_Ex0"] =  kup
            NV_rates["hf_Eym1 -> hf_Exm1"] = kup
            NV_rates["hf_Ey1 <- hf_Ex1"] = kdown
            NV_rates["hf_Ey0 <- hf_Ex0"] =  kdown 
            NV_rates["hf_Eym1 <- hf_Exm1"] = kdown

            NV_rates = tidyup_ratedict(smallNV, NV_rates)
            laser_rates =  tidyup_ratedict(smallNV, laser_rates)
            all_rates = fuse_dictionaries(laser_rates, NV_rates)
            all_rates = tidyup_ratedict(smallNV, all_rates)
            c_ops_on = transition_operators(smallNV, all_rates)
            c_ops_off = transition_operators(smallNV, NV_rates)
            if self.spinoptions["nitrogen"] or self.spinoptions["furtherspins"]:
                c_ops_on_real = []
                c_ops_off_real = []
                for c_op in c_ops_on:
                    ropu  = [r.unit()*r.shape[0] for r in rop]
                    c_ops_on_real.append(reverse_subsystem(c_op, ropu, reverse))
                for c_op in c_ops_off:
                    c_ops_off_real.append(reverse_subsystem(c_op, ropu, reverse))
                return c_ops_on_real, c_ops_off_real
            else:
                return c_ops_on, c_ops_off


#####################
# NV Transition rates
#####################

# Optical decay rates. 
def decay_rates(T = 298):
    """ Optical decay rates for the NV center at a given temperature. See
    DOI: 10.1103/PhysRevB.108.085203 for in-depth explanation.

    :param T: Temperature in Kelvin.
    :returns: A dictionary with optical decay rates for fluorescent emission, spin-selective inter-system crossing and decay from the shelving state to the ground state.

    """
    rates = {}
    # Optical emission from the ES to the GS.
    rates["optical_emission"]  = 55.7e6 # [Hz]
    # Average intersystem crossing rate (ES to SS) for ms pm 1.
    rates["ISC_ms1"] = 98.7e6 # [Hz]
    # Average intersystem crossing rate (ES to SS) for ms 0.
    rates["ISC_ms0"] = 8.2e6 # [Hz]
    # Spin-selective decay from SS to GS.
    deltaE = 16.6e-3*elementary_charge #  [J] SS emitted phonon energy 
    rs = 2.26 # SS branching ratio 
    tau = 320e-9 # s # SS decay time t T = 0K 
    if T > 0 :
        tau = tau * (1-_np.exp(-deltaE/(kB*T)))
    s = 1/tau
    s0 = rs*s/(1+rs)
    s1 = s-s0
    rates["ssdecay_ms0"] = s0
    rates["ssdecay_ms1"] = s1
    return rates 

# Rates for Phonon-Mediated Hopping.
def phonon_rates(T, E_perp):
    """ Rates for phonon-mediated hopping between the orbital branches of the NV center excited state. See
    DOI: 10.1103/PhysRevB.108.085203 for in-depth explanation.

    :param T: Temperature in Kelvin.
    :param E_perp: excited state in-plane strain/el. field in Hz. 

    :returns: Up- and down rates for phonon-mediated hopping in Hz.

    """

    def DebyeIntegrand(x, T, Delta):
        x_perp = hmod*Delta/(kBmod*T)
        if x <= x_perp:
            quotient  = _np.float64(0) 
        else:
            dividend = 0.5 * _np.exp(-x) * x * (x-x_perp) * (x**2 + (x-x_perp)**2)
            divisor = _np.exp(-x_perp) - _np.exp(-x) - _np.exp(-x_perp-x) + _np.exp(-2*x)
            if divisor == 0.:
                quotient = _np.float64(0.)
            else:
                quotient = dividend/divisor
        return quotient
    # Define required constants that are not globally available.
    # Convert units of some globally available constants.
    eta = 176*1e15 # s/(eV)**3; ES electron phonon coupling strength
    kBmod = kB/elementary_charge  # Boltzman constant, unit: eV/K
    hmod = h/elementary_charge # Plank constant to convert energy in GHz to eV. Units: eV/Hz
    hbarmod = hmod/(2*_np.pi)
    Delta = 2*E_perp
    # Calculate the one phonon downwards rate. 
    if Delta == 0:
        kdown1 = 0  # spectral density goes to 0 
    elif T > 0 :
        kdown1 = 32*eta*hmod**3*E_perp**3*(1/(_np.exp(2*E_perp*hmod/(kBmod*T))-1) + 1)
    else:
        kdown1 = 32*eta*hmod**3*E_perp**3
    if kdown1 < 1e0:
        kdown1 = 0 
    # Calculate the two phonon downwards rate.
    if T == 0:
        kdown2 = 0
    else:
        phononcutoffenergy = 168e-3 # eV
        cutoff = phononcutoffenergy/(kBmod*T)
        x_perp = 2*hmod*E_perp/(kBmod*T)
        Integral = _sc.integrate.quad(DebyeIntegrand, x_perp, cutoff, args = (T, Delta))[0]
        kdown2 = (64/_np.pi) * hbarmod * eta**2*kBmod**5*T**5  * Integral   # fac4? 
        khopplimit = 250e12 # Hz 
        kdown2 = _np.min([kdown2, khopplimit])
        if kdown2 < 1e0:
            kdown2 = 0 
    # Calculate total downwards and upwards rates from here. 
    kdown = kdown1+ kdown2
    if T > 0:
        kup = kdown * _np.exp(-2*hmod*E_perp/(kBmod*T))
    else:
        kup = 0 
    return kup, kdown 


###########################################################
# INTERACTIONS with nuclear spins 
###########################################################


def coupl2geom(apara,aperp,angular=False,ynuc=yC13,yel=ye):
	"""Converts the parallel and perpendicular coupling frequencies as obtained from correlation spectroscopy to r and theta assuming a point dipole model.
	
	Returns [r, theta]
	
	Parameters:
	apara : parallel coupling frequency
	aperp : perpendicular coupling frequency

	Optional parameters:
	angular : by default False, assuming frequencies in Hz as input, set to True for inputting angular frequencies
	ynuc    : gyromagnetic ratio for the nuclear spin, by default built-in value for C13
	yel     : gyromagnetic ratio for the electron spin,  by default built-in value for electron
	"""
    # Convert to angular
	y_e = _np.abs(yel)
	y_n = _np.abs(ynuc)
	if not angular:
		apara = 2*_np.pi*apara
		aperp = 2*_np.pi*aperp
    
	# Formulas taken from DOI 10.1103/PhysRevLett.116.197601
	theta = 1/2*(-3*apara/aperp+_np.sqrt(9*apara**2/aperp**2 + 8))
	theta = _np.arctan(theta)
	r = 1e-7 * y_e*y_n*hbar*(3*_np.cos(theta)**2-1)/(apara)
	r = r**(1/3)
    
	return [r, theta]

def auto_pairwise_coupling(spin_system,approx=False,only_to_NV=False,output=True):
    """ Returns the combined Hamiltonian for all pairwise dipolar couplings of a spin system using the tabulated isotropic gyromagnetic
    ratios of their type and positions.

    This routine requires that the dictionaries of the spin members specify their spin type and position.

    :param System spin_system: An instance of the System class.
    
    :Keyword Arguments:
    * *approx* (``bool``) -- If True, a secular approximation is performed. If False, the full coupling is considered.
    * *only_to_NV* (``bool``) -- If True, only couplings to NV center electronic spins are considered.
    * *output* (``bool``) -- If True, output is printed.

    """
    if approx:
        approx = 'secular'
    else:
        approx = 'Full'
    NV_name = 'S'
    NV_present = False
    Nitrogen_name = 'N'
    Nitrogen_present = False
    for spin in spin_system.system:
        spin_type = spin['type'] 
        if spin_type in ['NV', 'NV-', 'NV+']:
            NV_name =  spin['name']
            NV_present = True
        elif spin_type == 'NV-Nitrogen':
            Nitrogen_present = True
            Nitrogen_name = spin['name']

    H = spin_system.id * 0
    
    # Couplings to the NV center
    if NV_present:
        for spin in spin_system.system:
            spin_name = spin['name']
            spin_type = spin['type']
            #print(spin)
            if (spin_type not in ['NV', 'NV-', 'NV+', 'NV-Nitrogen']):
                # This is not the NV!
                #print('Coupling...')
                spin_op =  getattr(spin_system, spin_name + '_vec')
                if spin_type == 'blind':
                    continue
                else:
                    distance, theta, phi = spin['pos']
                if spin_type in ['NV', 'NV-', 'NV+','electron']:
                    y = ye
                elif spin_type == "15N":
                    y = yN15
                elif spin_type == "1H":
                    y = yH1
                elif spin_type == "13C":
                    y = yC13
                elif spin_type == "19F":
                    y = yF19
                elif spin_type == "blind":
                    y = 0
                else:
                    raise ValueError("Unknown spin type " + str(spin_type))
                #H += dipolar(distance,theta,phi,ye,y,NV_op,spin_op,approx=approx)
                H += dipolar_coupling(spin_system, NV_name, spin_name, ye, y, distance,theta,phi, mode = 'spher', approx = approx)
                if output:
                    print("---------------")
                    print(spin_name, "to NV")
                    mat = dipolar_spatial(ye,y,distance,theta,phi)
                    #print(mat)
                    apara = mat[2,2]
                    aperp = _np.sqrt(mat[2,0]**2 + mat[2,1]**2)
                    print('apara =', _np.round(_np.abs(w2f(apara))*1e-3,3), 'kHz')
                    print('aperp =',_np.round(_np.abs(w2f(aperp))*1e-3,3), 'kHz')
    
    # Nitrogen
    if Nitrogen_present:
        # Sometimes weired offset of 1e-6 in apara
        apara = f2w(3.03e6) #*2*np.pi
        aperp = f2w(3.65e6) # *2*np.pi
        H += apara * getattr(spin_system, NV_name + 'z') * getattr(spin_system, Nitrogen_name + 'z') 
        H += aperp * (getattr(spin_system, NV_name + 'x') * getattr(spin_system, Nitrogen_name + 'x') + getattr(spin_system, NV_name + 'y') * getattr(spin_system, Nitrogen_name + 'y'))
    
    if only_to_NV:
        return H

    for spin1, spin2 in combinations(spin_system.system,2): #pairwise(spin_system.system):
        spin_name1 = spin1['name']
        spin_type1 = spin1['type']
        spin_name2 = spin2['name']
        spin_type2 = spin2['type']
        if (spin_name1 != spin_name2) and not (spin_type1 in ['NV', 'NV-', 'NV+', 'NV-Nitrogen'])  and not (spin_type2 in ['NV', 'NV-', 'NV+', 'NV-Nitrogen']):
            pos1_cart = spher2cart(spin1['pos'])
            pos2_cart = spher2cart(spin2['pos'])
            delta = _np.array(pos2_cart) - _np.array(pos1_cart)
            delta_spherical = cart2spher(delta)

            spin_op1 =  getattr(spin_system, spin_name1 + '_vec')
            if spin_type1 == "15N":
                y1 = yN15
            elif spin_type1 == "1H":
                y1 = yH1
            elif spin_type1 == "13C":
                y1 = yC13
            elif spin_type1 == "19F":
                y1 = yF19
            elif spin_type1 == "blind":
                y1 = 0
            elif spin_type1 == "electron":
                y1 = ye
            else:
                raise ValueError("Unknown spin type")

            spin_op2 =  getattr(spin_system, spin_name2 + '_vec')
            if spin_type2 == "15N":
                y2 = yN15
            elif spin_type2 == "1H":
                y2 = yH1
            elif spin_type2 == "13C":
                y2 = yC13
            elif spin_type2 == "19F":
                y2 = yF19
            elif spin_type2 == "blind":
                y2 = 0
            elif spin_type1 == "electron":
                y2 = ye
            else:
                raise ValueError("Unknown spin type")

            #H += dipolar(delta_spherical[0],delta_spherical[1],delta_spherical[2],y1,y2,spin_op2,spin_op1,approx=approx_nucnuc)
            #H += dipolar(distance,theta,phi,ye,y,NV_op,spin_op,approx=approx)
            #H += dipolar_coupling(spin_system, NV_name, spin_name, ye, y, distance,theta,phi, mode = 'spher', approx = 'Full')
            H += dipolar_coupling(spin_system, NV_name, spin_name, y1, y2, delta_spherical[0],delta_spherical[1],delta_spherical[2], mode = 'spher', approx = 'Full')

            if output:
                print("---------------")
                print(spin_name1, " with ", spin_name2)
                print("\u0394r =", _np.round(delta_spherical[0]*1e10,2), "A, \u0394\u03b8 =", _np.round(_np.rad2deg(delta_spherical[1]),2),"°,  \u0394\u03D5 =", _np.round(_np.rad2deg(delta_spherical[2]),2),'°')
                print(w2f(dipolar_spatial(y1,y2,delta_spherical[0],delta_spherical[1],delta_spherical[2])))        

    return H


####################################################
# State initialization
####################################################

def gen_rho0(spinsystem, NV_name='S'):
    rho_arr = []
    method = spinsystem.method
    for i in range(len(spinsystem.system)):
        spin = spinsystem.system[i]
        spin_type = spin['type']
        if spin_type in ['NV','NV-','NV+']:
            op0 = getattr(spinsystem, NV_name + 'p')[0]
            op0 = op0.ptrace(i)/2
            rho_arr.append(op0.unit())
        else:
            if 'pol' in spinsystem.system[i].keys():
                rho_arr.append(pol_spin(spinsystem.system[i]['pol'],method=method).unit())
            else:
                 rho_arr.append(pol_spin(1,method=method).unit())
    rho0 = getattr(getattr(backends,method), 'tensor')(rho_arr).unit()
    return getattr(getattr(backends,method), 'tidyup')(rho0)


####################################################
## NV axes , coordinate transformations
####################################################

def get_NVaxes(orientation = (0,0,1), axisind = (1,1,1)):
    """Calculates the NV axis in the laboratory frame (z axis orthogonal to diamond surface) for the specified diamond surface termination and NV axis.
    Input: 
    - Orietntation: Diamond surface ermination, most commonly (0,0,1)
    - Axisind: NV axis specification, e.g. (1,1,1)
    Outout:
    - np.array specifying the NV axis"""

    NV  = _np.array(axisind)     # NV axis in diamond frame of reference 
    NV = NV/_np.linalg.norm(NV)  # Normalize 
    miller = _np.array(orientation) # Surface normal for desired miller plane in diamond frame 
    miller = miller / _np.linalg.norm(miller) # Normalize 
    # Get angles (theta) between the axis and the surface normal, 
    if NV[0] == miller[0] and NV[1] == miller[1] and NV[2] == miller[2]:
        theta = 0
    else:
        theta = _np.arccos(_np.dot(NV, miller) / (_np.linalg.norm(NV) * _np.linalg.norm(miller)))
    phi = 0    # Phi always zero bc we always want NV to lie in xz plane  
    # Calculate NV coordinate system 
    xaxis = _np.array([_np.cos(theta)*_np.cos(phi),_np.cos(theta)*_np.sin(phi), -_np.sin(theta)]) 
    yaxis = _np.array([_np.sin(phi), _np.cos(phi), 0]) 
    zaxis = _np.array([_np.sin(theta)*_np.cos(phi), _np.sin(theta)*_np.sin(phi),_np.cos(theta)]) 

    return _np.array([xaxis, yaxis, zaxis])

def lab2NVcoordinates(*args, orientation = (0,0,1), axisind = (1,1,1),  input_coords='cart', output_coords='cart'):
    """ Transforms vector from NV coordinate system to laboratory system, arbitrary diamond surface termination and any NV axis"""
    # Input handling
    if len(args) == 3:
        vector_lab = _np.array([args[0], args[1], args[2]]).transpose()
    elif len(args) == 1:
        vector_lab = args[0]
    else:
        raise ValueError("Wrong number of coordinates provided")
    # Pack
    if len(vector_lab.shape) == 1:
        vector_lab = _np.array([vector_lab])
    # Spher2cart
    if input_coords == 'spher':
        vector_lab = spher2cart(vector_lab)
    # Get NV axes  
    NV_axes = get_NVaxes(orientation, axisind)
    # Transform
    vector_NV = _np.dot(vector_lab, NV_axes)
    if output_coords == 'spher':
        vector_NV = cart2spher(vector_NV)
    # De-Pack
    if _np.shape(vector_NV)[0] == 1: # de-pack 
            return vector_NV[0]
    else:
        return vector_NV


###########################################################
# DYNAMICAL DECOUPLING SEQUENCES
###########################################################

def XY8(H, tau, spinsystem, *rho, N=8, phi=0, NV_name='S', c_ops=[], dthet=0):
    """ Simulation of XY8 sequence
     
    .x..y..x..y..y..x..y..x. (. denotes tau/2, x and y are the pi pulse axes)
    
    Args:
      H          : Hamiltonian for free evolution.
      tau        : interpulse delay
      spinsystem : Spinsystem to take the rotation operators out of.

    Other Parameters:
      *rho    : Density matrix or state to which the sequence is applied. If none is given, the Propagator will be returned (only for direct method, see below)
      N       : Number of pulses (must be multiple of 8)
      phi     : change rotation axis, by default 0.
      NV_name : lookup prefix for the rotation operators in spinsystem (will look up '[NV_name]x_red' and corresp. y)
      method  : 'direct': use evol function and assume no relaxation or timedependence of the Hamiltonian
      method  : 'mesolve': Pass the collaps operators c_ops to the Master equation solver of Qutip and use this solver for the timeevolution between the pi pulses (slow!).
      args    : Arguments for the mesolve routine
      options : Options for the mesolve routine
      dthet   : Pulse error in z rotation
    
    References:
      https://doi.org/10.1016/0022-2364(90)90331-3
    """
    # For robustness see DOI 10.1103/PhysRevLett.106.240501

    S_x_lookup = getattr(spinsystem, NV_name + 'op_x_red')
    S_y_lookup = getattr(spinsystem, NV_name + 'op_y_red')

    # If phi != 0 is specified, rotate the reference frame and use "effective" x & y operators
    S_x = _np.cos(phi)*S_x_lookup + _np.sin(phi)*S_y_lookup
    S_y = _np.cos(phi)*S_y_lookup + _np.sin(phi)*S_x_lookup
    
    U_xrot = rot(S_x,_np.pi+dthet)
    U_yrot = rot(S_y,_np.pi+dthet)
    U_t = evol(H,tau)
    U_t2 = evol(H,tau/2)

    method = backends.get_backend(H)
    isketfun = getattr(getattr(backends,method), 'isket')
    
    # time order:  .x..y..x..y..y..x..y..x.
    # order inverse to time order

    if len(c_ops) == 0:
        UCP = U_t2*U_xrot*U_t*U_yrot*U_t*U_xrot*U_t*U_yrot*U_t*U_yrot*U_t*U_xrot*U_t*U_yrot*U_t*U_xrot*U_t2
        UCP = UCP**int(N/8)
        if len(rho) == 0:
            return UCP
        elif len(rho) == 1:
            if isketfun(*rho):
                return UCP*rho[0]
            else:
                return UCP*rho[0]*UCP.dag()
        else:
            raise ValueError('More than one rho was provided. This is not supported.')
    else:
        # Propagate with Lindbladian
        import qutip as qu
        U_xrot_super = qu.to_super(U_xrot)
        U_yrot_super = qu.to_super(U_yrot)
        L = qu.liouvillian_ref(H,c_ops)
        U_t2_super = (L*tau/2).expm()
        U_t_super = U_t2_super*U_t2_super
        UCP_super = U_t2_super*U_xrot_super*U_t_super*U_yrot_super*U_t_super*U_xrot_super*U_t_super*U_yrot_super* \
                    U_t_super*U_yrot_super*U_t_super*U_xrot_super*U_t_super*U_yrot_super*U_t_super*U_xrot_super*U_t2_super
        UCP_super = UCP_super**int(N/8)
        if len(rho) == 0:
            return UCP_super
        elif len(rho) == 1:
            rho_vec = qu.operator_to_vector(rho[0])
            rho_vec = UCP_super*rho_vec
            return qu.vector_to_operator(rho_vec)
        else:
            raise ValueError('More than one rho was provided. This is not supported.')


###########################################################
# Helper functions for NV Measurement and Plotting
###########################################################

def meas_NV(rho, spinsystem, NV_name='S'):
    """Perform measurement operation on NV.
    
    :param rho: Density matrix of the state to measure.
    :param spinsystem: Spinsystem class object containing all operators
    :param NV_name: NV center spin name. If None, search for all NVs

    :returns: expectation value, rho after measurement.
    """
    op0 = getattr(spinsystem, NV_name + 'p')[0]
    val0 = expect(op0,rho)
    #nucs, reverse = subsystem(spinsystem, rho, NV_name,keep=False)
    #op0, _ = subsystem(spinsystem, op0, NV_name, keep=True)
    #rho = reverse_subsystem(spinsystem, [op0,nucs],reverse)
    # temporary fix
    if len(rho.dims[0]) > 1:
        sel = _np.arange(len(rho.dims[0])-1)+1
        rho_elec = ptrace(op0,0)
        rho_nucs = ptrace(rho,sel)
        rho = tensor([rho_elec,rho_nucs])
    rho = rho.unit()
    return _np.real(val0), rho

def exp2cts(x,contrast,ref_cts,ref=0.5):
    """Calculates photoluminescence counts for a given spin state initialization fidelity. 
    
    :param x: Expectation value for ms(0) sublevel.
    :param contrast: ODMR contrast of the NV.
    :param ref_cts: Photoluminescence counts of the NV center (cts/s) in the ms(0) state.

    :returns: Photoluminescence counts (cts/s)
    """
    m = contrast*ref_cts/(1-contrast+contrast*ref)
    c = -ref_cts*(-1 + contrast)/(1-contrast+contrast*ref)
    return m*x+c

def normalize_data(data,upper_ref,lower_ref):
    """Normalizes  data such that the upper reference correponds to 1, the lower reference trace to 0.

    :param data: Data trace to normalize
    :param upper_ref: mean value of the upper reference
    :param lower_ref: mean value of the lower reference

    :returns: Normalized data trace.
    """
    return (data-lower_ref)/(upper_ref-lower_ref)

def Wcp(f,alpha,n,tau):
	"""Filter function of the ideal CPMG sequence.
	
	:param f: frequency of signal
	:param alpha: Phase of signal
	:param n: number of pulses
	:param tau : inter-pulse spacing
	
	:returns: The filter function.

	See https://doi.org/10.1103/RevModPhys.89.035002 page 18
	"""
	if f == 0:
		return 0
	tmp = _np.sin(_np.pi*f*n*tau)/(_np.pi*f*n*tau)
	tmp = tmp*(1-1/_np.cos(_np.pi*f*tau))
	tmp = tmp*_np.cos(alpha+_np.pi*f*n*tau)
	return tmp

# def plot_NV_trace(*args,**kwargs):
#     import matplotlib.pyplot as plt
#     plt.axhline(1,color='k')
#     plt.axhline(0,color='k')
#     plt.axhline(0.5,color='k',alpha=0.5,linestyle='--')
#     plt.plot(*args,**kwargs)