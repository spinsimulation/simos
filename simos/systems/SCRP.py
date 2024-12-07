import numpy as _np
import scipy as _sc
from ..constants import ye
from ..qmatrixmethods import *
from ..core import subsystem
from ..core import reverse_subsystem
from ..utils import lor
from ..utils import dlor
from .. import backends 
from ..propagation import rot, evol, rotate_operator
import copy

######################################################################
# Please note  - this is an SCRP specific submodule. 
# Functions are not necessarily applicable to other spin/level systems.
######################################################################

############################################################
# Initial States of SCRPs
############################################################

def initial_state(system, e1, e2, alpha, beta, delta, theta = 0, phi = 0):
    """ Initial density matrix of a spin-correlated radical pair using the parameterization introduced by Völker et al. in 
    J. Chem. Phys. 158, 161103 (2023).
    
    :param System system: An instance of the System-class.
    :param str e1: Name of the first electron spin of the SCRP in the system.
    :param str e2: Name of the second electron spin of the SCRP in the system. 
    :param alpha: Alpha parameter.
    :param beta: Beta parameter.
    :param delta: Delta parameter.
    :param theta: Angle between z-axis of laboratory frame and molecular axis.
    :param phi: Angle between x-axis of laboratory frame and molecular axis projection onto laboratory frame xy plane. 

    :returns: The density matrix of the initial state. 
    """
    # STEP 1  generate rho0 in the molecular frame of reference.
    # Symbolic required?
    if all(backends.get_backend(a) != 'sympy' for a in [alpha, beta, delta, theta, phi]):
        mode = "numeric"
    else:
        mode = "symbolic"  
    sin_fun = backends.get_calcmethod("sin", mode)
    cos_fun = backends.get_calcmethod("cos", mode)
    sqrt_fun = backends.get_calcmethod("sqrt", mode)
    exp_fun = backends.get_calcmethod("exp", mode)
    I = backends.get_calcmethod("I", mode)
    # Get the required operators (in the subspace of the radical pair).
    p1  = getattr(system, e1+"p")
    p2 =  getattr(system, e2 +"p")
    ab, rest, reverse =  subsystem(system, (p1[0.5]*p2[-0.5]).unit(),  [e1, e2], keep= True)
    ba, rest , reverse = subsystem(system, (p1[-0.5]*p2[0.5]).unit(),  [e1, e2], keep= True)
    # Construct the state (rho0) in the molecular frame of reference. 
    ket0 = 1/sqrt_fun(2)*(cos_fun(alpha) + exp_fun(I*beta)*sin_fun(alpha))* dm2ket(ab.unit())  - 1/sqrt_fun(2)*(cos_fun(alpha) - exp_fun(I*beta)*sin_fun(alpha)) * dm2ket(ba.unit())
    rho0 = ket2dm(ket0)
    scaler = delta * _np.ones(rho0.shape)
    _np.fill_diagonal(scaler, _np.ones(rho0.shape[0]))
    rho0 = tidyup(data(rho0)*data(scaler), method = system.method, dims = rho0.dims)
    rho0full = reverse_subsystem(rho0, rest, reverse)
    # STEP 2 rotate rho0 to the laboratory frame of reference.
    rho0full = rotate_operator(system, rho0full, -1*phi, 0, -1*theta)
    return rho0full


################################################
# Helper Functions for EPR Experiments on SCRPs
###############################################

def cwEPR_fieldsweep(system, selection, initial_state, Hfield, Hrest, omega, broadening = 0, axis = None, weights = None):
    """ Calculates the fieldswept EPR spectrum of a radical pair. Eigenfields are obtained 
    by diagonalizing the Hamiltonian in a higher dimensional Liouville-type space, the amplitudes of the individual 
    resonance fields are calculated as products of transition polarizations and transition probabilities. 

    :param system: An instance of the System class.
    :param e1: Name of the first electron spin.
    :param e2: Name of the second electron spin. 
    :param initial_state: Density matrix.
    :param Hamiltonian: System Hamiltonian.
    :param broadening: Line broadening (1/s).
    
    :returns: Axis and spectrum.
    """
    
    # This method does not work for sympy.
    if system.method == "sympy":
        raise ValueError("Sorry, this method is not compatible with the sympy backend.")
    # Define operator for EPR excitation.
    # Use initial state polarisation to weight the excitation operator.
    ex = 0*system.id
    for member in selection:
        ex = ex +  getattr(system, member +"x") 
    excite , _ , _ = subsystem(system, ex.unit() , selection, keep= True)
    rho0   , _ , _ = subsystem(system, initial_state.unit(),  selection, keep= True)
    subdim = rho0.shape[0]
    fac = subdim/system.dim
    idsub = identity(subdim, method = system.method)
    # Loop over ensemble components.
    fields = _np.zeros(len(Hfield)*system.dim**2, dtype = complex)
    intensities = _np.zeros(len(Hfield)*system.dim**2, dtype = complex)
    for i in range(len(Hfield)):
        F , _, _  = subsystem(system, Hrest[i],  selection, keep= True)
        G, _ , _ = subsystem(system, Hfield[i],  selection, keep= True)
        F = fac*F
        G = fac*G
        # Calculate the Eigenfields (and their eigenvectors).
        B = spre(G) - spost(G)
        A = spre(omega*idsub) - (spre(F) - spost(F))
        vals, vecs = _sc.linalg.eig(data(A), data(B)) 
        # Calculate the Intensities.
        intensitiestmp = _np.zeros(len(vals), dtype = complex )
        for int_ind, val in enumerate(vals):
            # Only take into account fields that are > 0 und < infty. 
            if _np.abs(val) < _np.infty and  val > 0 :
                # Get Hamiltonian at this field and eigenvalues. 
                Htot = F + G*val
                truevals, truevecs =  Htot.eigenstates()
                # Set up the transformation matrix.
                T_toeigen = _np.zeros(Htot.shape, dtype = complex) 
                for ind_vec, vec in enumerate(truevecs):
                    T_toeigen[ind_vec,:] = data(vec)[:,0]
                # Sum up all transition probabilities 
                ex_red_eigen  = excite.transform(T_toeigen)
                rho0_red_eigen  = rho0.transform(T_toeigen)
                amps = []
                for ind_val in range(len(truevals)):
                    for ind_val2 in range(ind_val+1, len(truevals)):
                        popdiff = (rho0_red_eigen[ind_val2, ind_val2]) - (rho0_red_eigen[ind_val, ind_val])
                        #popdiff = (rho0_red_eigen[ind_val2, ind_val2])**2 - (rho0_red_eigen[ind_val, ind_val])**2                       
                        erlaubt = (ex_red_eigen[ind_val, ind_val2])**2
                        if (_np.abs((truevals[ind_val2]-truevals[ind_val])- omega)) < 0.00001*omega:
                            amps.append(-erlaubt*popdiff)
                intensitiestmp[int_ind] = _np.sum(amps)
            else:
                intensitiestmp[int_ind] = 0
        fields[i*subdim**2:i*subdim**2+subdim**2] = vals
        intensities[i*subdim**2:i*subdim**2+subdim**2] = intensitiestmp 
    # Calculate the spectrum.
    if axis is None:
        minval = _np.min([_np.abs(field) for field in fields  if _np.abs(field) < _np.infty])
        maxval = _np.max([_np.abs(field) for field in fields  if _np.abs(field) < _np.infty])        
        axis = _np.linspace(0.5*minval, 2*maxval, 1000)
    spectrum = _np.zeros(_np.shape(axis), dtype = complex)
    for ind_B, B in enumerate(fields):
        #if _np.abs(intensities[ind_B]) >  1e-16 and _np.abs(B) < _np.infty: 
        if _np.abs(B) < _np.infty: 
            spectrum = spectrum + lor(axis, intensities[ind_B], B, broadening)
    #if _np.sum(_np.abs(spectrum)) > 1e-10:
    #    spectrum = spectrum/_np.sum(_np.abs(spectrum))
    return axis, spectrum, fields, intensities, rho0 

def cwEPR_frequencysweep(system, e1, e2, initial_state, Hamiltonian, broadening = 0, axis = None, weights = None):
    """ Calculates the frequency swept EPR spectrum of a radical pair. Frequencies are obtained 
    by diagonalizing the Hamiltonian, the amplitudes of the individual transitions are calculated as products 
    of transition polarizations and transition probabilities. 

    :param system: An instance of the System class.
    :param e1: Name of the first electron spin.
    :param e2: Name of the second electron spin. 
    :param initial_state: Density matrix.
    :param Hamiltonian: System Hamiltonian.
    :param broadening: Line broadening (1/s).
    
    :returns: Axis and spectrum.
    """
    # This method does not work for sympy.
    if system.method == "sympy":
        raise ValueError("Sorry, this method is not compatible with the sympy backend.")
    # Define operator for EPR excitation.
    ex =  getattr(system, e1+"x") + getattr(system, e2 +"x")
    # Reduce our system to the radical pair
    rho0_red, _ , reverse = subsystem(system, initial_state,  [e1, e2], keep= True)
    ex_red, _ , _ = subsystem(system, ex , [e1, e2], keep= True)
    # Prepare data storage.
    amps = []
    ws = []
    # Loop over list of Hamiltonians and process. 
    for ind_H, H in enumerate(Hamiltonian):
        H_red, _ , _ =  subsystem(system,  H,  [e1, e2], keep= True) 
        # Calculate Eigenvalues and Eigenvectors of the system Hamiltonian. 
        vals, vecs = H_red.eigenstates()
        # Set up the transformation matrix.
        T_toeigen = _np.zeros(H_red.shape, dtype = complex) # set up the transformation matrix using the eigenvectors 
        for ind_s, s in enumerate(vecs):
            T_toeigen[ind_s,:] = data(s)[:,0]
        # Multiply the frequencies with the transition probabilities.
        ex_red_eigen  = ex_red.transform(T_toeigen)
        rho0_red_eigen  = rho0_red.transform(T_toeigen)
        for ind_val in range(len(vals)):
            val = vals[ind_val]
            for ind_val2 in range(ind_val+1, len(vals)):
                val2 = vals[ind_val2]
                popdiff = rho0_red_eigen[ind_val2, ind_val2] - rho0_red_eigen[ind_val, ind_val]
                erlaubt =(ex_red_eigen[ind_val, ind_val2])**2
                if weights is None:
                    amps.append(popdiff*erlaubt)
                else:
                    amps.append(popdiff*erlaubt*weights[ind_H])
                ws.append(val2-val)
    # Set up frequency axis and calculate the spectrum.
    if axis is None:
        axis = _np.linspace(_np.min(ws)/100, _np.max(ws)+2*broadening, 1000)
    broadening = _np.max([_np.abs(broadening), axis[1]-axis[0]])
    spectrum = _np.zeros(_np.shape(axis), dtype = complex)
    for ind_w, w in enumerate(ws):
        if _np.abs(amps[ind_w]) >  1e-16: 
            spectrum = spectrum + lor(axis, amps[ind_w], w, broadening)
    if _np.sum(_np.abs(spectrum)) > 1e-10:
        spectrum = spectrum/_np.sum(_np.abs(spectrum))
    return axis, spectrum, ws, amps

def cwEPR_transient(system, e1, e2, initial_state, fields, Hfield, Hrest, omega, tax, Bax, broadening = 0 ):
    raise NotImplementedError

def ESEEM(system, e1, e2, rho0, H, taus, flip1 = _np.deg2rad(45), flip2 = _np.deg2rad(180), phase1 = "x", phase2 = "x", detector = "y", c_ops = []):
    results = _np.array(len(taus), dtype = complex)
    pul1op = getattr(system, e1+phase1)+getattr(system, e2+phase1)
    pul2op = getattr(system, e1+phase2)+getattr(system, e2+phase2)
    measop = getattr(system, e1+detector) + getattr(system, e2+detector)
    # Excite
    rho0 = rho0.copy()
    rho0 = rot(pul1op, flip1, rho0)
    # Loop over pulse delays.
    for ind_t, t in enumerate(taus):
        # First free evolution period.
        rho = rho0.copy()
        rho = evol(H, t, rho, c_ops = c_ops)
        # Pi pulses.
        rho =  rot(pul2op, flip2, rho)
        # Second free evolution period
        rho =  rot((getattr(system, e1+phase2)+getattr(system, e2+phase2)), flip2, rho)
        # Measure
        results[ind_t] = expect(measop, rho)
    return results


#################################
# Semiclassical Helpers 
#################################

def schultenwoynes_randomfields(As, Is, N):
    """ Returns a list of field vectors for a series of nuclear spins with given hyperfine coupling
    tensors. 
    
    :params As: A list of hyperfine coupling tensors.
    :params Is: A list of spin quantum numbers.
    :params N: Number of random fields.
    
    :returns: N field vectors to be used with Schulten Woynes semiclassical method. """
    # Length of nuclear spin vectors.
    lIs = _np.array([_np.sqrt(i*(i+1)) for i in Is])
    hfivec = _np.empty((N,3))
    hfinorm = _np.empty(N)
    for num_config in range(N):
        tot_matrix = _np.empty((3))
        # Loop over all spins. For each spin, get
        # a random orientation and calculate the resulting hyperfine field
        # at the nuclear site.  
        for spin in range(len(Is)): 
            pos = _sc.stats.uniform_direction.rvs(3)*lIs[spin]
            tot_matrix += _np.matmul(_np.transpose(pos), As[spin])
        hfivec[num_config, :] = tot_matrix
        hfinorm[num_config] = _np.linalg.norm(tot_matrix)
    return hfivec, hfinorm


#################################
# Singlet and Triplet Yield
#################################


def ST_yields(system, e1, e2, rhos, dt, **kwargs):

    # Prepare Singlet and Triplet States
    system = copy.copy(system) 
    system.add_ghostspin("C", [e1, e2])
    S =  system.C_1p[0].unit()
    states = [system.C_1p[0].unit(), system.C_3p[0].unit(), system.C_3p[1].unit(), system.C_3p[-1].unit()]
    names = ["S", "T0", "T1", "Tm1"]
    PST = {}
    # Calculate and store yields. 
    if "k" in kwargs:
        decay = True
    else:
        decay = False
    for ind_s, state in enumerate(states):
        PST[names[ind_s]] = _yield(state, rhos, dt, kwargs.get("k"), decay_accounting = decay)
    return PST

def _yield(state, rhos, dt, k, decay_accounting = False):
    srho = _np.shape(rhos)
    # Define a time axis.
    tax = _np.arange(0, srho[-1])*dt
    # Calculate  population in all density matrices.
    P = _np.empty(srho, dtype = complex)
    phi = _np.empty(srho[0], dtype = complex)
    if len(srho) == 2:
        for i in range(srho[0]):
            for j in range(srho[1]):
                P[i,j] = expect(state, rhos[i,j])
            if decay_accounting:
                P[i, :] =  _np.multiply(P[i,:], _np.exp(-k*tax)) 
            phi[i] = k*_np.sum(P[i,:]*dt, axis = 0)
    elif len(srho) == 1:
        for j in range(srho[1]):
            P[j] = expect(state, rhos[j])
        if decay_accounting:
            P =  _np.multiply(P, _np.exp(-k*tax)) 
        phi = k*_np.sum(P*dt, axis = 0)
    else:
        raise ValueError("Invalid shape for rho input.")  
    # Return time evolution and total yield. 
    return P , phi


