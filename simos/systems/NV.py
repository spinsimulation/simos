import numpy as _np
from ..constants import ye, yC13, eps_z, eps_xy, eps_0, hbar, yN15, yH1, yF19
from ..trivial import spher2cart, cart2spher
from ..propagation import evol, rot
from ..qmatrixmethods import applySuperoperator, ptrace, expect, tensor
from .. import backends
from ..core import subsystem, reverse_subsystem
from ..states import pol_spin
from itertools import combinations
from ..coherent import dipolar_coupling, dipolar_spatial
from ..trivial import f2w, w2f
###########################################################
# INTERACTIONS
###########################################################

def efield(*args, eps_rel = 1, mode='cart'):
    """ Calculates the electric field of individual charges as well as total electric field in [V/m]
    Vectorized
    """
    # Argument structure handling
    if len(args) == 4:
        coords = _np.array([args[0], args[1], args[2]]).transpose()
        q = args[3]
    elif len(args) == 2:
        coords = _np.array(args[0])
        q = args[1]
    else:
        raise ValueError("Wrong number of coordinates provided")
    if not isinstance(q, _np.ndarray):
         q = _np.array(q)
    # Spher2cart if necessary
    if mode =='spher':
        coords = spher2cart(coords)
    if len(coords.shape) == 1: # pack if single vector provided 
        coords = _np.array([coords])
    # Calculate Efield
    normalizer = _np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)
    pre = q/(4*_np.pi*eps_0*eps_rel)
    E =  _np.einsum('a, ab -> ab', pre/normalizer**3, coords)
    Etot = _np.sum(E, axis = 0)
    return E, Etot

def stark_interaction(spinsystem, Name, efield):
    """ Calculates stark Hamiltonian in angular frequencies
    """
    # Get Spin Operators
    Sz = getattr(spinsystem, Name + 'z')
    Sx = getattr(spinsystem, Name + 'x')
    Sy = getattr(spinsystem, Name+ 'y') 
    # Assemble stark Hamiltonian
    Ex = efield[0]
    Ey = efield[1]
    Ez = efield[2]
    Hstarkaxial = 2*_np.pi*eps_z*Ez*Sz**2 
    Hstarknonaxial = - 2*_np.pi*eps_xy*Ex*(Sx**2-Sy**2) +  2*_np.pi*eps_xy*Ey*(Sx*Sy + Sy*Sx)
    Hstark = Hstarkaxial + Hstarknonaxial

    return Hstark

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
    """Generates the initial state of the NV-Nuclear spin system
    For the NV center, the initial state is the |0> state of the NV center.
    The nuclear spins are initialized to polarized states specified by 'pol' in the spinsystem.
    """
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
                rho_arr.append(pol_spin(spinsystem.system[i]['pol']))
            else:
                 rho_arr.append(pol_spin(1).unit())
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
# Helper functions
###########################################################

def meas_NV(rho, spinsystem, NV_name='S'):
    """Perform measurement operation on NV
    
    rho        : Density matrix of the state to measure
    spinsystem : Spinsystem class object containing all operators
    NV_name    : NV center spin name. If None, search for all NVs

    Returns expectation value, rho after measurement
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

def plot_NV_trace(*args,**kwargs):
    import matplotlib.pyplot as plt
    plt.axhline(1,color='k')
    plt.axhline(0,color='k')
    plt.axhline(0.5,color='k',alpha=0.5,linestyle='--')
    plt.plot(*args,**kwargs)

def exp2cts(x,contrast,ref_cts,ref=0.5):
    m = contrast*ref_cts/(1-contrast+contrast*ref)
    c = -ref_cts*(-1 + contrast)/(1-contrast+contrast*ref)
    return m*x+c

def _noise_gen(xi):
    return _np.random.poisson(xi) #+ xi
add_shot_noise = _np.vectorize(_noise_gen)

def normalize_data(data,upper_ref,lower_ref):
    """Normalizes the data such that the upper reference correponds to 1, the lower reference to 0
    
    Parameters:
    data      : Data trace to normalize
    upper_ref : mean value of the upper reference
    lower_ref : mean value of the lower reference
    """
    return (data-lower_ref)/(upper_ref-lower_ref)

def Wcp(f,alpha,n,tau):
	"""Filter function of the ideal CPMG sequence
	
	Parameter:
	f     : frequency of signal
	alpha : Phase of signal
	n     : number of pulses
	tau   : inter-pulse spacing
	
	
	See https://doi.org/10.1103/RevModPhys.89.035002 page 18
	"""
	if f == 0:
		return 0
	tmp = _np.sin(_np.pi*f*n*tau)/(_np.pi*f*n*tau)
	tmp = tmp*(1-1/_np.cos(_np.pi*f*tau))
	tmp = tmp*_np.cos(alpha+_np.pi*f*n*tau)
	return tmp