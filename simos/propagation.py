import qutip as _qu
import scipy as _sc
import numpy as _np
import scipy.linalg as la
import types
from .constants import *
from .qmatrixmethods import expm
from .core import *
from . import backends
import types
import sys
from typing import Union
from .core import globalclock
import warnings


###########################################################
# TIME PROPAGATION OPERATIONS
###########################################################

def rot(op, phi,*rho):
    """ Rotation of a spin on the Bloch sphere. Returns rho' if rho is given. 

    This operation corresponds to an infinitely fast pulse of the spin operator along the axis defined by the operator op by an angle phi.

    :param op: Operator for the rotation (use proper normalization and appropriate projection)
    :param phi: angle to rotate in radians
    :param *rho: Input state (optional)
    :return: Rotated state. If rho is not specified, then the rotation operator is returned.
    """
    # rotation oerator is (-1j*phi*op).expm()
    # We use the evol implementation
    return evol(op,phi,*rho)

def evol(H, t, *rho, c_ops=[], proj=[],wallclock='global'):
    """  Evolve state rho under static conditions for time t and return rho' if rho is given.

    :param H: Hamiltonian in angular frequencies for the evolution. Can be a quantum object or None if c_ops are given.
    :param t: Evolution time.
    :param *rho: Input state (optional)
    :return: Evolved state. If rho is not specified, the evolution operator is returned. In cases where c_ops are given, the vector superoperator is returned.

    :Keyword Arguments:
        * *c_ops* (``list``) --
            List of collapse operators. If c_ops is given, the evolution is following the Lindblad master equation.
        * *wallclock* --
            If 'global', the global wallclock is used to track the time. If None, no time tracking is performed. If a wallclock object is given, it is used to track the time.
        * *proj* (``list``) --
            List of projectors for space reduction. Currently not implemented, reserved for future use.
    """
    # Fälle zum Unterscheiden:
    # H & t -> Unitäre Evolution
    # H & t & cops -> Lindblad Evolution
    if proj != []:
        raise NotImplementedError('Projection is not yet implemented.')
    if H is not None:
        method = backends.get_backend(H)
    else:
        method = backends.get_backend(c_ops[0])
    
    if method == "sympy":
        I  = backends.get_calcmethod('I', 'symbolic')
    else:
        I  = backends.get_calcmethod('I', 'numeric')

    if wallclock == 'global':
        clock = globalclock
    else:
        clock = wallclock
    if clock is not None:
        clock.inc(t)
    
    if len(c_ops) == 0 and H is not None:
        # Hilbert space evolution
        #print(I,t,H)
        Ut = (-I*t*H)
        if hasattr(H,'expm'):
            #print("Attribute called")
            Ut = Ut.expm()
        else:
            #print("Function called")
            Ut = expm(Ut)
        if len(rho) == 0:
            # Propagator return
            return Ut
        elif len(rho) == 1:
            rho = rho[0]
            # check if rho is a ket or a density matrix
            isket = (rho.shape[1] == 1)
            if isket:
                return Ut*rho
            else:
                return Ut*rho*Ut.dag()
        else:
            raise ValueError('More than one rho was given. This is not supported.')
    elif len(c_ops) > 0:
        liouvillian = getattr(getattr(backends,method), 'liouvillian')
        L = liouvillian(H,c_ops)
        
        UL = (L*t)
        if hasattr(UL,'expm'):
            UL = UL.expm()
        else:
            #expm_fun = getattr(getattr(backends,method), 'expm')
            #UL = expm_fun(UL)
            UL = expm(UL)
        
        if len(rho) == 0:
            # Propagator return
            return UL
        elif len(rho) == 1:
            rho = rho[0]
            # check if rho is a ket or a density matrix
            isket = (rho.shape[1] == 1)
            applySuperoperator = getattr(getattr(backends,method), 'applySuperoperator')
            if isket:
                ket2dm = getattr(getattr(backends,method), 'ket2dm')
                rho = ket2dm(rho)
            return applySuperoperator(UL,rho)



def prop(H0,dt,*rho,c_ops=[],H1=None,carr1=None,c_ops2=[],carr2=None,proj=[],engine='cpu',wallclock = 'global',**kwargs):
    """  Evolve state for time-dependent conditions for given parametrization and return rho' if rho is given.

    :param H0: Static background Hamiltonian in angular frequencies for the evolution. Can be a quantum object or None. Note, at least one of H0, H1, c_ops or c_ops2 must be given.
    :param dt: Time step for the modulation array(s) for the evolution.
    :param *rho: Input state (optional)
    :return: Evolved state. If rho is not specified, the evolution operator is returned. In cases where c_ops are given, the vector superoperator is returned.

    :Keyword Arguments:
        * *c_ops* (``list``) --
            List of collapse operators. If c_ops is given, the evolution is following the Lindblad master equation.
        * *H1* (``list``) --
            List of time-dependent Hamiltonians. If H1 is given, carr1 must be given as well. The number of Hamiltonians in H1 and the number of modulation arrays in carr1 must be equal.
        * *carr1* (``list``) --
            List of coefficients for the time-dependent Hamiltonians. If carr1 is given, H1 must be given as well. The number of arrays in carr1 and the number of Hamiltonians in H1 must be equal. Discretization of the time-dependent Hamiltonians is done by the time step dt.
        * *c_ops2* (``list``) --
            List of time-dependent collapse operators. If c_ops2 is given, carr2 must be given as well. The number of collapse operators in c_ops2 and the number of modulation arrays in carr2 must be equal.
        * *carr2* (``list``) --
            List of coefficients for the time-dependent collapse operators. If carr2 is given, c_ops2 must be given as well. The number of arrays in carr2 and the number of collapse operators in c_ops2 must be equal. Discretization of the time-dependent collapse operators is done by the time step dt.
        * *engine* (``str``) --
            Engine to be used for the evolution. Can be 'cpu', 'qutip' or 'parament'. If 'qutip' is chosen, mesolve from qutip is used for the evolution. This only works if the qutip backend is used. If 'parament' is chosen, the parament library is used for the evolution. This only works if the parament backend is installed and a GPU is available to use.
        * *wallclock* --
            If 'global', the global wallclock is used to track the time. If None, no time tracking is performed. If a wallclock object is given, it is used to track the time.
        * *proj* (``list``) --
            List of projectors for space reduction. Currently not implemented, reserved for future use.
        * *kwargs* --
            Additional keyword arguments for the evolution. Will be passed to the used engine. For the qutip engine, e.g., this can be the options for the mesolve function.

    This operation corresponds to an infinitely fast pulse of the spin operator along the axis defined by the operator op by an angle phi.    
    """

    # Time-dependent Hamiltonian
    # Methods are: direct, mesolve (only for qutip), parament, 
    if proj != []:
        raise NotImplementedError('Projection is not yet implemented.')
    if H1 is None:
        H1 = []
    if carr1 is None:
        carr1 = []
    if c_ops2 is None:
        c_ops2 = []
    if carr2 is None:
        carr2 = []

    # method and backend compatibility
    backend = backends.get_backend(H0)
    dims = None
    if hasattr(H0,'dims'):
        dims = H0.dims
    if engine == 'qutip' and backend != 'qutip':
        raise ValueError('The qutip engine is only available for qutip backends.')
    if backend == 'sympy':
        raise NotImplementedError('The sympy backend is (and will) not supported by the prop function.')
    if backend == 'sparse':
        raise NotImplementedError('The sparse backend is currently not implemented by the prop function.')
    tidyup_fun = getattr(getattr(backends, backend), 'tidyup')
    data_fun =  getattr(getattr(backends, backend), 'data')
    isket = getattr(getattr(backends, backend), 'isket')
    # Inconsistent input
    # Raise error if H1 is given but carr1 is not and vice versa
    # Raise error if c_ops2 is given but carr2 is not and vice versa
    # Not given means []
    if (H1 is not []) and (carr1 is []):
        raise ValueError('If H1 is given, carr1 must be given as well.')
    if (H1 is []) and (carr1 is not []):
        raise ValueError('If carr1 is given, H1 must be given as well.')
    if (c_ops2 is not []) and (carr2 is []):
        raise ValueError('If c_ops2 is given, carr2 must be given as well.')
    if (c_ops2 is []) and (carr2 is not []):
        raise ValueError('If carr2 is given, c_ops2 must be given as well.')
    
    # Consistency in input: H1 and c_ops2 must be lists
    if not isinstance(H1,(list,tuple)):
        H1 = [H1]
    if not isinstance(c_ops2,(list,tuple)):
        print("c_ops2",c_ops2,type(c_ops2))
        c_ops2 = [c_ops2]

    if not isinstance(carr1,(list,tuple)):
        carr1 = [carr1]
    if not isinstance(carr2,(list,tuple)):
        carr2 = [carr2]

    # Asset that len(H1) == len(carr1) and len(c_ops2) == len(carr2)
    if len(H1) != len(carr1):
        print(H1,_np.shape(H1),type(H1))
        print(carr1,_np.shape(carr1),type(carr1))
        raise ValueError('The number of Hamiltonians in H1 and coefficients in carr1 must be equal.',len(H1),"!=",len(carr1))
    if len(c_ops2) != len(carr2):
        print(c_ops2)
        print(carr2)
        raise ValueError('The number of collapse operators in c_ops2 and coefficients in carr2 must be equal.',len(c_ops2),"!=",len(carr2))
    
    
    # Static case
    if H1 == [] and c_ops2 == []:
        return evol(H0,dt,*rho,c_ops=c_ops,proj=proj,wallclock=wallclock)
    
    N = max([len(H1),len(carr1),len(c_ops2),len(carr2)])
    
    if wallclock == 'global':
        clock = globalclock
    else:
        clock = wallclock
    if clock is not None:
        clock.inc(dt*N)    
    
    # Only time varying Hamiltonian and qutip
    if H1 != [] and c_ops2 == []:
        if engine == 'qutip':
            if len(rho) == 0:
                raise ValueError('For the qutip engine, a state vector or density matrix must be provided.')
            else:
                qu = getattr(getattr(backends,'qutip'),'_qu')
                PTS = _np.shape(carr1)[1]
                times = _np.arange(PTS)*dt
                # Make out of the lists H1 and carr1 a list like [[H1[0],carr1[0]],[H1[1],carr1[1]],...]
                coeff = [[H1[i],carr1[i]] for i in range(len(H1))]
                
                #print("coeff",_np.shape(_np.array(coeff)))
                #print("times",_np.shape(_np.array(times)))
                qevo = qu.QobjEvo([H0] + coeff, tlist=times, order=1)
                result =  qu.mesolve(qevo,rho[0],times,c_ops=c_ops,**kwargs)
                return result.states[-1]
            
    if engine == 'RK45':
        # ensure that rhoin is a density matrix
        if len(rho) == 0:
            raise ValueError('For the RK45 engine, a state vector or density matrix must be provided.')
        if isket(rho[0]):
            ket2dm = getattr(getattr(backends,backend), 'ket2dm')
            rhoin = ket2dm(rhoin)
        else:
            rhoin = rho[0]
        rhoout = solve_lindblad(data_fun(rhoin),dt,H0=data_fun(H0),H1=[data_fun(H1[i]) for i in range(len(H1))],c_ops=[data_fun(c_ops[i]) for i in range(len(c_ops))],carr1=carr1,c_ops2=[data_fun(c_ops2[i]) for i in range(len(c_ops2))],carr2=carr2)
        return tidyup_fun(rhoout,dims=dims)

    def build_timedependent_liouvillian(H0,c_ops,carr1,c_ops2,carr2,backend):
        liouvillian = getattr(getattr(backends,backend), 'liouvillian')
        if c_ops2 == None:
            c_ops2 = []
            carr2 = []
        if len(c_ops2) == 0:
            return liouvillian(H0,c_ops),[liouvillian(H1[i],[]) for i in range(len(H1))],carr1
        else:
            lindbladian = getattr(getattr(backends,backend), 'lindbladian')
            L = liouvillian(H0,c_ops)
            L1 = [liouvillian(H1[i],[]) for i in range(len(H1))]
            L2 = [lindbladian(c_ops2[i]) for i in range(len(c_ops2))]
            return L,L1+L2,_np.concatenate((carr1,carr2))
        
    # Check if we have to go to Liouvillian
    if (len(c_ops) > 0) or (len(c_ops2) > 0):
        Hm0,Hm1,cmarr = build_timedependent_liouvillian(H0,c_ops,carr1,c_ops2,carr2,backend)
    else:
        Hm0, Hm1, cmarr = H0, H1, carr1

    if engine == 'cpu':
        Ut = propagate_cpu(Hm0,Hm1,cmarr,dt,**kwargs)
        Ut = tidyup_fun(Ut,dims=dims)
    elif engine == 'parament':
        import parament
        handle = parament.Parament()
        H0 = data_fun(H0)
        H1 = [data_fun(H1[i]) for i in range(len(H1))]
        handle.set_hamiltonian(H0,*H1)
        Ut = handle.equiprop(dt,carr1)
        Ut = tidyup_fun(Ut,dims=dims)
        handle.destroy()
    
    if len(rho) == 0:
        # Propagator return
        return Ut
    elif len(rho) == 1:
        rho = rho[0]
        # check if it is a Liouvillian
        if (len(c_ops) > 0) or (len(c_ops2) > 0):
            applySuperoperator = getattr(getattr(backends,backend), 'applySuperoperator')
            if isket(rho):
                ket2dm = getattr(getattr(backends,backend), 'ket2dm')
                rho = ket2dm(rho)
            return applySuperoperator(Ut,rho)
        #Check if it is a ket
        isket = (rho.shape[1] == 1)
        if isket:
            return Ut*rho
        else:
            return Ut*rho*Ut.dag()


  
def generate_simpson(cin,dt=None,magnus=True,second_order=True):
    """Generate the Simpson coefficients for the time-dependent Hamiltonian. Compute also modulation arrays for the Magnus expansion.

    :param cin: Coefficients for the time-dependent Hamiltonian.
    :param dt: Time step for the modulation array(s) for the evolution. If dt is None, the time step is set to 1.
    :param magnus: If True, the modulation arrays for the Magnus expansion commutators are computed.
    :param second_order: If True, the modulation arrays for the second order Magnus expansion commutators are computed. These are commutators between modulations array functions itself. For many modulation terms, this can be computationally expensive.    
    """
    # eq 14 parament paper
    c = _np.atleast_2d(cin)
    amps = c.shape[0]
    n = c.shape[1]
    amps_new = amps
    if magnus:
        amps_new = amps*2
    else:
        second_order = False
    if second_order:
        amps_new = amps*2 + amps*(amps-1)//2
        print(amps_new)
    len_new = (n - 3) // 2 + 1
    out = _np.zeros((amps_new,len_new),dtype=_np.complex128)
    idx =  _np.arange(len_new)
    for j in range(amps):
        out[j,:] = (c[j,2 * idx]  + 4*c[j,2 * idx +1] + c[j,2*idx+2]) / 6
      
    if magnus:
        for j in range(amps):
            out[j+amps,:] = 1j*(c[j,idx+2] - c[j,idx]) / 6 * dt
    
    if second_order:
        l = 0
        for j in range(amps):
            for k in range(j):
                out[2*amps+l] = 1j*(c[k,idx]*c[j,idx+2]  - c[k,idx+2]*c[j,idx])*dt/6
                l+=1

    return out

def generate_magnus_commutators(H0,H1,magnus=True,second_order=True):
    """Generate the commutators for the Magnus expansion.

    :param H0: Static background Hamiltonian in angular frequencies for the evolution.
    :param H1: List of time-dependent Hamiltonians.
    :param magnus: If True, the Magnus expansion is used.
    :param second_order: If True, the second order Magnus expansion is used (commutators between the modulation Hamiltonians). Can be computationally expensive for many modulation terms.
    """
    H1 = _np.atleast_3d(H1)
    amps = H1.shape[0]
    newn = amps
    if magnus:
        newn = 2*amps
    else:
        second_order = False
    if second_order:
        newn = amps*2 + amps*(amps-1)//2
    H1out = _np.zeros((newn,H1.shape[1],H1.shape[2]),dtype=_np.complex128)
    H1out[:amps,:,:] = H1[:,:,:]
    if magnus:
        H1out[amps:2*amps,:,:] = (H0 @ H1) - (H1 @ H0)
    
    if second_order:
        l = 0
        for j in range(amps):
            for k in range(j):
                H1out[2*amps+l,:,:] = H1[k,:,:] @ H1[j,:,:]  - H1[j,:,:] @ H1[k,:,:]
                l+=1
    return H1out

def propagate_cpu(H0,H1,carr,dt,maxsize=1,**kwargs):
    """Propagate the state vector under the time-dependent Hamiltonian using the CPU.

    :param H0: Static background Hamiltonian in angular frequencies for the evolution. Can be a quantum object or None.
    :param H1: List of time-dependent Hamiltonians.
    :param carr: List of coefficients for the time-dependent Hamiltonians.
    :param dt: Time step for the modulation array(s) for the evolution.
    :param maxsize: Maximum size of the memory in GB. If the memory exceeds this size, the batch size is reduced.
    :param **kwargs: Additional keyword arguments. The keyword 'magnus' can be used to set the Magnus expansion. If 'magnus' is True, the Magnus expansion is used. If 'second_order' is True, the second order Magnus expansion is used. If 'second_order' is False, the first order Magnus expansion is used. If 'magnus' is False, the Magnus expansion is not used.
    """
    backend = backends.get_backend(H0)
    tidyup_fun = getattr(getattr(backends, backend), 'tidyup')
    data_fun =  getattr(getattr(backends, backend), 'data')
    dims = None
    if hasattr(H0,'dims'):
        dims = H0.dims
    
    H0 = data_fun(H0)
    H1 = [data_fun(H1[i]) for i in range(len(H1))]

    H1 = _np.atleast_3d(H1)
    carr = _np.atleast_2d(carr)
    mem = maxsize*1024*1024*1024
    perMat = _np.prod(H0.shape)*128/8
    batchmax = int(mem // perMat)
    carr = _np.flip(carr,axis=-1)

    # Check if magnus is set
    if 'magnus' in kwargs:
        magnus = kwargs['magnus']
    else:
        magnus = False
    if 'second_order' in kwargs:
        second_order = kwargs['second_order']
    else:
        second_order = True
    
    if magnus:
        H1 = generate_magnus_commutators(H0,H1,magnus,second_order)
        carr = generate_simpson(carr,dt,magnus,second_order)
        dt = dt*2

    U = _np.eye(H0.shape[0],dtype=_np.complex128)
    pts = len(carr)
    for i in range(0,pts,batchmax):
        currlen = int(_np.clip(batchmax,0,pts-i))
        Hc = _np.tensordot(carr[i:i+currlen],H1,axes=(0,0))
        Hc = Hc + H0  
        #print(Hc,Hc.shape,type(Hc))
        Uc = la.expm(-1j*dt*Hc)
        
        #Uc = _np.flip(Uc,axis=0)
        U = functools.reduce(_np.dot,Uc) @ U
    tidyup_fun(U,dims=dims)
    return U

def euler_rodrigues(op,tol=1e-10):
    """Return (alpha, rotation_axis) of a given propagator

    This routine uses the Euler-Rodrigues formula to return the the rotation axis and the rotation angle of a given spin 1/2 propagator. The operator has to be 2x2.

    Parameters:
    op : Propagator, must be a 2x2 qutip.Qobj

    Optional parameters:
    tol : Tolerance for the consistency check
    """
    # Convert to SU(2)
    method = backends.get_backend(op)
    op_data = getattr(getattr(backends,method), 'data')(op)
    op = op / _np.sqrt(_np.linalg.det(op_data))
    op_data = getattr(getattr(backends,method), 'data')(op)

    # Check whether this is an appropriate rotation matrix
    if _np.absolute(_np.real(op_data[0,0])-_np.real(op_data[1,1])) > tol:
        raise ValueError("x0 is not consistent, either " + str(_np.real(op_data[0,0])) + " or " + str(_np.real(op_data[1,1])))
    if _np.absolute(_np.imag(op_data[0,1])-_np.imag(op_data[1,0])) > tol:
        raise ValueError("x1 is not consistent, etiher " + str(_np.imag(op_data[0,1])) + " or " + str(_np.imag(op_data[1,0])))
    if _np.absolute(_np.real(op_data[0,1])+_np.real(op_data[1,0])) > tol:
        raise ValueError("x2 is not consistent, either " + str(_np.real(op_data[0,1])) + " or " + str(-_np.real(op_data[1,0])))
    if _np.absolute(_np.imag(op_data[0,0])+_np.imag(op_data[1,1])) > tol:
        raise ValueError("x3 is not consistent, either " + str(_np.imag(op_data[0,0])) + " or " + str(-_np.imag(op_data[1,1])))
    
    # Extract Euler parameters
    x0 = _np.real(op_data[0,0])
    x1 = -_np.imag(op_data[0,1])
    x2 = -_np.real(op_data[0,1])
    x3 = -_np.imag(op_data[0,0])

    # Do the mapping    
    alpha = _np.arccos(x0)*2
    kx = x1/_np.sin(alpha/2)
    ky = x2/_np.sin(alpha/2)
    kz = x3/_np.sin(alpha/2)

    # Constrain to 0..pi rotations
    if alpha >= _np.pi:
        alpha = 2*_np.pi - alpha
        kx = -1*kx
        ky = -1*ky
        kz = -1*kz
    
    return alpha, _np.array([kx,ky,kz])


#############################
# Coordinate Transformation #
#############################

def rotate_operator(system, dm, *args):
    """  Global system rotation, i.e. spin space coordinate transformation of a system operator.  
    The relative orientation of the two coordinate systems can be specified in various formalisms.

    :param system: An instance of the System class.
    :param dm: The operator that will be rotated (Note: if a ket is supplied, it is automatically transformed to a dm and a dm is also returned).
    :param *args: For qutip, numpy and sparse backends this can be a scipy.spatial.transform.Rotation object, for sympy a Sympy Quaternion. For all backends, the input can also be the euler angles (zxy convention) of the rotation.  

    :returns: The density matrix after the coordinate transformation.
    
    """
    # If a ket is provided, make a dm.
    isket_fun = getattr(getattr(backends, system.method), 'isket')
    ket2dm_fun = getattr(getattr(backends, system.method), 'ket2dm')
    dm2ket_fun = getattr(getattr(backends, system.method), 'dm2ket')
    if isket_fun(dm):
        dm = ket2dm_fun(dm)
        ketflag = True
    else:
        ketflag = False
    # Connstruct the global x, y and z operators for spin rotations.
    xop = 0*system.id
    yop = 0*system.id
    zop = 0*system.id
    for level in system.system:
        if level["val"] > 0:
            xop += getattr(system, level["name"]+"x")
            yop += getattr(system, level["name"]+"y")
            zop += getattr(system, level["name"]+"z")
    # Exctract the euler angles from the input.
    if len(args) == 3:
        # Euler angles provided as three separate args.
        euler = [args[0], args[1], args[2]]
    elif (len(args) == 1 and _np.shape(args[0]) != ()):
        # Euler angles provided as single arg.
        euler = args[0]
    elif len(args) == 1:
        if backends.get_backend(args[0]) != "sympy":
            euler = backends.get_calcmethod("toeuler", "numeric")(args[0], "zxy")
        else:
            euler = backends.get_calcmethod("toeuler", "symbolic")(args[0], "zxy")
    # Rotate the state.
    dm = rot(zop, euler[0], dm)
    dm = rot(xop, euler[1], dm)
    dm = rot(yop, euler[2], dm)
    # Return dm or ket.
    if ketflag:
        return dm2ket_fun(dm)
    else:
        return dm

from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

def _precompute_interpolators(carr, dt, interpolation_kind='linear'):
    """
    Prepare interpolated functions or arrays for time-dependent coefficients.

    Parameters:
        carr (list): List of time-dependent coefficient arrays or callables.
        dt (float): Time step for fixed-point evaluation.
        interpolation_kind (str): Kind of interpolation ('linear', 'cubic', etc.).

    Returns:
        list: List of callable functions (either interpolators or original callables).
    """
    interpolators = []
    N = len(carr[0])
    t = _np.arange(0, N * dt, dt)
    for c in carr:
        interpolators.append(interp1d(t, c, kind=interpolation_kind, bounds_error=False, fill_value="extrapolate"))
    return interpolators

def _compute_hamiltonian(t, H0, H1, carr1):
    """
    Efficiently compute the time-dependent Hamiltonian H(t).

    Parameters:
        t (float): Current time.
        H0 (ndarray): Static Hamiltonian.
        H1 (ndarray): Stack of time-dependent Hamiltonians (3D array).
        carr1 (list of callables): Time-dependent coefficients.

    Returns:
        ndarray: The Hamiltonian H(t).
    """
    if H1 is None:
        return H0
    if len(H1) == 0:
        return H0
    coefficients = _np.array([carr(t) for carr in carr1])  # Evaluate all coefficients
    weighted_sum = _np.tensordot(coefficients, H1, axes=1)  # Weighted sum of H1

    if H0 is not None:
        weighted_sum += H0
    return weighted_sum

def _compute_lindblad_terms(t, rho, c_ops, c_ops2, carr2):
    """
    Efficiently compute Lindblad terms for the equation.

    Parameters:
        t (float): Current time.
        rho (ndarray): Current density matrix.
        c_ops (ndarray): List of static collpase operators.
        c_ops2 (ndarray): List of time-dependent collapse operators.
        carr2 (list of callables): Time-dependent coefficients for c_ops2.

    Returns:
        ndarray: Sum of Lindblad terms.
    """
    
    if c_ops2 is None:
        L_total = c_ops
    elif len(c_ops2) == 0:
        L_total = c_ops
    else:
        coefficients = _np.array([carr(t) for carr in carr2])
        L_time_dependent = _np.tensordot(coefficients, c_ops2, axes=1)
        L_total = _np.concatenate((c_ops, L_time_dependent), axis=0)

    # Compute Lindblad terms in batch
    L_rho_Ldag = _np.einsum('aij,jk,akl->ail', L_total, rho, _np.conj(L_total))  # L_k @ rho @ L_k^†
    Ldag_L = _np.einsum('aij,ajk->aik', _np.conj(L_total), L_total)  # L_k^† @ L_k
    anti_commutator = _np.einsum('aij,jk->aik', Ldag_L, rho) + _np.einsum('jk,aik->aij', rho, Ldag_L)
    return _np.sum(L_rho_Ldag - 0.5 * anti_commutator, axis=0)

def lindblad_rhs(t, rho_flat, H0, H1, carr1, c_ops, c_ops2, carr2):
    """
    Compute the right-hand side of the Lindblad equation.

    Parameters:
        t (float): Current time.
        rho_flat (ndarray): Flattened density matrix (1D array).
        H0, H1, carr1, c_ops, c_ops2, carr2: Hamiltonian and collapse operator parameters.

    Returns:
        ndarray: Flattened time derivative of the density matrix.
    """
    n = int(_np.sqrt(len(rho_flat)))
    rho = rho_flat.reshape((n, n))

    # Compute H(t)
    H_t = _compute_hamiltonian(t, H0, H1, carr1)

    # Compute -i[H(t), rho]
    L = -1j * (H_t @ rho - rho @ H_t)

    if not (len(c_ops) == 0 and len(c_ops2) == 0):
       # Compute Lindblad terms
       L += _compute_lindblad_terms(t, rho, c_ops, c_ops2, carr2)
    
    return L.flatten()

def solve_lindblad(rho0, dt, H0=None, H1=None,carr1=None,c_ops=None,c_ops2=None,carr2=None):
    """
    Solve the Lindblad equation

    Parameters:
        rho0 (ndarray): Initial density matrix.
        dt (float): Time step for fixed-point evaluation or end-time when list of functions is provided.
        H0, H1, carr1, c_ops, c_ops2, carr2: Parameters for the Lindblad RHS.
        t_array (ndarray, optional): Array of time points corresponding to fixed-point coefficient arrays.

    Returns:
        ndarray: Density matrix at final time t_eval[-1].
    """

    if carr1 is None:
        carr1 = []
    if carr2 is None:
        carr2 = []
    if c_ops is None:
        c_ops = []
    if c_ops2 is None:
        c_ops2 = []
    
    if H1 is None:
        H1 = []
    
    tend1 = dt
    if len(carr1) > 0:
        # check if carr1[0] is a callable
        if not callable(carr1[0]):
            tend1 = len(carr1[0]) * dt
            carr1 = _precompute_interpolators(carr1, dt)
    if len(carr2) > 0:
        # check if carr2[0] is a callable
        if not callable(carr2[0]):
            carr2 = _precompute_interpolators(carr2, dt)
            tend2 = len(carr2[0]) * dt
            assert tend1 == tend2, "Time-dependent coefficients must have the same length."
    
    t_span = (0, tend1)


    # Convert inputs for batch processing

    rho0_flat = rho0.flatten()
    
    result = solve_ivp(
        lindblad_rhs, t_span, rho0_flat, t_eval=None, vectorized=False,
        args=(H0, H1, carr1, c_ops, c_ops2, carr2),
        method='RK45'
    )
    return result.y[:, -1].reshape(rho0.shape)  # Return final density matrix


#############################
# Lab Frame Evolution       #
#############################

def rabi2b1(w_rabi,S, gyromagnetic_ratio,ms=None):
    """Convert an angular rabi frequency to the corresponding B1 field under the rotating wave approximation ms=None means highest ms"""
    # https://doi.org/10.1016/0009-2614(90)85493-V
    # S = (spin_multiplicity - 1)/2
    # w_rabi = y*b1/2 * sqrt((S+1)*S-ms*(ms-1))
    if ms == None:
        ms = S
    fact = _np.sqrt(S*(S+1)-ms*(ms-1))
    #print(fact)
    return _np.abs((w_rabi/fact*2)/gyromagnetic_ratio)
    
def square_pulse(H0,Hrf,f,phase,amplitude,dur,*rho,pts_per_cycle=100,wallclock='global'):
    # Handle the wallclock
    if wallclock == 'global':
        clock = globalclock
    else:
        clock = wallclock

    if clock is not None:
        phase = phase + clock.phase(f,is_angular=False)
        
    
    # Add wallclock phase
    
    U = _square_propagator_cpu(H0,Hrf,f,phase,amplitude,dur,pts_per_cycle)

    if clock is not None:
        clock.inc(dur)

    # Return the state if rho is given
    if len(rho) == 0:
        return U
    elif len(rho) == 1:
        rho = rho[0]
        # check if rho is a ket or a density matrix
        isket = (rho.shape[1] == 1)
        if isket:
            return U*rho
        else:
            return U*rho*U.dag()
    else:
        raise ValueError('More than one rho was given. This is not supported.')
    
def _square_propagator_cpu(H0,Hrf,f,phase,amplitude,dur,pts_per_cycle=100):
    num_full_osc = int(_np.floor(dur*f))
    dt = 1/(f*pts_per_cycle)
    cycle = _np.linspace(0,1/f-dt,pts_per_cycle)
    cycle = amplitude*_np.sin(2*_np.pi*f*cycle+phase)
    Ucycle = prop(H0,dt,H1=[Hrf],carr1=cycle,wallclock=None)
    
    Ucycle = Ucycle**num_full_osc
    
    # Now add the remaining part of non-complete oscillation
    rf_remain = _np.arange(num_full_osc/f, dur, dt) 
    
    if len(rf_remain) == 0:
        return Ucycle
    
    cycle = amplitude*_np.sin(2*_np.pi*f*rf_remain+phase)
    Uremain = prop(H0,dt,H1=[Hrf],carr1=cycle,wallclock=None)
    return Uremain*Ucycle

def cos2_pulse(H0,Hrf,f,phase,amplitude,dur,*rho,pts_per_cycle=100,wallclock='global'):
    #
    raise NotImplementedError('This function is not yet implemented.')

# Legacy arb_pulse
arb_pulse = prop