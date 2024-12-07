import numpy as _np
from .qmatrixmethods import data
from . import backends

##############################################################
# Helper-Functions for Spherical Tensor Notation
# Note: SimOS core routines do currently not utilize spherical
# tensors extensively.
##############################################################

### Spatial spherical tensors (rank 0, 1,2) <-> 3x3 matrices ###

def mat2spher(mat): 
    """ Transforms a 3x3 matrix into a spherical tensor (ranks 0,1,2).
    
    :param mat: 3x3 matrix
    :returns dict: Spherical tensor representation of input matrix as a dictionary with keys 0, 1, 2 for the tensor ranks."""
    
    # Ensure that we work with data. If a quantum object is provided, extract data.
    if hasattr(mat, "dims"):
        mat = data(mat)
    # Probe for symbolic vs numeric
    if backends.get_backend(mat) != 'sympy':
        calcmode = "numeric"
    else:
        calcmode = "symbolic"  
    trace_fun = backends.get_calcmethod("trace", calcmode)
    multelem_fun = backends.get_calcmethod("multiply", calcmode)
    pow_fun = backends.get_calcmethod("pow", calcmode )
    I = backends.get_calcmethod("I", calcmode)
    sqrt_fun = backends.get_calcmethod("sqrt", calcmode)
    # Initialise the dictionary for the spherical tensor notation.
    sphten = {}
    # Rank 0 component.
    sphten[0] = {}
    sphten[0][0] = -1* pow_fun(sqrt_fun(3), -1) * trace_fun(mat)
    # Rank 1 components.
    sphten[1] = {}
    sphten[1][-1] = pow_fun(-2, -1)*(mat[2,0]-mat[0,2]+I*(mat[2,1]-mat[1,2]))
    sphten[1][0] =  pow_fun(-sqrt_fun(2), -1)*I*(mat[0,1]-mat[1,0])
    sphten[1][1] =  pow_fun(-2, -1)*(mat[2,0]-mat[0,2]-I*(mat[2,1]-mat[1,2]))
    # Rank 2 components.
    sphten[2]= {}
    sphten[2][-2]= pow_fun(2, -1)*(mat[0,0]-mat[1,1]+I*(mat[0,1]+mat[1,0]))
    sphten[2][-1]= pow_fun(-2, -1)*(mat[0,2]+mat[2,0]+I*(mat[1,2]+mat[2,1]))
    sphten[2][0]=  pow_fun(sqrt_fun(6), -1)*(2*mat[2,2]-mat[0,0]-mat[1,1])
    sphten[2][1]=  pow_fun(2, -1)*(mat[0,2]+mat[2,0]-I*(mat[1,2]+mat[2,1]))
    sphten[2][2]=  pow_fun(2, -1)*(mat[0,0]-mat[1,1]-I*(mat[0,1]+mat[1,0]))
    return sphten

def spher2mat(spher, selrank = [0, 1, 2]):
    """ Transforms a spherical tensor (max. rank 2) into a 3x3 matrix.
    
    :param spher: Dictionary representing a spherical tensor.
    :returns mat: 3x3 matrix.
    """
    # Probe for symbolic vs numeric
    if backends.get_backend(spher) != 'sympy':
        calcmode = "numeric"
    else:
        calcmode = "symbolic"  
    array_fun = backends.get_calcmethod("array", calcmode)
    I = backends.get_calcmethod("I", calcmode)
    sqrt_fun = backends.get_calcmethod("sqrt", calcmode)
    pow_fun = backends.get_calcmethod("pow", calcmode )
    multelem_fun = backends.get_calcmethod("multiply", calcmode)
    # Initialise matrix.
    mat = array_fun([[0,0,0], [0,0,0], [0,0,0]])
    # Build matrices for all rank contributions (0,1,2).
    for l in spher.keys():
        # Rank 0 component.
        if l == 0 and l in selrank:
            mat = mat - pow_fun(sqrt_fun(3), -1) * multelem_fun(spher[0][0],array_fun([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
        # Rank 1 components.
        if l == 1 and l in selrank:
            mat = mat - multelem_fun(spher[1][-1]*pow_fun(2, -1) , array_fun([[0, 0, -1], [0, 0, I], [1, -I, 0]]))
            mat = mat + multelem_fun(spher[1][0]*pow_fun(sqrt_fun(2), -1),  array_fun([[0, I, 0], [-I, 0, 0], [0, 0,0]]))
            mat = mat - multelem_fun(spher[1][1]*pow_fun(2, -1),  array_fun([[0, 0, -1], [0, 0, -I], [1, I, 0]]))
        # Rank 2 components.
        if l == 2 and l in selrank:
            mat = mat + multelem_fun(spher[2][-2]*pow_fun(2,-1),array_fun([[1, -I, 0], [-I, -1, 0], [0, 0, 0]]))
            mat = mat + multelem_fun(spher[2][-1]*pow_fun(-2,-1),array_fun([[0, 0, 1], [ 0, 0, -I], [1, -I, 0]]))
            mat = mat + multelem_fun(spher[2][0]*pow_fun(sqrt_fun(6), -1), array_fun([[-1, 0, 0], [0, -1, 0], [0, 0, 2]]))            
            mat = mat + multelem_fun(spher[2][1]*pow_fun(2, -1), array_fun([[0, 0, 1], [0, 0, I], [1, I, 0]]))
            mat = mat + multelem_fun(spher[2][2]*pow_fun(2,-1), array_fun([[1, I, 0], [I, -1, 0], [0, 0, 0]]))
    return mat


### Spin spherical tensor operators ###


def spherspin(spinsystem, spinname:str): 
    """Spherical tensor operators for a spin of an existing system.

    :param spinsystem: An instance of the System class.
    :param str spinname: Name of a spin of the system.
    :returs dict: Spherical tensor operators of the spin,  a dictionary with keys for the tensor ranks."""

    p = getattr(spinsystem, spinname+"p")
    val = [i for i in p.keys()][-1]
    N = int(2*val+1) # multiplicity 
    T = {}  # preallocate storage
    id = getattr(spinsystem, spinname+"id")
    if spinsystem.method == "sympy":
        calcmode = "symbolic"
    else:
        calcmode = "numeric"
    pow_fun = backends.get_calcmethod("pow", calcmode)
    sqrt_fun = backends.get_calcmethod("sqrt", calcmode)
    if N == 1:
        T[0] = {} 
        T[0][0] =  id
    else:
        Lp = getattr(spinsystem, spinname+"plus")
        Lm = getattr(spinsystem, spinname+"minus")
        for k in range(N): # loop over ranks (k)
            T[k] = {} # every rank again holds a dictionary 
            if k == 0: # if rank 0, return a unit matrix 
                T[k][0] =  id
            else:
                T[k][k] =  pow_fun(-1, k)* pow_fun(2, -k*pow_fun(2,-1))* Lp**k  #((-1)**k)*(2**(-k/2))*Lp**k  get the top state  
                for q in _np.arange(k-1, -k-1, -1): # apply sequential lowering using Racah's commutation rule
                    pre =  pow_fun( sqrt_fun((k+(q+1))*(k-(q+1)+1)), -1) #1/_np.sqrt((k+(q+1))*(k-(q+1)+1)) 
                    T[k][q] = pre * (Lm*T[k][q+1] - T[k][q+1]*Lm)
    return T 


import itertools
def spherbasis(system): 
    """Complete basis of spherical tensor operators for all spins of a system.
    
    :params system: An instance of the system class.
    :returns dict: Complete spherical tensor operator basis of the input system, as a dictionary with keys for the individual spins and combinations thereof.""" 
    # Retrieve backend specific methods. 
    data_fun = getattr(getattr(backends, system.method), 'data')
    cg_fun = backends.get_calcmethod("cg", "numeric")
    # Define coupling of spherical tensor operators.                  
    def couplesphten(T1, T2):
        # allocate storage 
        cst = []
        # drop zero coeff, these wont couple  #T1.pop(0) #T2.pop(0)
        # get ranks
        for l1 in T1.keys():
            for l2 in T2.keys():
                Lvals = _np.unique(_np.arange(_np.abs(l1-l2), l1+l2+1))
                for L in Lvals:
                    pstdict = {}
                    for M in  _np.arange(-L, L+1, 1):
                        pst = 0*system.id
                        for m1 in _np.arange(-l1, l1+1, 1):
                            for m2 in _np.arange(-l2, l2+1, 1):
                                pst +=  cg_fun(int(l1), int(m1), int(l2), int(m2), int(L), int(M)) *  T1[l1][m1]*T2[l2][m2]
                        pstdict[M] = pst
                    cst.append({L: pstdict})
        return cst
    # Get all spins of the system.
    spin_idx = _np.where([i["val"] > 0 for i in system.system])[0]
    spin_names = [system.system[i]["name"] for i in spin_idx]
    spin_ids = [getattr(system, spin_name+"id") for spin_name in spin_names]
    spin_pos = [_np.where(data_fun(id) > 10e-10)[0] for id in spin_ids]
    # Sort into subsystems.
    indices = _np.argsort(_np.array([len(pos) for pos in spin_pos]))
    subsystem_pos = []
    subsystem_names = []
    for ind in indices:
        added = False
        for pos, name in zip(subsystem_pos, subsystem_names):
            if all(i in spin_pos[ind] for i in pos):
                name += spin_names[ind] # this spin is part of this spinsystem
                added = True
        if not added:
            subsystem_pos.append(spin_pos[ind])
            subsystem_names.append([spin_names[ind]])
    Tsingle  =  {}
    #allT = []
    # Get single spin operators, strap 00 component. 
    for spin_name in spin_names:
        T = spherspin(system, spin_name)
        T.pop(0)
        Tsingle[spin_name] = [T]
        #allT += [{i: T[i]} for i in T.keys()] # flatten dict 
    # Multi spin operators , separately for each subsystem. 
    for subsystem in subsystem_names: # loop over subsystem_names
        nspins = len(subsystem) # how many spins are in subsystem
        sorted_subsystem = sorted(subsystem) 
        for n in range(nspins):
            loopnames = list(itertools.combinations(sorted_subsystem, n+2))
            for loopname in loopnames:
                t1 = Tsingle["".join([i for i in loopname[0:-1]])]
                t2 = Tsingle["".join([i for i in loopname[-1]])]
                cst = []
                for i in t1:
                    for j in t2:
                        cst += couplesphten(i, j) 
                Tsingle["".join([i for i in loopname[:]])] = cst
                #print("added", "".join([i for i in loopname[:]]))
    # Add system id.
    Tsingle["id"] = [{0: { 0: system.id}}]
    # Calc length.
    L = 0
    for Tname in Tsingle.keys():
        Tspin = Tsingle[Tname] # das ist eine liste 
        keys = [[i for i in j.keys()] for j in Tspin]
        Lsub = [[2*i+1 for i in j]  for j in keys]
        Ltot = _np.sum(_np.array(Lsub))
        L += Ltot
    return Tsingle