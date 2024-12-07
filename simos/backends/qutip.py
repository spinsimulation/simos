import qutip as _qu 
import numpy as _np
import scipy.linalg as _la
from typing import Union
import warnings
import types
from .numpy import differentiation_matrix as _differentiation_matrix

###############################################################
# Utility functions to convert data types and change elements:#
###############################################################

def tidyup(arg, dims=None, atol = None):
    # Type conversion with dims setting for non-quQobj.
    if not isinstance(arg, _qu.Qobj):
        arg = _qu.Qobj(arg, dims=dims)
    # Otherwise set dims, but do not overwrite existing with None.
    else:
        if dims is not None:
            arg.dims = dims
    return arg

def data(arg : _qu.qobj.Qobj):
    # Return the data of the quantum object as a numpy array.
    return arg.full()

def changeelement(matrix : _qu.Qobj, ind1 : int, ind2 :int, newvalue : Union[int, float, complex, list, tuple, _np.ndarray, _qu.Qobj]): 
    # Transform to numpy ndarray.
    dims = matrix.dims
    matrix = matrix.full()
    # If the values to be filled in are also a qutip object, also transform those.
    if isinstance(newvalue, _qu.Qobj):
        newvalue = newvalue.full()
    # Single value case.
    if len(_np.shape(newvalue)) == 0:
        matrix[ind1, ind2] = newvalue
    # Matrix case.
    else:
        newvalue = _np.array(newvalue)
        if len(newvalue.shape) > 2:
            raise ValueError("Changelement only accepts scalars, arrays or 2D matrices.")
        if len(newvalue.shape) == 1:
            newvalue = _np.array([newvalue])
        matrix[ind1:ind1+newvalue.shape[0], ind2:ind2+newvalue.shape[1]] = newvalue
    return tidyup(matrix, dims = dims)


###############################################################
#################     Math / Generic        ###################
###############################################################

def expect(oper,state):
    return _qu.expect(oper,state)

def tensor(operators):
    # Function computes the tensor product of the operators in the list operators.
    return _qu.composite(*operators)

def directsum(*args):
    if len(args) > 1:
        operators  = [*args]
    else:
        operators = args[0]
    dimtot = sum([i if isinstance(i, int) else _np.shape(i)[0]  for i in operators])
    dimtot2 = sum([i if isinstance(i, int) else _np.shape(i)[1]  for i in operators])
    M = _np.zeros((dimtot,dimtot2), dtype = complex)
    pre = 0
    pre2 = 0
    for operator in operators:
        if isinstance(operator, int):
            pre += operator
            pre2 += operator
        else:
            M[pre:pre+_np.shape(operator)[0], pre2:pre2+_np.shape(operator)[1]] = operator.full()
            pre += _np.shape(operator)[0]
            pre2 += _np.shape(operator)[1]
    M = _qu.Qobj(M)   
    return M 

def ddrop(operator, idx):
    operator = operator.full()
    for axis in (0,1):
        operator = _np.delete(operator, idx, axis = axis)
    operator = _qu.Qobj(operator, dims = [[_np.shape(operator)[0]], [_np.shape(operator)[1]]])
    return operator

def dpermute(operator, neworder):
    if len(operator.shape) != 2 or operator.shape[0] != operator.shape[1]:
        raise ValueError("Dpermute can only be applied to 2-dimensional, square quantum objects.")
    operator = operator.full()
    operator =  operator[:, neworder]
    operator = operator[neworder,:]
    return _qu.Qobj(operator, dims = [[_np.shape(operator)[0]], [_np.shape(operator)[1]]])

def block_diagonal(L_list):
    L_list = [Li.full() for Li in L_list]
    Lnew = _la.block_diag(*L_list)
    dims = [[len(L_list),L_list[0].shape[0]],[len(L_list),L_list[0].shape[1]]]
    obj =  tidyup(Lnew,dims=dims)
    return obj

###############################################################
#################    Hilbert Space          ###################
###############################################################

# Kets, bras, operators and their conversion. 

def isket(obj):
    return obj.isket

def isbra(obj):
    return obj.isbra

def isoper(obj):
    return obj.isoper

def ket2dm(ket):
    return _qu.ket2dm(ket)

def dm2ket(dm, tol=1e-6): # NOTE: tol is smaller in qutip than in other backends because qutip .unit() function is not good enough to handle smaller tols 
    if dm.tr() < 1-tol or dm.tr() > 1+tol:
        raise ValueError('Density matrix is not normalized.')
    if (dm**2).tr() < (1-tol):
        raise ValueError('Density matrix is not a pure state.')
    if _np.max(dm.diag()) >= 1-tol:
        return _qu.Qobj(dm.diag(), dims = [dm.dims[0], [1 for i in dm.dims[1]]])
    else:
        # Find the decmposition into the pure state that produces the density matrix.
        evals, evecs = dm.eigenstates()
        idx = _np.argmax(evals)
        return _qu.Qobj(evecs[idx], dims = [dm.dims[0], [1 for i in dm.dims[1]]])

# Part 2: Construction of standard kets, bras and operators:

def ket(m, n):
    return _qu.basis(m,n)

def bra(m, n):
    return _qu.basis(m,n).trans()

def jmat(*args):
    return _qu.jmat(*args)

def identity(N, dims = None):
    out = _qu.qeye(N)
    if dims is None:
        dims = list(out.shape)
    out.dims = dims
    return out

def diags(v : Union[_qu.Qobj, list, tuple, _np.ndarray], k : int = 0, dims = None):
    # Return the kth diagonal if a 2D quantum object is provided.
    if len(_np.shape(v)) == 2:
        if not isinstance(v, _qu.Qobj):
            raise ValueError("Diagonal is only extracted from data type qutip.Qobj. For extraction from  list, tuple, array etc. directly utilize the diags function of the numpy backend.")
        if dims is not None:
            warnings.warn("Dims keyword argument will be ignored.")
        v = v.full()
        out = _np.diag(v, k)
        return out 
    # Otherwise a quantum object is consructed from the input diagonal.
    else:
        out = _qu.qdiags(v, k)
        if dims is not None:
            out.dims = dims
        return out


###############################################################
#################    Liouville Space         ##################
###############################################################

def dm2vec(op):
    return _qu.operator_to_vector(op)

def vec2dm(vec):
    # Note that adjusment to standard vector_to_operators from qutip was necessary 
    # to enable FokkerPlanck functionality.
    if vec.superrep != 'super':
        s = max(vec.shape)
        s = int(_np.sqrt(s))
        dims = [[[s], [s]], [1]]
        vec = _qu.Qobj(vec.full(), superrep='super', dims=dims)
    return _qu.vector_to_operator(vec)

def issuper(obj):
    return obj.issuper

def spre(obj):
    return _qu.spre(obj)

def spost(obj):
    return _qu.spost(obj)

def liouville(H):
    return _qu.liouvillian(H)

def liouvillian(*args,**kwargs):
    return _qu.liouvillian(*args,**kwargs)

def applySuperoperator(S, rho):
    rhovec = _qu.operator_to_vector(rho)
    rhovec = S*rhovec
    return _qu.vector_to_operator(rhovec)


###############################################################
#################    Fokker Planck Space      #################
###############################################################
def differentiation_matrix(N,m,method='optimal',boundary='periodic'):
    return _qu.Qobj(_differentiation_matrix(N,m,method=method,boundary=boundary))

def concatenate(*states, dims=None):
    if len(states) == 1:
        states = states[0]
    #qutip_type = states[0].type
    return _qu.Qobj(_np.concatenate([state.full() for state in states], axis=0),dims=dims)

def split(state, n:int, new_dims=None, is_hilbert=False):
    d = state.shape[0]
    if d % n != 0:
        raise ValueError("State dimension must be divisible by n.")
    L = d//n
    state = state.full()
    #return [_qu.Qobj(state[i*L:(i+1)*L],dims=new_dims) for i in range(n)]
    segment = [state[i*L:(i+1)*L] for i in range(n)]
    if not is_hilbert:
        # make segments square matrices
        L2 = int(_np.sqrt(L))
        segment = [_np.reshape(segment[i],(L2,L2)).T for i in range(n)]
    return [_qu.Qobj(segment[i],dims=new_dims) for i in range(n)]

def fp2dm(v,is_hilbert,dims=None,weights=None):
    dof = v.dims[0][0]

    rhos = split(v, dof,is_hilbert=is_hilbert)

    if weights is None:
        weights = _np.ones(dof)/dof
    
    if isinstance(weights, types.GeneratorType):
        weights = list(weights)

    if len(weights) != dof:
        raise ValueError("Number of weights must match number of operators.")    

    if is_hilbert:
        rhos = [ket2dm(rhos[i])*weights[i] for i in range(len(rhos))]
    else:
        rhos = [(rhos[i])*weights[i] for i in range(len(rhos))]

    rhoAVG = sum(rhos).unit()    
    rhoAVG = tidyup(rhoAVG, dims = dims)
    return rhoAVG

def dm2fp(rho,dof,try_ket=True):
    if try_ket:
        try:
            rho = dm2ket(rho)
        except ValueError:
            pass
    if not isket(rho):
        rho = _qu.operator_to_vector(rho)
    rhoF = concatenate([rho]*dof,dims = [[dof,rho.shape[0]],[1,1]])
    return rhoF

def ket2fp(rho,dof):
    if not isket(rho):
        raise ValueError("Input state must be a ket.")
    rhoF = concatenate([rho]*dof,dims = [[dof,rho.shape[0]],[1,1]])
    return rhoF


# def applyFokkeroperator(F,rho,weights=None):
#     if weights is None:
#         weights = _np.ones(len(F))/len(F)
#     if len(weights) != len(F):
#         raise ValueError("Number of weights must match number of operators.")
#     F_shape = F[0].shape
#     F_dims = F[0].dims
#     rho_shape = rho.shape
#     rho_dims = rho.dims

#     # Decide whether we are dealing with a series of
#     # Hilbert space propagators
#     # or a series of Liouville space propagators.
#     # or a fokker-Hilbert space propagator.
#     # or a fokker-Liouville space propagator.
#     if F_shape == rho_shape:
#         print("Hilbert space series")
#         # Hilbert space propagators.
#         if len(weights) != len(F):
#             raise ValueError("Number of weights must match number of operators.")
#         if isket(rho):
#             rho = ket2dm(rho)
#         rhos = [F[i]*rho*weights[i]*F[i].dag() for i in range(len(F))]
#     elif F_shape[0] == rho_shape[0]*rho_shape[1]:
#         print("Liouville space series")
#         # Liouville space propagators.
#         if len(weights) != len(F):
#             raise ValueError("Number of weights must match number of operators.")
#         if isket(rho):
#             rho = ket2dm(rho)
#         rho = _qu.operator_to_vector(rho)
#         rhos = [_qu.vector_to_operator(F[i]*rho*weights[i]) for i in range(len(F))]
#     else:
#         # Fokker-space operator
#         dof = F_dims[0][0]
#         if len(F) != 1:
#             raise ValueError("Fokker space propagator must be a single operator.")
        
#         if len(weights) != dof:
#             raise ValueError("Number of weights must match number of operators.")
        
#         if F_dims[0][1] == rho.shape[0]: # we are in the Hilbert space case
#             if not isket(rho):
#                 rho = dm2ket(rho)
#             rhoF = concatenate([rho]*dof,dims = [[dof,rho.shape[0]],[1,1]])
#             rhoF = F[0]*rhoF
#             rhos = split(rhoF, dof)
#             rhos = [ket2dm(rhos[i])*weights[i] for i in range(len(rhos))]
#         else: # we are in the Liouville space case
#             rho = _qu.operator_to_vector(rho)
#             rhoF = concatenate([rho]*dof,dims = [[dof,rho.shape[0]],[1,1]])
#             rhoF = F[0]*rhoF
#             rhos = split(rhoF, dof)
#             rhos = [_qu.vector_to_operator(rhos[i])*weights[i] for i in range(len(rhos))]
    
#     rhoAVG = sum(rhos).unit()    
#     rhoAVG = tidyup(rhoAVG, dims = rho_dims)
#     return rhoAVG