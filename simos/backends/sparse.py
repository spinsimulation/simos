import scipy.sparse as _sp
import numpy as _np
import functools
from typing import Union
import warnings
from .. trivial import flatten
from .numpy import differentiation_matrix as _differentiation_matrix


class Qcsc_matrix(_sp.csc_matrix): 
    """ Scipy sparse analogue of a qutip object."""
    
    def __init__(self, arg1, dims  = None, shape = None, dtype = None,copy=False):
        _sp.csc_matrix.__init__(self, arg1, shape = shape, dtype = dtype,copy=copy)
        if dims is None:
            self.dims = [[self.shape[0]], [self.shape[1]]]
        else:
            self.set_dims(dims)

    # Change Class Methods
    def __mul__(self,other):
        """
        MULTIPLICATION with QNdarray on LEFT [ ex. Qobj*4 ]
        """
        if isinstance(other, (Qcsc_matrix,_sp.coo_matrix,_sp.csr_matrix,_sp.csc_matrix)):
            return tidyup(self@other, dims = [self.dims[0], other.dims[1]])
        elif isinstance(other, list):
            return [tidyup(_sp.csc_matrix(self)*i) for i in other]
        elif isinstance(other, (int, float, complex, _np.integer, _np.floating, _np.complexfloating)):
            return tidyup(_sp.csc_matrix(self)*other, dims = self.dims)
        else:
            raise TypeError("Incompatible object for multiplication")

    def __rmul__(self,other):
        """
        MULTIPLICATION with QNdarray on RIGHT [ ex. 4*Qobj ]
        """
        if isinstance(other, Qcsc_matrix):
            return tidyup(self@other, dims = self.dims)
        elif isinstance(other, list):
            return [tidyup(_sp.csc_matrix(self)*i) for i in other]
        elif isinstance(other, (int, float, complex, _np.integer, _np.floating, _np.complexfloating)):
            return tidyup(_sp.csc_matrix(self)*other, dims = self.dims)
        else:
            raise TypeError("Incompatible object for multiplication")
    
    def __add__(self,other):
        """
        ADDITION of a QNdarray [ ex. Qobj+4 ]
        """
        if isinstance(other, (Qcsc_matrix,_sp.coo_matrix,_sp.csr_matrix,_sp.csc_matrix)):
            return tidyup(_sp.csc_matrix(self)+other, dims = self.dims)
        elif isinstance(other, (int, float, complex, _np.integer, _np.floating, _np.complexfloating)):
            return tidyup(_sp.csc_matrix(self)+other, dims = self.dims)
        else:
            raise TypeError("Incompatible object for addition")
        
    def __radd__(self,other):
        """
        ADDITION of a QNdarray [ ex. 4+Qobj ]
        """
        if isinstance(other, (Qcsc_matrix,_sp.coo_matrix,_sp.csr_matrix,_sp.csc_matrix)):
            return tidyup(_sp.csc_matrix(self)+other, dims = self.dims)
        elif isinstance(other, (int, float, complex, _np.integer, _np.floating, _np.complexfloating)):
            return tidyup(_sp.csc_matrix(self)+other, dims = self.dims)
        else:
            raise TypeError("Incompatible object for addition")
        
    def __sub__(self,other):
        """
        SUBTRACTION of a QNdarray [ ex. Qobj-4 ]
        """
        if isinstance(other, (Qcsc_matrix,_sp.coo_matrix,_sp.csr_matrix,_sp.csc_matrix)):
            return tidyup(_sp.csc_matrix(self)-other, dims = self.dims)
        elif isinstance(other, (int, float, complex, _np.integer, _np.floating, _np.complexfloating)):
            return tidyup(_sp.csc_matrix(self)-other, dims = self.dims)
        else:
            raise TypeError("Incompatible object for subtraction")
        
    def __rsub__(self,other):
        """
        SUBTRACTION of a QNdarray [ ex. 4-Qobj ]
        """
        if isinstance(other, (Qcsc_matrix,_sp.coo_matrix,_sp.csr_matrix,_sp.csc_matrix)):
            return tidyup(other-_sp.csc_matrix(self), dims = self.dims)
        elif isinstance(other, (int, float, complex, _np.integer, _np.floating, _np.complexfloating)):
            return tidyup(other-_sp.csc_matrix(self), dims = self.dims)
        else:
            raise TypeError("Incompatible object for subtraction")
    
    def __pow__(self,other):
        """
        EXPONENTIATION of a QNdarray [ ex. Qobj**4 ]
        """
        if isinstance(other, (int, float, complex, _np.integer, _np.floating, _np.complexfloating)):
            return tidyup(_sp.linalg.matrix_power(self, other), dims = self.dims)
        else:
            raise TypeError("Incompatible object for power")

    def set_dims(self, dims : list):
        if len(dims) == 2 and _np.prod(dims[0]) == self.shape[0] and _np.prod(dims[1]) == self.shape[1]:
            self.dims = dims
        else:
            raise ValueError("Dims cannot be set to the specified value.")
        
    # Add Class Methods 
    def dag(self):
        out = self.conjugate().transpose() # native conjugate, transpose functions 
        out = Qcsc_matrix(out, dims = [self.dims[1], self.dims[0]])
        return out

    def conj(self):
        out = self.conjugate() # native conjugate 
        out = Qcsc_matrix(out, dims = self.dims)
        return out

    def trans(self): 
        out = self.transpose() #  native transpose function
        out = Qcsc_matrix(out, dims = [self.dims[1], self.dims[0]])
        return out

    def tr(self, atol = 1e-10): # TEST
        out = self.trace()
        if _np.abs(_np.imag(out)) <= atol: 
            return float(_np.real(out))  # check if real valued , if yes do type conversion
        else: 
            return out

    def diag(self): 
        out = self.diagonal(k=0)
        return out
    
    def unit(self):
        out = (self/_sp.linalg.norm(self))**2 
        out.dims = self.dims
        # nuclear norm
        #norm = ((self.dag()*self).sqrt()).tr()
        #out = self/norm
        return out
    
    def transform(self, U):  
        U = tidyup(U)
        out = U.dot(self.dot(U.transpose()))
        out.dims = self.dims
        return out

    def expm(self, eliminate_zeros = True):
        dims = self.dims
        # Tidyup to eliminate zeros
        self = _sp.linalg.expm(self)
        if eliminate_zeros:
            self.eliminate_zeros()
        return Qcsc_matrix(self,dims=dims)

    def ptrace(self, sel):
        # Currently ptrace is only working for operators. Catch that.
        if not isoper(self):
            raise ValueError("Sparse backend only supports partial trace for operators.")
        # handle type of sel
        if isinstance(sel, int):
            sel = [sel]
        # prepare dimensions 
        dims = self.dims
        N = len(self.dims[0])# wie viele spins sind da
        dimflat = flatten(dims)
        i = _np.arange(0, N, 1) # all axis 1 
        j = i + N # all axis 2 
        # determine over which spins one has to trace -> all spin whose index is not in the selection list (sel is those which are being kept)
        sel.sort() 
        all = list(_np.arange(0,N, 1))
        notsel = [x for x in all if x not in sel] # x ist automatisch sortiert 
        # reshape the array 
        objt = _np.reshape(self.toarray(), dimflat)
        # trace out 
        for ind_A, A in enumerate(notsel): # loop over the components which we want to trace out
            axis1 = i[A]-ind_A
            axis2 = j[A]-2*ind_A
            objt = _np.trace(objt, axis1 = axis1, axis2 = axis2)
        # final reshape 
        dimnew = [ [dims[0][s] for s in sel], [dims[1][s] for s in sel]] # dimension of obj after ptrace 
        objt = _np.reshape(objt, [_np.prod(dimnew[0]), _np.prod(dimnew[1])]) # reshape properly, only necessary if more than one componend has not been traced over 
        return tidyup(objt, dims = dimnew)
    
    def eigenstates(self):
        raise NotImplementedError("Not yet implemented.")
    
###############################################################
# Utility functions to convert data types and change elements:#
###############################################################

def tidyup(arg, dims=None, eliminate_zeros = True):
    if not isinstance(arg, Qcsc_matrix):
        arg = Qcsc_matrix(arg, dims = dims)
    else:
        if dims is not None:
           arg.set_dims(dims)
    if eliminate_zeros:
        arg.eliminate_zeros()
    return arg

def data(arg):
    if isinstance(arg, _np.ndarray):
        return(arg)
    elif isinstance(arg, Qcsc_matrix):
        return arg.toarray()
    else:
        raise ValueError("Invalid input data type.")

def changeelement(matrix, ind1 :int, ind2 :int, newvalue: Union[int, float, complex, list, tuple, _np.ndarray, Qcsc_matrix, _sp.spmatrix, _sp.sparray] ): 
    dims = matrix.dims
    # Transform to lil_matrix to allow for most efficient change of the sparsity structure.
    matrix = matrix.tolil() 
    # Set a single value.
    if len(_np.shape(newvalue)) == 0:
        matrix[ind1, ind2] = newvalue
    # Set an array or a matrix.
    else:
        if len(_np.shape(newvalue)) > 2:
            raise ValueError("Changelelement only accepts scalars, arrays or 2D matrices.")
        if len(_np.shape(newvalue)) == 1:
            newvalue = _np.transpose(_np.array([newvalue]))
        matrix[ind1:ind1+_np.shape(newvalue)[0], ind2:ind2+_np.shape(newvalue)[1]] = newvalue
    matrix= matrix.tocsc()
    return tidyup(matrix, dims = dims)


###############################################################
#################     Math / Generic        ###################
###############################################################

def expect(oper, state):
    if isbra(state):
        state = state.dag()
    if isket(state):
        return (state.dag()*oper*state)
    else:
        return (oper*state).tr()

def tensor(operators): 
    """ Function computes the tensor product of the operators in the list operators."""
    dims = [flatten([i.dims[0] for i in operators]), flatten([i.dims[1] for i in operators]) ]
    out = functools.reduce(_sp.kron, operators)
    return tidyup(out, dims = dims)

def directsum(operator, pre, post):
    dimtot = pre + operator.shape[0] + post
    M = _sp.lil_matrix((dimtot, dimtot), dtype = complex)
    M[pre:pre+operator.shape[0], pre:pre+operator.shape[1]] = operator
    M = Qcsc_matrix(M, dims = None)   
    return M 

def block_diagonal(L_list):
    Lnew = _sp.block_diag(L_list)
    dims = [[len(L_list),L_list[0].shape[0]],[len(L_list),L_list[0].shape[1]]]
    obj =  tidyup(Lnew,dims=dims)
    return obj

###############################################################
#################    Hilbert Space          ###################
###############################################################

# Kets, bras, operators and their conversion. 

def isket(obj):
    return (obj.shape[1] == 1)

def isbra(obj):
    return (obj.shape[0] == 1)

def isoper(obj):
    return (isket(obj) is False and isbra(obj) is False and len(obj.dims[0]) == 1 and len(obj.dims[1]) == 1)

def ket2dm(ket): 
    if isket(ket):
        bra = ket.dag()
        return tidyup(ket.dot(bra), dims = [ket.dims[0], ket.dims[0]])
    elif isbra(ket):
        bra = ket.dag()
        return tidyup(bra.dot(ket), dims = [ket.dims[1], ket.dims[1]])
    else:
        raise ValueError("Invalid input shape.")
    
def dm2ket(dm, tol=1e-10):
    if dm.tr() < 1-tol or dm.tr() > 1+tol:
        raise ValueError('Density matrix is not normalized.')
    if (dm**2).tr() < (1-tol):
        raise ValueError('Density matrix is not a pure state.')
    if _np.max(dm.diag()) >= 1-tol:
        return Qcsc_matrix(_np.transpose([dm.diag()]), dims = [dm.dims[0], [1 for i in dm.dims[1]]])
    else:
        # Find the decmposition into the pure state that produces the density matrix.
        if dm.shape[0] == 2:
            evals, evecs = _np.linalg.eig(dm.toarray())
            idx = _np.argmax(evals)
        else:
            evals, evecs = _sp.linalg.eigs(dm, k=1, which = 'LR')
            idx = 0
            if len(evals) >1:
                evec = evecs[idx]
            else:
                evec = evecs
        return tidyup(evec, dims = [dm.dims[0], [1 for i in dm.dims[1]]])


# Part 2: Construction of standard operators.

def ket(m, n):
    vec = _sp.lil_matrix((m, 1))
    vec[n, 0] = 1
    return Qcsc_matrix(vec)

def bra(m, n):
    vec = _sp.lil_matrix((1, m))
    vec[0, n] = 1
    return Qcsc_matrix(vec)

def identity(N, dims = None):
    """Return identity matrix"""
    out = _sp.identity(N,dtype=complex)
    return Qcsc_matrix(out, dims = dims)

def jmat(j,op_spec):
    """Return the spin operator for a spin-j system.
    j is a non-negative integer or half-integer specifying the spin.
    op_spec is a string specifying the desired operator:
    """
    if (_np.fix(2 * j) != 2 * j) or (j < 0):
        raise TypeError('j must be a non-negative integer or half-integer')
    N = int(2*j+1)

    m = _np.arange(j, -j - 1, -1, dtype=complex)
    pm_data = (_np.sqrt(j * (j + 1.0) - (m + 1.0) * m))[1:]

    if op_spec == '+':
        return tidyup(_np.diag(pm_data, k=1))
    elif op_spec == '-':
        return tidyup(_np.diag(pm_data, k=-1))
    elif op_spec == 'x':
        return tidyup(0.5 * (_np.diag(pm_data, k=1) + _np.diag(pm_data, k=-1)))
    elif op_spec == 'y':
        return tidyup(-0.5 * 1j * (_np.diag(pm_data, k=1)  - _np.diag(pm_data, k=-1)))
    elif op_spec == 'z':
        data = _np.array([j-k for k in range(N)], dtype=complex)
        return tidyup(_np.diag(data))
    else:
        raise TypeError('Invalid operator specification')


def diags(v : Union[Qcsc_matrix, list, tuple, _np.ndarray], k  : int =0, dims = None):
    # Return the kth diagonal if a 2D quantum object is provided.
    if len(_np.shape(v)) == 2:
        if not isinstance(v, Qcsc_matrix):
            raise ValueError("Diagonal is only extracted from data type Qcsc_matrix. For extraction from  list, tuple, array etc. directly utilize the diags function of the numpy backend.")
        if dims is not None:
            warnings.warn("Dims keyword argument will be ignored.")
        out = v.diagonal(k=k).toarray()
        return out  
    # # If a 1D input was provided, build a diagonal matrix.   
    else:
        out = _sp.diags(v,k) 
        return tidyup(out, dims = dims)

###############################################################
#################    Liouville Space         ##################
###############################################################

def operator_to_vector(op):
    # stack the columns of the matrix
    dims = [op.dims, [1]]
    return tidyup(op.T.reshape((-1, 1)),dims=dims)

def vector_to_operator(vec):
    s = max(vec.shape)
    s = int(_np.sqrt(s))
    # unstack the columns of the matrix
    dims = [[s], [s]]
    return tidyup(vec.reshape((s, s), order='F'),dims=dims)

def issuper(obj):
    return( len(obj.dims[0]) == 2 and len(obj.dims[1]) == 2)

def spre(H):
    ident = _sp.eye(H.shape[0])
    pre =  _sp.kron(ident,H)
    pre = Qcsc_matrix(pre,  dims = [[H.dims[0], H.dims[0]], [H.dims[1], H.dims[1]]])
    return pre

def spost(H):
    # Kronecker product of 1 otimes H.T
    ident = _sp.eye(H.shape[0])
    post = _sp.kron(H.T,ident)
    post = Qcsc_matrix(post, dims = [[H.dims[0], H.dims[0]], [H.dims[1], H.dims[1]]])
    return post

def sprepost(H : Qcsc_matrix ,factor=-1):
    """Computes 1 \\otimes H + factor * H.T \\otimes 1"""
    L = spre(H) + spost(H)*factor
    return L

def liouville(H : Qcsc_matrix):
    return -1j*sprepost(H,factor=-1)

def lindbladian(a : Qcsc_matrix ,b=None):
    if b is None:
        b = a
    #bdag = b.dag()
    #adag = a.dag()
    #part1 = Qcsc_matrix(_sp.kron(bdag,a),  dims = [[a.dims[0], a.dims[0]], [a.dims[1], a.dims[1]]])
    #return part1 - sprepost(adag*(b),factor=1)/2
    ad_b = a.dag() * b
    L = spre(a) * spost(b.dag()) - spre(ad_b)/2 - spost(ad_b)/2
    return L

def liouvillian(H,cops=None):
    if isinstance(cops,list) or isinstance(cops,tuple):
        if len(cops)==0:
            cops = None
    if cops is None:
        L =  liouville(H)
        L.eliminate_zeros()
        return L
    if H is not None:
        L = liouville(H)
    else:
        L = lindbladian(cops[0])
        cops = cops[1:]
    for c in cops:
        L += lindbladian(c)
    L.eliminate_zeros()
    return L

def applySuperoperator(L, rho):
    shape = rho.shape
    N = rho.shape[0]
    rho2 = rho.transpose().reshape(N**2,1)
    
    rho2 = L.dot(rho2)
    out =  rho2.reshape(*shape)
    return tidyup(out, dims = rho.dims)

###############################################################
#################    Fokker Planck Space      #################
###############################################################
def differentiation_matrix(N,m,method='optimal',boundary='periodic'):
    return Qcsc_matrix(_differentiation_matrix(N,m,method=method,boundary=boundary))

def concatenate(*states, dims=None):
    if len(states) == 1:
        states = states[0]
    return Qcsc_matrix(_sp.vstack(states),dims=dims)

def split(state, n:int, new_dims=None):
    d = state.shape[0]
    if d % n != 0:
        raise ValueError("State dimension must be divisible by n.")
    L = d//n
    return [Qcsc_matrix(state[i*L:(i+1)*L],dims=new_dims) for i in range(n)]


