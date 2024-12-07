import numpy as _np
import functools
import scipy.linalg as _la
from scipy.special import gammaln as _gammaln
from .. import core
from typing import Union
import warnings 
from .. trivial import flatten

class QNdarray(_np.ndarray): 
    """ Numpy analogue of a qutip object."""

    def __new__(cls, input_array, dims=None):
        if isinstance(input_array, QNdarray):
            if dims is None:
                return input_array
            else:
                input_array.set_dims(dims)
                return input_array
        else:
            # Input array is an already formed ndarray instance.
            # We first cast to be our class type.
            # Ensure that a QNdarray is always 2-dimensional.
            indim = len(_np.shape(input_array))
            if indim > 2: 
                raise ValueError("QNdarray has to be 2-dimensional.")
            if indim == 0: 
                obj =  _np.asarray([[input_array]]).view(cls)
            elif indim == 1: 
                obj =  _np.asarray(_np.transpose([input_array])).view(cls)
            else : 
                obj = _np.asarray(input_array).view(cls)
            # Note: At this point, array_finalize has been called already.
            # Add the new attribute to the created instance.
            # Set the new attributes to the value passed.     
            if dims is None:
                obj.dims = [[obj.shape[0]], [obj.shape[1]]]
            else:
                obj.set_dims(dims)
            # Finally, we must return the newly created object:
            return obj

    def __array_finalize__(self, obj):
        # See InfoArray.__array_finalize__ for comments.
        if obj is None: 
            return
        indim = len(_np.shape(self))
        # If the following comment is included, one ensures that the QNdarray can never be > 2D.
        # However, some native numpy functions (like np.kron) can then no longer applied on a QNdarray.
        # Note that this would also break compatibility with SimOS code which heavily relies on native numpy functions.
        # If indim > 2 : raise ValueError("QNdarray has to be 2-dimensional.")
        if indim == 2:
            if hasattr(self, 'dims'): 
                self.dims = self.dims
            else: 
                self.dims = getattr(obj, 'dims', [[self.shape[0]], [self.shape[1]]])
        
    # Change Class Methods:
    def __mul__(self,other):
        """
        MULTIPLICATION with QNdarray on LEFT [ ex. Qobj*4 ]
        """
        if isinstance(other, QNdarray):
            out = self@other
            out.set_dims([self.dims[0], other.dims[1]])
            return out
        elif isinstance(other, list):
            return [self*i for i in other]
        else:
            return super().__mul__(other)

    def __rmul__(self,other):
        """
        MULTIPLICATION with QNdarray on RIGHT [ ex. 4*Qobj ]
        """
        if isinstance(other, QNdarray):
            out = other@self
            out.set_dims([other.dims[0], self.dims[1]])
            return out
        elif isinstance(other, list):
            return [self*i for i in other]
        else:
            return super().__rmul__(other)
    
    def __pow__(self,other):
        """
        EXPONENTIATION of a QNdarray [ ex. Qobj**4 ]
        """
        if isinstance(other, (int, float, complex, _np.integer, _np.floating, _np.complexfloating)):
            return _np.linalg.matrix_power(self, other)
        else:
            raise TypeError("Incompatible object for power")

    def __matmul__(self, other):
        if isinstance(other, QNdarray):
            return tidyup(_np.matmul(self, other), dims = [self.dims[0], other.dims[1]])
        else:
            raise TypeError("Incompatible object for multiplication")

    def __rmatmul__(self, other):
        if isinstance(other, QNdarray):
            return tidyup(_np.matmul(other,self), dims = [other.dims[0], self.dims[1]])
        else:
            raise TypeError("Incompatible object for multiplication")

    def set_dims(self, dims : list):
        if len(dims) == 2 and _np.prod(dims[0]) == self.shape[0] and _np.prod(dims[1]) == self.shape[1]:
            self.dims = dims
        else:
            raise ValueError("Dims cannot be set to the specified value.")
        
    # Add Class Methods:

    def dag(self):
        out = _np.transpose(_np.conjugate(self))
        out.dims = [self.dims[1], self.dims[0]]
        return out 

    def conj(self):
        out = _np.conjugate(self)
        return out 

    def trans(self):
        out = _np.transpose(self)
        out.dims = [self.dims[1], self.dims[0]]
        return out

    def tr(self, atol = 1e-10):
        out = _np.trace(self)
        # Check if the trace is real valued; if yes do type conversion.
        if _np.abs(_np.imag(out)) <= atol: 
            return float(_np.real(out))
        else: 
            return out

    def diag(self):
        out = _np.diag(self,k=0)
        return _np.array(out)

    def unit(self):
        out = self/_np.linalg.norm(self, 'nuc')
        return out
    
    def transform(self, U):
        out = _np.matmul(U, _np.matmul(self ,_np.transpose(U)))
        out[_np.abs(out) < 1e-10] = 0
        return out
    
    def expm(self):
        out = _la.expm(self)
        return QNdarray(out, dims = self.dims)
    
    def ptrace(self, sel):
        # Currently ptrace is only working for operators. Catch that.
        if not isoper(self):
            raise ValueError("Numpy backend only supports partial trace for operators.")
        # If a single element is provided, still make a list.
        if isinstance(sel, int):
            sel = [sel]
        # Prepare dimensions. 
        dims = self.dims
        N = len(dims[0]) # how many spins 
        dimflat = flatten(dims)
        i = _np.arange(0, N, 1) # all axis 1 
        j = i + N # all axis 2 
        # Determine over which spins one has to trace.
        # All spins whose index is not in the selection list (sel is those which are being kept).
        sel.sort() 
        all = list(_np.arange(0,N, 1))
        notsel = [x for x in all if x not in sel] # x ist automatisch sortiert 
        # Reshape the array. 
        objt = _np.reshape(_np.array(self), dimflat)
        # Trace out. 
        for ind_A, A in enumerate(notsel): # loop over the components which we want to trace out
            axis1 = i[A]-ind_A
            axis2 = j[A]-2*ind_A
            objt = _np.trace(objt, axis1 = axis1, axis2 = axis2)
        # Final reshape.
        dimnew = [ [dims[0][s] for s in sel], [dims[1][s] for s in sel]] # dimension of obj after ptrace 
        objt = _np.reshape(objt, [_np.prod(dimnew[0]), _np.prod(dimnew[1])]) # reshape properly, only necessary if more than one componend has not been traced over 
        return QNdarray(objt, dims = dimnew)
    
    def eigenstates(self,is_hermitian=True):
        if is_hermitian:
            vals, vecs = _la.eigh(self)
        else:
            vals, vecs = _la.eig(self)
        # Convert to QNdarray.
        new_dims = [self.dims[0], [1 for i in self.dims[1]]]
        vecs = [tidyup(vecs[:,i], dims = new_dims) for i in range(len(vals))]
        return vals, vecs

    def permute(self, order : list):
        # Prepare dimensions. 
        dims = self.dims
        N = len(dims[0]) # how many spins 
        dimflat = flatten(dims)
        i = _np.arange(0, N, 1) # all axis 1 
        j = i + N # all axis 2 
        # Reshape the array. 
        objt = _np.reshape(_np.array(self), dimflat)
        # Perform the permutation.
        oldorder = list(i)+list(j)
        neworder = list(order)+list([o+N for o in order])
        objt = _np.moveaxis(objt, oldorder, neworder)
        # Reshape back.
        objt = _np.reshape(objt, [_np.prod(dims[0]), _np.prod(dims[1])])
        # Settle new dimensions.
        newdims = [dims[0][i] for i in order]
        return QNdarray(objt, dims = [newdims, newdims])
    
###############################################################
# Utility functions to convert data types and change elements:#
###############################################################

def tidyup(arg , dims=None, atol = None):
    # Type conversion with dims setting for non-QNdarray.
    if not isinstance(arg, QNdarray):
        arg = QNdarray(arg, dims = dims)
    else:
    # Only set dims if no type conversion is necessary.
        if dims is not None:
            arg.set_dims(dims)
    return arg

def data(arg : Union[_np.ndarray, QNdarray]):
    # Return the data of the quantum object as a numpy array.
    if not isinstance(arg, (_np.ndarray, QNdarray)):
        raise ValueError("Invalid input data type.")
    else:
        return _np.array(arg)

def changeelement(matrix, ind1 : int, ind2 : int, newvalue : Union[int, float, complex, list, tuple, _np.ndarray, QNdarray]): 
    matrix = matrix.copy()
    # Set a scalar value. 
    if len(_np.shape(newvalue)) == 0:
        matrix[ind1, ind2] = newvalue
    # Set an array or a matrix. 
    else:
        newvalue = _np.array(newvalue)
        if len(newvalue.shape) > 2:
            raise ValueError("Changelement only accepts scalars, arrays or 2D matrices.")
        if len(newvalue.shape) == 1:
            newvalue = _np.array([newvalue])
        matrix[ind1:ind1+newvalue.shape[0], ind2:ind2+newvalue.shape[1]] = newvalue
    return matrix


###############################################################
#################     Math / Generic        ###################
###############################################################

def expect(oper, state):
    if isbra(state):
        state = state.dag()
    vector = isket(state)
    if vector:
        return _np.dot(state.conj().T, _np.dot(oper, state))
    else:
        return _np.trace(_np.matmul(oper,state))

def tensor(*args):
    if len(args) > 1:
        operators = [*args]
    else:
        operators = args[0]
    out = functools.reduce(_np.kron, operators)
    out.dims = [flatten([i.dims[0] for i in operators]), flatten([i.dims[1] for i in operators])]
    return out 

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
            M[pre:pre+_np.shape(operator)[0], pre2:pre2+_np.shape(operator)[1]] = operator
            pre += _np.shape(operator)[0]
            pre2 += _np.shape(operator)[1]
    M = QNdarray(M, dims = None)   
    return M 

def ddrop(operator, idx):
    for axis in (0,1):
        operator = _np.delete(operator, idx, axis = axis)
    operator = QNdarray(operator, dims = [[_np.shape(operator)[0]], [_np.shape(operator)[1]]])
    return operator

def dpermute(operator, neworder):
    if len(operator.shape) != 2 or operator.shape[0] != operator.shape[1]:
        raise ValueError("Dpermute can only be applied to 2-dimensional, square quantum objects.")
    operator =  operator[:, neworder]
    operator = operator[neworder,:]
    operator.dims = [[_np.shape(operator)[0]], [_np.shape(operator)[1]]]
    return operator

def block_diagonal(L_list):
    Lnew = _la.block_diag(*L_list)
    dims = [[len(L_list),L_list[0].shape[0]],[len(L_list),L_list[0].shape[1]]]
    obj =  tidyup(Lnew,dims=dims)
    return obj

def cg(j1, m1, j2, m2, J, M):
    # Match the notation to Varshalovich, Section 8.2.1
    cg = 0
    c = J
    gam = M
    a = j1
    alp = m1
    b = j2
    bet = m2

    if not ((J % 1 == 0 or J % 1 == 0.5) and (M % 1 == 0 or _np.abs(M % 1) == 0.5) and 
            (j1 % 1 == 0 or j1 % 1 == 0.5) and (m2 % 1 == 0 or _np.abs(m2 % 1) == 0.5) and 
            (j2 % 1 == 0 or j2 % 1 == 0.5) and (m1 % 1 == 0 or _np.abs(m1 % 1) == 0.5)):
        raise ValueError('All arguments must be integer or half-integer.')
    
    if not (isinstance(J, (int, float, _np.int32, _np.float64)) and isinstance(M, (int, float, _np.int32, _np.float64)) and 
            isinstance(j1, (int, float, _np.int32, _np.float64)) and isinstance(m1, (int, float, _np.int32, _np.float64)) and 
            isinstance(j2, (int, float, _np.int32, _np.float64)) and isinstance(m2, (int, float, _np.int32, _np.float64))):
        raise ValueError('All arguments must be numeric.')

    if not (j1 >= 0 and j2 >= 0 and J >= 0):
        return 0

    # Run zero tests (Stage I)
    prefactor_is_nonzero = (a + alp >= 0) and (a - alp >= 0) and \
                           (b + bet >= 0) and (b - bet >= 0) and \
                           (c + gam >= 0) and (c - gam >= 0) and \
                           (gam == alp + bet)

    # Proceed if appropriate
    if prefactor_is_nonzero:

        # Run zero tests (Stage II)
        delta_is_nonzero = (a + b - c >= 0) and (a - b + c >= 0) and \
                           (-a + b + c >= 0) and (a + b + c + 1 >= 0)

        # Proceed if appropriate
        if delta_is_nonzero:

            # Run zero tests (Stage III)
            lower_sum_limit = max([alp - a, b + gam - a, 0])
            upper_sum_limit = min([c + b + alp, c + b - a, c + gam])
            sum_is_nonzero = (upper_sum_limit >= lower_sum_limit)
            # Proceed if appropriate
            if sum_is_nonzero:
                # Compute Equation 8.2.1(5) using log-factorials
                for z in _np.arange(lower_sum_limit, upper_sum_limit + 1):
                    cg += ((-1) ** (b + bet + z)) * _np.sqrt(2 * c + 1) * \
                          _np.exp(_gammaln(1+c + b + alp - z) + _gammaln(1+a - alp + z) - _gammaln(1+z) - \
                                 _gammaln(1+c - a + b - z) - _gammaln(1+c + gam - z) - _gammaln(1+a - b - gam + z) + \
                                 (_gammaln(1+a + b - c) + _gammaln(1+a - b + c) + _gammaln(1-a + b + c) + \
                                  _gammaln(1+c + gam) + _gammaln(1+c - gam) - _gammaln(1+a + b + c + 1) - \
                                  _gammaln(1+a + alp) - _gammaln(1+a - alp) - _gammaln(1+b + bet) - \
                                  _gammaln(1+b - bet)) / 2)

    return cg


def WignerD(J, alpha, beta, gamma):
    ops = core.spinops({"val": J, "type": "default"}, method = "numpy", prefix = "")
    return  _np.array((-1j*ops["z"]*alpha).expm()*(-1j*ops["y"]*beta).expm()*(-1j*ops["x"]*gamma).expm())


###############################################################
#################    Hilbert Space          ###################
###############################################################

# Kets, bras, operators and vectoroperators and their conversion. 

def isket(obj):
    return (obj.shape[1] == 1)

def isbra(obj):
    return (obj.shape[0] == 1)

def isoper(obj):
    return (isket(obj) is False and isbra(obj) is False and len(_np.shape(obj.dims[0])) == 1 and len(_np.shape(obj.dims[1])) == 1)

def ket2dm(ket): 
    if isket(ket):
        return QNdarray(_np.dot(ket, ket.dag()), dims = [ket.dims[0], ket.dims[0]])
    if isbra(ket):
        return QNdarray(_np.dot(ket.dag(), ket), dims = [ket.dims[1], ket.dims[1]])
    else:
        raise ValueError("Invalid input shape.")
    
def dm2ket(dm, tol=1e-10):
    if dm.tr() < 1-tol or dm.tr() > 1+tol:
        raise ValueError('Density matrix is not normalized.')
    if (dm**2).tr() < (1-tol):
        raise ValueError('Density matrix is not a pure state.')
    if _np.max(dm.diag()) >= 1-tol:
        return QNdarray(_np.transpose([dm.diag()]), dims = [dm.dims[0], [1 for i in dm.dims[1]]])
    else:
        # Find the decmposition into the pure state that produces the density matrix
        evals, evecs = dm.eigenstates()
        idx = _np.argmax(evals)
        return QNdarray(evecs[idx],dims = [dm.dims[0], [1 for i in dm.dims[1]]])

# Part 2: Construction of standard kets, bras and operators:

def ket(m, n):
    vec = bra(m,n)
    return vec.trans()

def bra(m, n):
    vec = _np.atleast_2d(_np.zeros(m))
    vec[0,n] = 1
    return QNdarray(vec)

def jmat(j,op_spec):
    if (_np.fix(2 * j) != 2 * j) or (j < 0):
        raise TypeError('j must be a non-negative integer or half-integer')
    N = int(2*j+1)

    m = _np.arange(j, -j - 1, -1, dtype=complex)
    pm_data = (_np.sqrt(j * (j + 1.0) - (m + 1.0) * m))[1:]

    if op_spec == '+':
        return QNdarray(_np.diag(pm_data, k=1))
    elif op_spec == '-':
        return QNdarray(_np.diag(pm_data, k=-1))
    elif op_spec == 'x':
        return QNdarray(0.5 * (_np.diag(pm_data, k=1) + _np.diag(pm_data, k=-1)))
    elif op_spec == 'y':
        return QNdarray(-0.5 * 1j * (_np.diag(pm_data, k=1)  - _np.diag(pm_data, k=-1)))
    elif op_spec == 'z':
        data = _np.array([j-k for k in range(N)], dtype=complex)
        return QNdarray(_np.diag(data))
    else:
        raise TypeError('Invalid operator specification')

def identity(N, dims = None):
    out = _np.eye(N,dtype=complex)
    return QNdarray(out, dims = dims)

def diags(v : Union[QNdarray, tuple, list, _np.ndarray], k : int = 0, dims = None):
    # utilize native numpy function 
    out = _np.diag(v,k)
    # if a diagonal was extracted, ignore dims 
    if len(_np.shape(v)) == 2:
        if dims is not None:
            warnings.warn("Dims keyword argument will be ignored.")
        return _np.array(out) 
    # if quantum object was build, set dims 
    else:
        return QNdarray(out, dims = dims)


###############################################################
#################    Liouville Space         ##################
###############################################################

def dm2vec(op):
    # stack the columns of the matrix
    dims = [op.dims, [1]]
    return tidyup(op.ravel('F')[:, None],dims=dims)

def vec2dm(vec):
    s = max(vec.shape)
    s = int(_np.sqrt(s))
    # unstack the columns of the matrix
    return tidyup(vec.reshape([s, s], order='F'),dims= [[s], [s]])
    
def issuper(obj):
    return( len(obj.dims[0]) == 2 and len(obj.dims[1]) == 2)

def spre(H):
    pre = _np.zeros((H.shape[0]**2, H.shape[0]**2),dtype=_np.complex128)
    for i in range(H.shape[0]):
        pre[i*H.shape[0]:(i+1)*H.shape[0], i*H.shape[0]:(i+1)*H.shape[0]] = H
    pre = QNdarray(pre, dims = [[H.dims[0], H.dims[0]], [H.dims[1], H.dims[1]]])
    return pre

def spost(H):
    post = _np.zeros((H.shape[0]**2, H.shape[0]**2),dtype=_np.complex128)
    for i in range(H.shape[0]):
        post[i::H.shape[0], i::H.shape[0]] = H.T
    post = QNdarray(post, dims = [[H.dims[0], H.dims[0]], [H.dims[1], H.dims[1]]])
    return post

def sprepost(H,factor=-1):
    """Computes 1 \\otimes H + factor * H.T \\otimes 1"""
    L = _np.zeros((H.shape[0]**2, H.shape[0]**2),dtype=_np.complex128)
    for i in range(H.shape[0]):
        L[i*H.shape[0]:(i+1)*H.shape[0], i*H.shape[0]:(i+1)*H.shape[0]] += H
        L[i::H.shape[0], i::H.shape[0]] += factor*H.T
    L = QNdarray(L, dims = [[H.dims[0], H.dims[0]], [H.dims[1], H.dims[1]]])
    return L

def liouville(H):
    return -1j*sprepost(H,factor=-1)

def lindbladian(a,b=None):
    if b is None:
        b = a
    #bdag = b.dag()
    #adag = a.dag() 
    ad_b = a.dag() * b
    #L = QNdarray(_np.kron(bdag,a), dims = [[a.dims[0], a.dims[0]], [a.dims[1], a.dims[1]]]) 
    #L = L - 0.5*sprepost(_np.dot(adag,b),factor=1)
    L = spre(a) * spost(b.dag()) - 0.5 * spre(ad_b) - 0.5 * spost(ad_b)
    return L

def liouvillian(H,cops=None):
    if isinstance(cops,list) or isinstance(cops,tuple):
        if len(cops)==0:
            cops = None
    if cops is None:
        return liouville(H)
    if H is not None:
        L = liouville(H)
    else:
        L = lindbladian(cops[0])
        cops = cops[1:]
    for c in cops:
        L += lindbladian(c)
    return L

def applySuperoperator(L, rho):
    rho2 = rho.ravel('F')[:, None]
    shape = rho.shape
    rho2 = L.dot(rho2)
    return rho2.reshape(shape, order='F')

###############################################################
#################    Fokker Planck Space      #################
###############################################################

def _cot(x):
    return _np.tan(_np.pi/2-x)

def _csc(x):
    return 1/_np.sin(x)

def differentiation_matrix(N,m,method='optimal',boundary='periodic'):
    if boundary == 'periodic':
        if method not in ['optimal','fft']:
            raise ValueError('Invalid method')
        h = 2*_np.pi/N
        kk = _np.arange(0,N)+1
        n1 = int(_np.floor((N-1)/2))
        n2 = int(_np.ceil((N-1)/2))
        if m==1 and method == 'optimal':
            if N%2 == 0:
                topc = _cot(_np.arange(1,n2+1)*h/2)
                topc = _np.concatenate((topc,-_np.flipud(topc[0:n1])))
                col1 = _np.concatenate((_np.zeros(1),0.5*(-1)**kk[:-1]*topc))
            else:
                topc = _csc(_np.arange(1,n2+1)*h/2)
                topc = _np.concatenate((topc,_np.flipud(topc[0:n1])))
                col1 = _np.concatenate((_np.zeros(1),0.5*(-1)**kk[:-1]*topc))
            row1 = -col1      
        elif m == 2 and method == 'optimal':
            if N%2 == 0:
                topc = _csc(_np.arange(1,n2+1)*h/2)**2
                topc = _np.concatenate((topc,_np.flipud(topc[0:n1])))
                col1 = _np.concatenate(([-_np.pi**2/(3*h**2)-1/6],-1/2*(-1)**kk[:-1]*topc))
            else:
                topc = _csc(_np.arange(1,n2+1)*h/2)*_cot(_np.arange(1,n2+1)*h/2)
                topc = _np.concatenate((topc,-_np.flipud(topc[0:n1])))
                col1 = _np.concatenate(([-_np.pi**2/(3*h**2)+1/12],-1/2*(-1)**kk[:-1]*topc))
            row1 = col1
        elif m == 0:
            col1=_np.concatenate(([1],_np.zeros(N-1)))
            row1=col1
        else:
            N1 = int(_np.floor((N-1)/2))
            N2 = (-N // 2) * ((m + 1) % 2) * _np.ones((N + 1) % 2)
            mwave = 1j * _np.concatenate((_np.arange(0, N1+1), N2, _np.arange(-N1, 0)))
            col1 = _np.real(_np.fft.ifft((mwave**m) * _np.fft.fft(_np.concatenate(([1], _np.zeros(N-1))))))
            if m%2 == 0:
                row1 = col1
            else:
                col1 = _np.concatenate(([0], col1[1:N]))
                row1 = -col1

        return _la.toeplitz(col1,row1)
    else:
        raise NotImplementedError('Only periodic boundary conditions are implemented')
    
def concatenate(*states, dims=None):
    if len(states) == 1:
        states = states[0]
    return tidyup(_np.concatenate([state for state in states], axis=0),dims=dims)

def split(state, n:int, new_dims=None,is_hilbert=False):
    d = state.shape[0]
    if d % n != 0:
        raise ValueError("State dimension must be divisible by n.")
    L = d//n
    if new_dims is None:
        if is_hilbert:
            new_dims = [[state.dims[0][-1]], [1]]
    #return [tidyup(state[i*L:(i+1)*L],dims=new_dims) for i in range(n)]
    segment = [state[i*L:(i+1)*L] for i in range(n)]
    if not is_hilbert:
        # make segments square matrices
        L2 = int(_np.sqrt(L))
        segment = [_np.reshape(segment[i],(L2,L2)).T for i in range(n)]
    return [tidyup(segment[i],dims=new_dims) for i in range(n)]


#### Numeric math methods ####


# Set up the transformation matrix for diagonalization using the eigenvectors. 
## Eigenvalues are sorted such that a proper rotation matrix (det(mat) = 1) is obtained.
## The obtained matrix T_toeigen can therefore be used as an input to create a rotation object of scipy. 
def matrix_eigenstates(matrix):
    if _la.ishermitian(matrix):
        vals_tmp, vecs_tmp = _la.eigh(matrix)
    else:
        vals_tmp, vecs_tmp = _la.eig(matrix)
    idx_tmp = _np.flip(_np.argsort(_np.abs(vals_tmp)))
    idx = [idx_tmp[1], idx_tmp[2], idx_tmp[0]]
    vals = vals_tmp[idx]
    T_toeigen = _np.zeros(_np.shape(matrix), dtype = complex) 
    for idx, i in enumerate(idx):
        T_toeigen[idx, : ] = vecs_tmp[:, i]   
    if (_np.abs(_np.linalg.det(T_toeigen)-1)) < 0.1:
        return vals,  T_toeigen
    else:
        idx =  [idx_tmp[2], idx_tmp[1], idx_tmp[0]]
        vals = vals_tmp[idx]
        T_toeigen = _np.zeros(_np.shape(matrix), dtype = complex) 
        for idx, i in enumerate(idx):
            T_toeigen[idx, : ] = vecs_tmp[:, i]   
        return vals, T_toeigen

def rotate(matrix, R, inverse = False):
    Rmat = R.as_matrix()
    if inverse:
        Rmat = Rmat.transpose()
    out = _np.matmul(Rmat, _np.matmul(matrix,Rmat.transpose()))
    out[_np.abs(out) < 1e-6] = 0
    return out  

def pow_fun(a,b):
    if b == -1:
        return 1/a
    else:
        return _np.power(a,b)