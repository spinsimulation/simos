import sympy as _sp
import numpy as _np
import functools
import sympy.physics.quantum as _spq
import copy
from typing import Union
import warnings
from .. trivial import flatten

class QSMatrix(_sp.Matrix):
    """ Sympy analogue of a qutip object. """

    def __new__(cls, *args, **kwargs):
        obj = _sp.Matrix.__new__(cls, *args, **kwargs)
        if "dims" in kwargs:
            if kwargs["dims"] is None:
                obj.dims = [[obj.shape[0]], [obj.shape[1]]]
            else:
                obj.set_dims(kwargs["dims"])
        else:
            obj.dims = [[obj.shape[0]], [obj.shape[1]]]
        return obj
    
    def __add__(self, other):
        if hasattr(self,'dims'):
            dims = self.dims
        else:
            dims = list(self.shape)
        if isinstance(other, (int,float,_sp.Symbol)):
            out = self.applyfunc(lambda x: x + other)
        else:
            out = super().__add__(other)
        out.dims = dims
        return out
    
    def __radd__(self, other):
        if hasattr(self,'dims'):
            dims = self.dims
        else:
            dims = list(self.shape)
        if isinstance(other, (int,float,_sp.Symbol)):
            out = self.applyfunc(lambda x: x + other)
        else:
            out = super().__radd__(other)
        out.dims = dims
        return
    
    def __sub__(self, other):
        if hasattr(self,'dims'):
            dims = self.dims
        else:
            dims = [[self.shape[0]], [self.shape[1]]]
        if isinstance(other, (int,float,_sp.Symbol)):
            out = self.applyfunc(lambda x: x - other)
        else:
            out = super().__sub__(other)
        out.dims = dims
        return out
    
    def __rsub__(self, other):
        if hasattr(self,'dims'):
            dims = self.dims
        else:
            dims = [[self.shape[0]], [self.shape[1]]]
        if isinstance(other, (int,float,_sp.Symbol)):
            out = self.applyfunc(lambda x: other - x)
        else:
            out = super().__rsub__(other)
        out.dims = dims
        return out
    
    def __mul__(self, other):
        if hasattr(self,'dims'):
            if hasattr(other,'dims'):
                dims = [self.dims[0], other.dims[1]]
            else:
                dims = self.dims
        else:
            dims = [[self.shape[0]], [self.shape[1]]]     
        out = super().__mul__(other)
        out.dims = dims
        return out
    
    def __rmul__(self, other):
        if hasattr(self,'dims'):
            if hasattr(other,'dims'):
                dims = [self.dims[1], other.dims[0]]
            else:
                dims = self.dims
        else:
                dims = [[self.shape[0]], [self.shape[1]]]        
        out = super().__rmul__(other)
        out.dims = dims
        return out

    def __pow__(self,other):
        if hasattr(self,'dims'):
            dims = self.dims
        else:
            dims = [[self.shape[0]], [self.shape[1]]]
        out = super().__pow__(other)
        out.dims = dims
        return out 

    def set_dims(self, dims : list):
        if len(dims) == 2 and _np.prod(dims[0]) == self.shape[0] and _np.prod(dims[1]) == self.shape[1]:
            self.dims = dims
        else:
            raise ValueError("Dims cannot be set to the specified value.")

    # Add Class Methods 
    def dag(self):
        out = self.adjoint()
        out.dims = [self.dims[1], self.dims[0]]
        return out

    def conj(self):
        out = self.conjugate()
        out.dims = self.dims
        return out 

    def trans(self):
        out = self.transpose()
        out.dims = [self.dims[1], self.dims[0]]
        return out

    def tr(self):
        out = self.trace()
        return out
    
    def diag(self):
        return diags(self, k = 0)
    
    def unit(self):
        #nuclear norm
        #norm = ((self.dag()*self).applyfunc(_sp.sqrt)).tr()
        norm = sum(self.singular_values())
        out = self/norm
        return out

    def transform(self, U):
        U = tidyup(U)
        out = U*(self*U.transpose())
        out.dims = self.dims
        #out = _np.matmul(U, _np.matmul(self ,_np.transpose(U)))
        #out = out[_np.abs(out) >= 1e-10]
        return out
    
    def expm(self,simplify = False):
        out = _sp.exp(self)
        if simplify:
            out = _sp.simplify(out)
        out = QSMatrix(out, dims = self.dims)
        return out
    
    def ptrace(self, sel):
        # Currently ptrace is only working for operators. Catch that.
        if not isoper(self):
            raise ValueError("Sympy backend only supports partial trace for operators.")
        # handle type of sel
        if isinstance(sel, int):
            sel = [sel]
        dims = self.dims
        N = len(dims[0])# wie viele spins sind da
        dimflat = flatten(dims)
        i = _np.arange(0, N, 1) # all axis 1 
        j = i + N # all axis 2 
        # determine over which spins one has to trace -> all spin whose index is not in the selection list (sel is those which are being kept)
        sel.sort() 
        all = list(_np.arange(0,N, 1))
        notsel = [x for x in all if x not in sel] # x ist automatisch sortiert 
        # reshape the array 
        objt = _sp.MutableDenseNDimArray(self).reshape(*dimflat)
        # trace out 
        for ind_A, A in enumerate(notsel): # loop over the components which we want to trace out
            axis1 = i[A]-ind_A
            axis2 = j[A]-2*ind_A
            objt = _sp.tensorcontraction(objt, (axis1, axis2)) 
        # final reshape 
        dimnew = [ [dims[0][s] for s in sel], [dims[1][s] for s in sel]] # dimension of obj after ptrace 
        objt = objt.reshape( int(_np.prod(dimnew[0])), int(_np.prod(dimnew[1]))) # reshape properly, only necessary if more than one componend has not been traced over 
        objt = QSMatrix(objt, dims = dimnew )
        #objt.dims = dimnew
        return objt

    def eigenstates(self):
        eigen = self.eigenvects()
        # eigen contains a tuple of eigenvalues, algebraic multiplicity and eigenvectors
        # unpack eigenvectors and eigenvalues to separate lists, double the eigenvalues for degeneracies
        eigvals = []
        eigvecs = []
        new_dims = [self.dims[0], [1 for i in self.dims[1]]]
        for eig in eigen:
            mult = eig[1]
            eigvals += [eig[0]]*mult
            # convert eigenvectors to QSMatrix with tidyup
            for i in range(mult):
                eigvecs.append(tidyup(eig[2][i], dims = new_dims))
        return eigvals, eigvecs

    def permute(self, order : list):
        # Prepare dimensions. 
        dims = self.dims
        N = len(dims[0]) # how many spins 
        dimflat = flatten(dims)
        i = _np.arange(0, N, 1) # all axis 1 
        j = i + N # all axis 2 
        # Reshape the array. 
        objt = data(self)
        objt = _np.reshape(objt, dimflat)
        # Perform the permutation.
        oldorder = list(i)+list(j)
        neworder = list(order)+list([o+N for o in order])
        objt = _np.moveaxis(objt, oldorder, neworder)
        # Reshape back.
        objt = _np.reshape(objt, [_np.prod(dims[0]), _np.prod(dims[1])])
        # Settle new dimensions.
        newdims = [dims[0][i] for i in order]
        return tidyup(objt, dims = [newdims, newdims])
    

    # The following methods are sympy-specific.
    # They do not exist for quantum objects of other backends.

    def simplify(self,*args, **kwargs):
        out = _sp.simplify(self, *args, **kwargs)
        return out
    
    def exp2trig(self):
        return _sp.expand(_sp.simplify(self),complex=True)
        
    def trig2exp(self):
        out = self.applyfunc(lambda x: x.rewrite(_sp.exp))
        return out
    
    def applyfunc_elementwise(self, func):
        out = self.applyfunc(func)
        return out
    

###############################################################
# Utility functions to convert data types and change elements:#
###############################################################

def tidyup(arg, dims=None):
    # Type conversion with dims setting for non-QSMatrix.
    if not isinstance(arg, QSMatrix):
        arg = QSMatrix(arg, dims = dims)
    # Otherwise set dims, but do not overwrite existing with None.
    else:
        if dims is not None:
           arg.set_dims(dims)
    return arg

def data(arg : QSMatrix):
    # Return the data of the quantum object as a numpy array.
    return _np.array(arg)

def changeelement(matrix, ind1, ind2, newvalue : Union[int, float, complex, list, tuple, _np.ndarray]): 
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
    return matrix


###############################################################
#################     Math / Generic        ###################
###############################################################

def expect(oper, state):
    ket = isket(state)
    if ket:
        return state.H*oper*state
    else:
        return _sp.trace(oper * state)

def tensor(*args):
    if len(args) > 1:
        operators = [*args]
    else:
        operators = args[0]
    out = functools.reduce(_spq.TensorProduct, operators)
    out.dims = [flatten([i.dims[0] for i in operators]), flatten([i.dims[1] for i in operators])]
    return out 


def directsum(*args):
    if len(args) > 1:
        operators  = [*args]
    else:
        operators = args[0]
    dimtot = sum([i if isinstance(i, int) else _np.shape(i)[0]  for i in operators])
    dimtot2 = sum([i if isinstance(i, int) else _np.shape(i)[1]  for i in operators])    
    M = _sp.zeros(dimtot, dimtot2)
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
    M = QSMatrix(M, dims = None)   
    return M 

def ddrop(operator, idx):
    operator = operator.copy()
    for i_ind, i in enumerate(idx):
        operator.row_del(i-i_ind)
    for i_ind, i in enumerate(idx):
        operator.col_del(i-i_ind)
    operator = QSMatrix(operator, dims = None)
    return operator

def dpermute(operator, neworder):
    if len(_np.shape(operator)) != 2 or _np.shape(operator)[0] != _np.shape(operator)[1]:
        raise ValueError("Dpermute can only be applied to 2-dimensional, square quantum objects.")
    operator = data(operator)
    operator =  operator[:, neworder]
    operator = operator[neworder,:]
    return QSMatrix(operator, dims = [[_np.shape(operator)[0]], [_np.shape(operator)[1]]])

def block_diagonal(L_list):
    Lnew = _sp.BlockDiagMatrix(*L_list)
    #dims = [[L_list[0].shape[0],len(L_list)], [L_list[0].shape[1],len(L_list)]]
    dims = [[len(L_list),L_list[0].shape[0]],[len(L_list),L_list[0].shape[1]]]
    obj =  tidyup(Lnew,dims=dims)
    return obj

def sympy_outer(A, B):
    return _sp.Matrix([[a*b for b in B] for a in A])

def cg(j1, m1, j2, m2, J, M):
    j1 = _sp.Rational(j1)
    m1 = _sp.Rational(m1)
    j2 = _sp.Rational(j2)
    m2 = _sp.Rational(m2)
    J = _sp.Rational(J)
    M = _sp.Rational(M)
    cg = _sp.physics.quantum.cg.CG(j1, m1, j2, m2, J, M).doit()
    return cg

###############################################################
#################    Hilbert Space          ###################
###############################################################  

# Kets, bras, operators and their conversion. 

def isket(obj):
    return (obj.shape[1] == 1)

def isbra(obj):
    return (obj.shape[0] == 1)

def isoper(obj):
    return (isket(obj) is False and isbra(obj) is False and len(_np.shape(obj.dims[0])) == 1 and len(_np.shape(obj.dims[1])) == 1)

def ket2dm(ket):
    if isket(ket):
        return QSMatrix(ket * ket.H, dims=[ket.dims[0], ket.dims[0]])
    if isbra(ket):
        return QSMatrix(ket.H * ket, dims=[ket.dims[1], ket.dims[1]])
    else:
        raise ValueError("Invalid input shape.")

def dm2ket(dm, tol = 1e-10):
    if isinstance(dm.tr().simplify(), _sp.Symbol) or isinstance((dm**2).tr().simplify(), _sp.Symbol):
        raise ValueError("Dm2ket cannot be performed because a value is symbolic")
    if dm.tr() < 1-tol or dm.tr() > 1+tol:
        raise ValueError('Density matrix is not normalized.')
    if (dm**2).tr() < (1-tol):
        raise ValueError('Density matrix is not a pure state.')
    if _np.max(dm.diag()) >= 1-tol:
        return tidyup(_np.transpose([dm.diag()]), dims = [dm.dims[0], [1 for i in dm.dims[1]]])
    else:
        # Find the decmposition into the pure state that produces the density matrix
        evals, evecs = dm.eigenstates()
        idx = _np.argmax(evals)
        ket = evecs[idx]
        ket.set_dims([dm.dims[0], [1 for i in dm.dims[1]]])
        return ket.unit()

# Part 2: Construction of standard kets, bras and operators:

def ket(m, n):
    vec = bra(m,n)
    return vec.trans()

def bra(m, n):
    vec = QSMatrix(_sp.zeros(1,m))
    vec[n] = 1
    return tidyup(vec)

def jmat(j,op_spec):
    if (2*j%1 != 0) or (j < 0):
        raise TypeError('j must be a non-negative integer or half-integer')
    N = int(2*j+1)
    # check if j is an integer
    if j%1 == 0:
        j = int(j)
    else:
        j = _sp.Rational(int(2*j),2)
    m = [j-k for k in range(N)]
    pm_data = [(_sp.sqrt(j * (j + 1) - (mi + 1) * mi)) for mi in m][1:]
    if op_spec == '+':
        return tidyup(diags(pm_data, k=1))
    elif op_spec == '-':
        return tidyup(diags(pm_data, k=-1))
    elif op_spec == 'x':
        return tidyup(_sp.Rational(1,2) * (diags(pm_data, k=1) + diags(pm_data, k=-1)))
    elif op_spec == 'y':
        return tidyup(-_sp.Rational(1,2) * _sp.I * (diags(pm_data, k=1) - diags(pm_data, k=-1)))
    elif op_spec == 'z':
        data = [j-k for k in range(N)]
        return tidyup(diags(data))
    else:
        raise TypeError('Invalid operator specification')

def identity(N, dims = None):
    out = _sp.eye(N)
    return QSMatrix(out, dims = dims)

def diags(v : Union[QSMatrix, list, tuple, _np.ndarray], k : int = 0, dims = None):
    # Return the kth diagonal if a 2D quantum object is provided.
    if len(_np.shape(v)) == 2:
        if not isinstance(v, QSMatrix):
            raise ValueError("Diagonal is only extracted from data type QSMatrix. For extraction from  list, tuple, array etc. directly utilize the diags function of the numpy backend.")
        if dims is not None:
            warnings.warn("Dims keyword argument will be ignored.")
        v = _np.array(v)
        out = _np.diag(v, k)
        return out 
    # If a 1D input was provided, build a diagonal matrix.
    else:
        if len(_np.shape(v)) != 1:
            raise ValueError("Diagonal matrix can only be constructed from one dimensional object.")
        else:
            out = _sp.eye(len(v) + abs(k))*0
            for i in range(len(v)):
                if k > 0:
                    out[i,i+k] = v[i]
                else:
                    out[i-k,i] = v[i]
        out = tidyup(out, dims  = dims)
    return out

###############################################################
#################    Liouville Space         ##################
###############################################################

def dm2vec(op):
    # stack the columns of the matrix
    dims = [op.dims, [1]]
    return tidyup(_sp.Matrix(op).transpose().reshape(op.shape[0] * op.shape[1], 1), dims=dims)

def vec2dm(vec):
    s = max(vec.shape)
    s = int(_np.sqrt(s))
    # unstack the columns of the matrix
    return tidyup(_sp.Matrix(vec).reshape(s, s).T, dims=[[s], [s]])

def issuper(obj):
    return( len(obj.dims[0]) == 2 and len(obj.dims[1]) == 2)

def spre(H):
    ident = _sp.eye(H.shape[0])
    pre =  _spq.TensorProduct(ident,H)
    pre =  QSMatrix(pre, dims = [[H.dims[0], H.dims[0]], [H.dims[1], H.dims[1]]])
    return pre

def spost(H):
    ident = _sp.eye(H.shape[0])
    ident = QSMatrix(ident)
    ident.dims = list(ident.shape)
    HT = QSMatrix(H.T)
    HT.dims = list(HT.shape)
    post = _spq.TensorProduct(HT,ident)
    post =  QSMatrix(post, dims = [[H.dims[0], H.dims[0]], [H.dims[1], H.dims[1]]])
    return post

def sprepost(H,factor=-1):
    L = spre(H) + spost(H)*factor
    return L

def liouville(H):
    return sprepost(H,factor=-1)*(-_sp.I)

def lindbladian(a,b=None):
    if b is None:
        b = a
    #bdag = b.conj().T
    #adag = a.conj().T 
    #part1 = QSMatrix(_spq.TensorProduct(bdag,a), dims = [[a.dims[0], a.dims[0]], [a.dims[1], a.dims[1]]])
    #return part1 - sprepost(adag*(b),factor=1)/2
    ad_b = a.dag() * b
    L = spre(a) * spost(b.dag()) - spre(ad_b)/2 - spost(ad_b)/2
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
    shape = rho.shape
    N = rho.shape[0]
    rho2 = rho.T.reshape(N**2,1)
    
    rho2 = L*(rho2)
    return rho2.reshape(*shape)

###############################################################
#################    Fokker Planck Space      #################
###############################################################

def _toeplitz(col1,row1):
    N = len(col1)
    return _sp.Matrix([[col1[(i-j) % N] for j in range(N)] for i in range(N)])

def differentiation_matrix(N,m,method='optimal',boundary='periodic'):
    if boundary == 'periodic':
        if method not in ['optimal']:
            raise ValueError('Invalid method')
        h = 2*_sp.pi/N
        kk = _np.arange(0,N)+1
        n1 = int(_np.floor((N-1)/2))
        n2 = int(_np.ceil((N-1)/2))
        if m==1 and method == 'optimal':
            if N%2 == 0:
                topc = _np.arange(1,n2+1)
                topc = _np.concatenate((topc,-_np.flipud(topc[0:n1])))
                topc = _sp.MutableDenseNDimArray(topc)
                topc = topc.applyfunc(lambda x: _sp.cot(x*h/2))
                pre = (-1)**kk[:-1]
                topc = [0] + [topc[i]*pre[i] for i in range(len(pre))]
                topc = _sp.MutableDenseNDimArray(topc)/2
                col1 = topc
            else:
                topc = _np.arange(1,n2+1)
                topc = _np.concatenate((topc,_np.flipud(topc[0:n1])))
                topc = _sp.MutableDenseNDimArray(topc)
                topc = topc.applyfunc(lambda x: _sp.csc(x*h/2))
                pre = (-1)**kk[:-1]
                topc = [0] + [topc[i]*pre[i] for i in range(len(pre))]
                topc = _sp.MutableDenseNDimArray(topc)/2
                col1 = topc
            row1 = -col1      
        elif m == 2 and method == 'optimal':
            if N%2 == 0:
                topc = _np.arange(1,n2+1)
                topc = _np.concatenate((topc,_np.flipud(topc[0:n1])))
                topc = _sp.MutableDenseNDimArray(topc)
                topc = topc.applyfunc(lambda x: _sp.csc(x*h/2)**2)
                pre = (-1)**kk[:-1]
                topc = [-_sp.pi**2/(3*h**2)-_sp.Rational(1,6)] + [-topc[i]*pre[i]/2 for i in range(len(pre))]
                topc = _sp.MutableDenseNDimArray(topc)
                col1 = topc
            else:
                topc = _np.arange(1,n2+1)
                topc = _sp.MutableDenseNDimArray(topc)
                topc = topc.applyfunc(lambda x: _sp.csc(x*h/2)*_sp.cot(x*h/2))
                topc = [e for e in topc] + [-e for e in reversed(topc[0:n1])]
                topc = _sp.MutableDenseNDimArray(topc)
                pre = (-1)**kk[:-1]
                topc = [-_sp.pi**2/(3*h**2)+_sp.Rational(1,12)] + [-topc[i]*pre[i]/2 for i in range(len(pre))]
                topc = _sp.MutableDenseNDimArray(topc)
                col1 = topc
            row1 = col1
        elif m == 0:
            col1=_np.concatenate([1],_np.zeros(N-1))
            col1 = _sp.MutableDenseNDimArray(col1)
            row1=col1
        else:
            raise NotImplementedError('Only m=0,1,2 are implemented')

        return _toeplitz(col1,row1)
    else:
        raise NotImplementedError('Only periodic boundary conditions are implemented')
    
   
def concatenate(*states, dims=None):
    if len(states) == 1:
        states = states[0]
    concatenated_matrix = _sp.Matrix.vstack(*[_sp.Matrix(state) for state in states])
    return tidyup(concatenated_matrix, dims=dims)

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




#### Symbolic math methods ####
def matrix_eigenstates(input):
    raise NotImplementedError

def rotate(input, R, inverse = False):
    input = _sp.Matrix(input)
    Rmat = R.to_rotation_matrix()
    if inverse:
        Rmat = Rmat.transpose()
    return _sp.MutableDenseNDimArray(Rmat*(input*Rmat.transpose()))

def sqrt_fun(input):
    if hasattr(input, "applyfunc"):
        return input.applyfunc(_sp.sqrt)
    else: 
        return _sp.sqrt(input)
    
def sin_fun(input):
    if hasattr(input, "applyfunc"):
        return input.applyfunc(_sp.sin)
    else:
        return _sp.sin(input)

def cos_fun(input):
    if hasattr(input, "applyfunc"):
        return input.applyfunc(_sp.cos)
    else:
        return _sp.cos(input)

def exp_fun(input):
    if hasattr(input, "applyfunc"):
        return input.applyfunc(_sp.exp)
    else:
        return _sp.exp(input)
    
def pow_fun(input, power):
    if hasattr(input, "applyfunc"):
        return input.applyfunc(lambda x: _sp.Pow(x,power))
    else:
        return _sp.Pow(input, power)

def arctan2_fun(input, input2):
    if len(_np.shape(input)) > 0:
        for ind_i in range(len(input)):
            input[ind_i ]= _sp.atan2(input[ind_i], input2[ind_i])
        return input
    else:
        return _sp.atan2(input, input2)

def multelem_fun(arg1, arg2):
    a = copy.deepcopy(arg1)
    b = copy.deepcopy(arg2)
    # if either a or b is a simple symbol
    if len(_np.shape(a)) == 0 or len(_np.shape(b)) == 0:
        return a * b
    # otherwise
    else: 
        if _np.shape(a) != _np.shape(b):  
            raise ValueError("Only objects with same shape can be multiplied.")
        else:
            result = _sp.matrices.dense.matrix_multiply_elementwise(_sp.Matrix(a), _sp.Matrix(b))
            if _np.shape(result) == _np.shape(a):
                return _sp.MutableDenseNDimArray(result)
            else:
                return _sp.MutableDenseNDimArray(result).transpose()[0]

def matmul_fun(arg1, arg2):
    a = _sp.Matrix(arg1)
    b = _sp.Matrix(arg2)
    return a*b

