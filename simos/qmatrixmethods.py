import numpy as _np
from .constants import *
from .core import *
from . import backends
from typing import Union
import warnings


###########################################################
# GLOBAL WRAPPER FUNCTIONS FOR BACKEND SPECIFIC METHODS
###########################################################

#############  UTILITY FUNCTIONS #############

def tidyup(arg, method = None, dims = None, **kwargs):
    """ Transforms an input argument to a quantum object data type of an available backend. The data types of the individual elements of the quantum object are NOT altered by this routine.  

    :param arg: Input object.
    :param str method: Backend to which type conversion should be performed. If omitted, the target backend is determined via the data type of the input argument.
    :param list dims: Specifies dims attribute of the output quantum object. If dims == None but the input argument is already a quantum object with valid dims, the old dims are retained.
    :param **kwargs:
        See below

    :Keyword Arguments:
        * *atol* (``float``) --
          Values smaller than atol are set to 0. Only used in 'scipy' and 'qutip' backends.     
    
    :return: Quantum object data type.
    """
    # Set target dimensions.
    if hasattr(arg, "dims") and dims is None:
        tdims = arg.dims
    else:
        tdims = dims 
    # Type conversion.
    curr_method = backends.get_backend(arg)
    # If explicitly asked for conversion to another backend. 
    if method is not None and curr_method is not method:
        targ = getattr(getattr(backends,curr_method), 'data')(arg)
        targ = getattr(getattr(backends, method), 'tidyup')(targ, dims = tdims, *kwargs)
    # If desired backend None or same as current. 
    else:
        targ = getattr(getattr(backends,curr_method), 'tidyup')(arg, dims = tdims, *kwargs)
    return targ 


def data(arg):
    """ Returns data of a quantum object as a numpy.ndarray.

    :param arg: Quantum object.
    :return: Numpy.ndarray holding data of input quantum object.
    """
    backend = backends.get_backend(arg)
    return getattr(getattr(backends,backend), 'data')(arg)


def changeelement(matrix, ind1 : int, ind2 : int, newvalue : Union[int, float, complex, list, tuple, _np.ndarray]):
    """
    Sets a single element or a block of elements in an operator or superoperator with a series of values.

    :param matrix: Input operator/superoperator that will be modified.
    :param ind1: Index specifying starting row for data change.
    :param ind2: Index specifying starting column for changed data.
    :param newvalue: Values to be inserted.

    :return: Updated operator/superoperator.

    """
    backend = backends.get_backend(matrix)
    matrix = tidyup(matrix, method = backend)
    return getattr(getattr(backends,backend), 'changeelement')(matrix, ind1, ind2, newvalue)


############# QUANTUM OBJECT ATTRIBUTES  #############

def dag(op):
    """ Dagger(adjoint) of a quantum object.
    
    :param: Quantum object.
    :return: Dagger(adjoint)."""
    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    return op.dag()

def conj(op):
    """ Conjugate of a quantum object.
    
    :param: Quantum object.
    :return: Conjugate."""
    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    return op.conj()

def trans(op):
    """ Transpose of a quantum object.
    
    :param: Quantum object.
    :return: Transpose."""
    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    return op.trans()

def tr(op):
    """Trace of a quantum object.

    :param: Quantum object.
    :return: Trace.
    
    """
    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    return op.tr()


def diag(op):
    """ Diagonal of a quantum object.
    
    :param: Quantum object.
    :return: Diagonal of the quantum object, returned as a numpy.ndarray."""
    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    return op.diag()


def unit(op):
    """ Normalizes a quantum object.
    
    :param: Quantum object.
    :return: Normalized quantum object."""
    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    return op.unit()


def transform(op,U):
    """ Coordinate Transformation of a quantum object.
    
    :param: Quantum object.
    :U: Transformation matrix.
    :return: Transformed quantum object."""
    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    return op.transform(U)

def expm(H):
    """ Matrix exponential of a quantum object.
    
    :param: Quantum object.
    :return: Matrix exponential."""
    backend = backends.get_backend(H)
    H = tidyup(H, method = backend)
    return H.expm()

def ptrace(op, sel , dims=None):
    """Partial trace of a quantum object.
    :param op: Operator.
    :param sel: Separable subsystem that is kept after the partial trace.
    :return: Partial trace.
    
    """
    backend = backends.get_backend(op)
    if backend != "qutip" and isoper(op) is False:
        raise ValueError("For all backends except qutip, partial trace currently only works on operators.") 
    op = tidyup(op, dims = dims)
    return op.ptrace(sel)

def eigenstates(op, isherm=True):
    """Eigenstates of a quantum object.
    
    :param op: Operator.
    :param isherm: If True, the operator is assumed to be Hermitian.
    :return: eigenvalues, eigenvectors.
    
    """
    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    if backend == 'numpy':
        return op.eigenstates(isherm)
    else:
        return op.eigenstates()

def permute(op, order : list):
    """ Permutes members of a quantum object. 
    
    :param op: Operator.
    :return: The operator with permuted members.
    
    """
    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    return op.permute(order)

###################### Math / Generic #######################

def expect(oper,state):
    """ Expectation value for an operator-state pair.

    :param oper: Operator.
    :param state: State.

    :return: Expectation value.
    
    """
    backend = backends.get_backend(state)
    return getattr(getattr(backends,backend), 'expect')(oper,state)

def tensor(operators : list):
    """ Tensor product of a series of operators.
    
    :param operators: A list of operators
    :return: Tensor product of the operators.
    """    

    backend = backends.get_backend(operators[0])
    for idx, item in enumerate(operators):
        operators[idx] = tidyup(item, method = backend)
    return getattr(getattr(backends,backend), 'tensor')(operators)


def directsum(*args):
    """ Direct sum of a series of operators. Empty, square operators may also be specified as 
    integers.
    
    :param *args: A series of operators, provided as a single list or multiple arguments.
    :return: Direct product of the operators.
    """       
    backend = backends.get_backend(operator)
    return getattr(getattr(backends,backend), 'directsum')(*args)

def ddrop(operator, idx : Union[int, list]):
    """ Drops dimension(s) from an operator (i.e. a quantum object).
    
    :param operator: A quantum object.
    :param idx: The dimension that will be dropped.
    
    :return: The operator without the dropped dimensions.
    """

    backend = backends.get_backend(operator)
    return getattr(getattr(backends,backend), 'ddrop')(operator, idx)

def dpermute(operator, order : list):
    """ Permutes dimensions of an operator.
    
    :param operator: A quantum object.
    :param order: A list specifying the new ordering of dimensions.

    :return: The operator with permuted dimensions. 
    """

    backend = backends.get_backend(operator)
    return getattr(getattr(backends,backend), 'dpermute')(operator, order)

def block_diagonal(L_list):
    """ Constructs a block-diagonal matrix from a list of operators.

    :param list L_list: A list of operators.

    :return: A block diagonal matrix, each block is an entry of L_list.
    """
    backend = backends.get_backend(L_list[0])
    L_list = [tidyup(L, method = backend) for L in L_list]
    return getattr(getattr(backends,backend), 'block_diagonal')(L_list)

###################### HILBERT SPACE  ########################

def isbra(state):
    """ Indicates if the quantum object is a bra.

    :param state: Quantum object.

    :return bool:
    
    """
    backend = backends.get_backend(state)
    return getattr(getattr(backends,backend), 'isbra')(state)

def isket(state):
    """ Indicates if the quantum object is a ket.

    :param state: Quantum object.

    :return bool:
    
    """
    backend = backends.get_backend(state)
    return getattr(getattr(backends,backend), 'isket')(state)

def isoper(state):
    """ Indicates if the quantum object is an operator.

    :param state: Quantum object.
    
    :return bool:
    
    """
    backend = backends.get_backend(state)
    return getattr(getattr(backends,backend), 'isoper')(state)

def ket2dm(ket):
    """ Constructs an density matrix from a state vector with an outer product. 
    
    :param ket: A fock state vector (ket or bra)

    :return: Density matrix.
    """

    backend = backends.get_backend(ket)
    ket = tidyup(ket, method = backend) 
    return getattr(getattr(backends,backend), 'ket2dm')(ket)

def dm2ket(dm, tol = 1e-6):
    """ Constructs a state vector (ket) from a density matrix. Does only work for pure states and normalized density matrices.
    
    :param dm: A density matrix.
    :parm tol: Tolerance used when probing whether the input density matrix is a pure state and normalized.

    :return: State vector (ket).
    """
    backend = backends.get_backend(dm)
    dm = tidyup(dm, method = backend)
    return getattr(getattr(backends,backend), 'dm2ket')(dm,tol)

def ket(N : int, m  : int ,method='qutip'):
    """ Pure state vector (ket) for the mth basis state of an N-dimensional Hilbert space.

    :param int N: Number of fock states in Hilbert space.
    :param int m: The index of the desired state, ranging from 0-(N-1)
    :param str method: The desired backend. Can be 'numpy',  'qutip','sparse' or 'sympy'.

    :return: The pure state vector (ket)
    """
    return getattr(getattr(backends,method), 'ket')(N,m)

def bra(N : int, m : int, method= 'qutip'):
    """ Pure state vector (bra) for the mth basis state of an N-dimensional Hilbert space.

    :param int N: Number of fock states in Hilbert space.
    :param int m: The index of the desired state, ranging from 0-(N-1)
    :param str method: The desired backend. Can be 'numpy', 'qutip','sparse' or 'sympy'.

    :return: The pure state vector (bra)
    
    """
    return getattr(getattr(backends,method), 'bra')(N, m)

def jmat(j , op_spec : str, method='qutip'):
    """Spin operator for a spin j.

    :param j:  A non-negative integer or half-integer specifying the spin
    :param str op_spec: A string specifying the desired operator, can be 'x', 'y', 'z','+','-'

    :return: The spin operator.
    """

    return getattr(getattr(backends,method), 'jmat')(j,op_spec)


def identity(N : int , dims = None , method ='qutip'):
    """ Identiy Matrix of an N-dimensional Hilbert space. 

    :param int N: Number of Fock states in Hilbert space.
    :param list dims: Structure of the Hilbert space, if None or omitted dims are set as [N].
    :param str method: The desired backend. Can be 'numpy',  'qutip','sparse' or 'sympy'.

    :return: The identity matrix.
    
    """
    return getattr(getattr(backends,method), 'identity')(N,dims)

def diags(v , k : int = 0, method='qutip', dims = None):
    """ Constructs an operator from a diagonal or extracts a diagonal.
    
    :param v: A sequence of elements to be placed along the selected diagonal or a quantum object from which diagonal is extracted.
    :param int k: The selected diagonal
    :param list dims: Structure of the Hilbert space, if None or omitted dims are set as [N].
    :param str method: The desired backend. Can be 'numpy',  'qutip','sparse' or 'sympy'.

    :return: Resulting matrix.
    """
    # If diags was called to extract a diagonal, use native backend.
    if len(_np.shape(v)) == 2: 
        backend = backends.get_backend(v)
        if backend != method:
            warnings.warn("Input data type suggests other backend than specified method that will be ignored.")
        method = backend
    return getattr(getattr(backends, method), 'diags')(v, k, dims)


################ Liouville Space ###########################ß

def dm2vec(op):
    """ Constructs a liouville space vector from a Hilbert space operator.
    
    :param op: An operator in Hilbert space.
    
    :return: State vector in Liouville space.
    """
    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    return getattr(getattr(backends,backend), 'dm2vec')(op)

def vec2dm(vec):
    """ Constructs  Hilbert space operator from a liouville space vector.
    
    :param op: State vector in Liouville space.An operator in Hilbert space.
    
    :return: An operator in Hilbert space.
    """
    backend = backends.get_backend(vec)
    vec = tidyup(vec, method = backend)
    return getattr(getattr(backends,backend), 'vec2dm')(vec)

def issuper(state):
    """ Indicates if the quantum object is a superoperator.

    :param state: Quantum object.

    :return bool:
    """
    backend = backends.get_backend(state)
    return getattr(getattr(backends,backend), 'issuper')(state)

def spost(op):
    """ Computes superoperator by post-multiplication by op.
    
    :param op: Operator.
    :return: Superoperator 1 \\otimes op. """

    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    return getattr(getattr(backends,backend), 'spost')(op)

def spre(op):
    """ Computes superoperator by pre-multiplication by op.
    
    :param op: Operator.
    :return: Superoperator  op.T \\otimes 1."""

    backend = backends.get_backend(op)
    op = tidyup(op, method = backend)
    return getattr(getattr(backends,backend), 'spre')(op)

def lindbladian(a, b=None):
    """ Lindbladian superoperator for a pair of collapse operators a,b, or a single collapse operator if only a is being provided.
    
    :param a: Collapse operator.
    :param b: Collapse operator, None per default.
    :return: Lindbladian superoperator. """
    
    backend = backends.get_backend(a)
    a = tidyup(a, method = backend)
    if b is not None:
        b = tidyup(b, method = backend)
    return getattr(getattr(backends,backend), 'lindbladian')(a, b)

def liouvillian(H, c_ops : list):
    """Liouvillian superoperator from a hamiltonian and a list of collapse operators.

    :param H: Hamiltonian operator. 
    :param list c_ops: List of collapse operators.

    :return: Liouville superoperator.
    """
    if H is not None:
        backend = backends.get_backend(H)
        H = tidyup(H, method = backend)
    else:
        try:
            backend = backends.get_backend(c_ops[0])
        except:
            raise ValueError("If no Hamiltonian is provided, list of collapse operators cannot be empty.")
    for idx, item in enumerate(c_ops):
        c_ops[idx] = tidyup(item, method = backend)
    return getattr(getattr(backends,backend), 'liouvillian')(H,c_ops)

def applySuperoperator(L,rho):
    """Apply superoperator to a state.

    :param L: Superoperator (Liouvillian).
    :param rho: State (density matrix).

    :return: Propagator.
    """
    backend = backends.get_backend(rho)
    return getattr(getattr(backends,backend), 'applySuperoperator')(L,rho)


################ Fokker Planck Space  #########################
def fp2dm(v,is_hilbert,dims=None,weights=None):
    """Converts a state vector in the Fokker-Planck space to a density matrix.

    In the Fokker-Planck space, the state vector is a vector of concatenated state vectors in Hilber space (pure quantum states) or vectorized density matrices (mixed quantum states).
    
    :param v: State vector in the Fokker-Planck space.
    :param bool is_hilbert: If True, the state vector is defined in the Hilbert space. If False, the state vector is defined in the Liouville space.
    :param list dims: dims to set. If None, dims are set as the Hilbert / Liouville part.
    :param weights: Weights for the density matrix. If None, the weights are set as 1.
    """
    backend = backends.get_backend(v)
    v = tidyup(v, method = backend)
    return getattr(getattr(backends,backend), 'fp2dm')(v,is_hilbert,dims,weights)

def dm2fp(rho,dof,try_ket=False):
    """Copies a density matrix to a state vector in the Fokker-Planck space.

    In the Fokker-Planck space, the state vector is a vector of concatenated state vectors in Hilbert space (pure quantum states) or vectorized density matrices (mixed quantum states).

    :param rho: Density matrix or state vector
    :param int dof: Number of degrees of freedom. This is the number of copies.
    :param bool try_ket: If True, the function tries to convert the density matrix to a state vector in Hilbert space.
    """

    backend = backends.get_backend(rho)
    rho = tidyup(rho, method = backend)
    return getattr(getattr(backends,backend), 'dm2fp')(rho,dof,try_ket)

def ket2fp(ket,dof):
    """Copies a state vector to a state vector in the Fokker-Planck space.

    In the Fokker-Planck space, the state vector is a vector of concatenated state vectors in Hilbert space (pure quantum states) or vectorized density matrices (mixed quantum states).

    :param ket: State vector
    :param int dof: Number of degrees of freedom. This is the number of copies.
    """

    backend = backends.get_backend(ket)
    ket = tidyup(ket, method = backend)
    return getattr(getattr(backends,backend), 'ket2fp')(ket,dof)