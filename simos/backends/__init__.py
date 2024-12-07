import importlib.util

# Numpy backend, always available
from . import numpy as numpy
from .. import constants
AVAILABLE = ['numpy']

######################################
# Conditional import of backends
######################################

if importlib.util.find_spec('qutip') is not None:
    from . import qutip as qutip
    AVAILABLE.append('qutip')

if importlib.util.find_spec('sympy') is not None:
    from . import sympy as sympy
    AVAILABLE.append('sympy')

if importlib.util.find_spec('scipy.sparse') is not None:
    from . import sparse as sparse
    AVAILABLE.append('sparse')

######################################
# Backend helper
######################################
def get_backend(obj):
    """ Determine which backend should be used based on object type. 
    :param obj: Object for which backend should be determined."""

    if isinstance(obj,(list, tuple)):
        # check if all elements are of the same type
        for item in obj:
            if get_backend(item) != get_backend(obj[0]):
                raise ValueError("All elements in the list must be of the same type")
        # Note for future - one might want to add dictionary checking here. Currently not here for performance reasons (not required and thus omitted.)
        return get_backend(obj[0])
    if 'numpy' in AVAILABLE:
        if isinstance(obj, numpy._np.ndarray):
            return 'numpy'
    if 'sparse' in AVAILABLE:
        if isinstance(obj, (sparse.Qcsc_matrix, sparse._sp.spmatrix, sparse._sp.sparray)):
            return 'sparse'
    if 'qutip' in AVAILABLE:
        if isinstance(obj, qutip._qu.Qobj):
            return 'qutip'    
    if 'sympy' in AVAILABLE:
        if isinstance(obj, (sympy._sp.Matrix, sympy._sp.ImmutableDenseMatrix, sympy._sp.MutableDenseMatrix, sympy._sp.ImmutableSparseMatrix, sympy._sp.MutableSparseMatrix, 
                            sympy._sp.Array, sympy._sp.ImmutableDenseNDimArray, sympy._sp.MutableDenseNDimArray, sympy._sp.ImmutableSparseNDimArray, sympy._sp.MutableSparseNDimArray, 
                            sympy._sp.core.symbol.Symbol, sympy._sp.Quaternion, sympy._sp.core.mul.Mul)):
                return 'sympy'
    # default to numpy
    return 'numpy'


###################################
# Symbolic vs numeric math methods
###################################
def get_calcmethod(function, mode):
    """ Handle backend specifity for simple math methods such as, for example, taking the square-root or exponent.
        Symbolic vs non symbolic: symbolic functions are returned for sympy backend, meant to operate on sympys MutableDenseNDimArray,
        non-symbolic functions are returend for all other backend, meant to operate on numpys Ndarray"""
    calcmethods = {}
    # eigenvalues/vectors
    calcmethods["eigenstates"] = {}
    calcmethods["eigenstates"]["symbolic"] = sympy.matrix_eigenstates
    calcmethods["eigenstates"]["numeric"] =  numpy.matrix_eigenstates 
    # rotation object / euler angles ; perform rotation 
    calcmethods["fromeuler"] = {} # build rotation object from euler angles
    calcmethods["fromeuler"]["symbolic"] = lambda x,y : sympy._sp.Quaternion.from_euler(y, seq = x)
    calcmethods["fromeuler"]["numeric"] = sparse._sc.spatial.transform.Rotation.from_euler
    calcmethods["frommatrix"] = {} # build rotation object from matrix 
    calcmethods["frommatrix"]["symbolic"] = sympy._sp.Quaternion.from_rotation_matrix
    calcmethods["frommatrix"]["numeric"] = sparse._sc.spatial.transform.Rotation.from_matrix
    calcmethods["toeuler"] = {} # get euler angles from rotation object
    calcmethods["toeuler"]["symbolic"] = lambda x, y: x.to_euler(y)
    calcmethods["toeuler"]["numeric"] = lambda x, y: x.as_euler(y)
    calcmethods["rotate"] = {}
    calcmethods["rotate"]["symbolic"] = sympy.rotate
    calcmethods["rotate"]["numeric"] = numpy.rotate
    # sqrt 
    calcmethods["sqrt"]= {}
    calcmethods["sqrt"]["symbolic"] = sympy.sqrt_fun 
    calcmethods["sqrt"]["numeric"]  = numpy._np.sqrt  
    # power 
    calcmethods["pow"]= {}
    calcmethods["pow"]["symbolic"] = sympy.pow_fun 
    calcmethods["pow"]["numeric"]  = numpy.pow_fun  
    # exp
    calcmethods["exp"]= {}
    calcmethods["exp"]["symbolic"] = sympy.exp_fun 
    calcmethods["exp"]["numeric"]  = numpy._np.exp  
    # sine 
    calcmethods["sin"]= {}
    calcmethods["sin"]["symbolic"] = sympy.sin_fun 
    calcmethods["sin"]["numeric"]  = numpy._np.sin  
    # cosine 
    calcmethods["cos"]= {}
    calcmethods["cos"]["symbolic"] = sympy.cos_fun
    calcmethods["cos"]["numeric"]  = numpy._np.cos  
    # atan2 
    calcmethods["arctan2"]= {}
    calcmethods["arctan2"]["symbolic"] = sympy.arctan2_fun 
    calcmethods["arctan2"]["numeric"]  = numpy._np.arctan2  
    # elemntwise multiply of 2 arrays with same shape 
    calcmethods["multiply"]= {}
    calcmethods["multiply"]["symbolic"] = sympy.multelem_fun 
    calcmethods["multiply"]["numeric"]  = numpy._np.multiply 
    # matrix multiplication
    calcmethods["matmul"] = {}
    calcmethods["matmul"]["symbolic"] = sympy.matmul_fun
    calcmethods["matmul"]["numeric"] = numpy._np.matmul
    # make to array, not an actual mathematical method, here for convenience
    calcmethods["array"]={}
    calcmethods["array"]["symbolic"] = sympy._sp.MutableDenseNDimArray
    calcmethods["array"]["numeric"] = numpy._np.array
    # symbols
    calcmethods["hbar"]={}
    calcmethods["hbar"]["symbolic"] = sympy._sp.physics.quantum.constants.hbar
    calcmethods["hbar"]["numeric"] = constants.hbar
    # symbols
    calcmethods["mu0"]={}
    calcmethods["mu0"]["symbolic"] = sympy._sp.Symbol("mu0")
    calcmethods["mu0"]["numeric"] = constants.mu_0
    # symbols
    calcmethods["pi"]={}
    calcmethods["pi"]["symbolic"] = sympy._sp.pi
    calcmethods["pi"]["numeric"] = numpy._np.pi
    # symbols
    calcmethods["I"]={}
    calcmethods["I"]["symbolic"] = sympy._sp.I
    calcmethods["I"]["numeric"] = 1j  
    # trace
    calcmethods["trace"] = {}
    calcmethods["trace"]["symbolic"] = lambda x: sympy._sp.Trace(sympy._sp.Matrix(x)).doit()
    calcmethods["trace"]["numeric"] = numpy._np.trace
    # cg coefficients 
    calcmethods["cg"] = {}
    calcmethods["cg"]["symbolic"] = sympy.cg
    calcmethods["cg"]["numeric"] =  numpy.cg
    # Wigner D matrix
    calcmethods["WignerD"] = {}
    calcmethods["WignerD"]["symbolic"] =  sympy._sp.physics.wigner.wigner_d
    calcmethods["WignerD"]["numeric"] =  numpy.WignerD
    return calcmethods[function][mode]