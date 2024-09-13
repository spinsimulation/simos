from numba.extending import get_cython_function_address
from numba import vectorize, njit,jit, complex128, complex64
import ctypes
import numpy as np
import qutip as qu
from numba.extending import typeof_impl

_PTR  = ctypes.POINTER

_dble = ctypes.c_double
_sngl = ctypes.c_float
_char = ctypes.c_char
_int  = ctypes.c_int

_ptr_select = ctypes.c_voidp
_ptr_dble = _PTR(_dble)
_ptr_sngl = _PTR(_sngl)
_ptr_char = _PTR(_char)
_ptr_int  = _PTR(_int)


# signature is:
# void zgees(
#  char *jobvs, 
#  char *sort, 
#  zselect1 *select, 
#  int *n, 
#  z *a, 
#  int *lda,
#  int *sdim, 
#  z *w, 
#  z *vs, 
#  int *ldvs, 
#  z *work, 
#  int *lwork, 
#  d *rwork, 
#  bint *bwork, 
#  int *info
# )
# bind to the real space variant of the function
addr = get_cython_function_address('scipy.linalg.cython_lapack', 'zgees')
functype = ctypes.CFUNCTYPE(None,
                            _ptr_int, # JOBVS
                            _ptr_int, # SORT
                            _ptr_dble, # SELECT
                            _ptr_int, # N
                            _ptr_dble, # A
                            _ptr_int, # LDA
                            _ptr_int, # SDIM
                            _ptr_dble, # W
                            _ptr_dble, # VS
                            _ptr_int, # LDVS
                            _ptr_dble, # WORK
                            _ptr_int, # LWORK
                            _ptr_dble, # RWORK
                            _ptr_int, # BWORK
                            _ptr_int, # INFO
                            )
zgees_fn = functype(addr)

# call signature is:
# void zgees(
#  char *jobvs, 
#  char *sort, 
#  zselect1 *select, 
#  int *n, 
#  z *a, 
#  int *lda,
#  int *sdim, 
#  z *w, 
#  z *vs, 
#  int *ldvs, 
#  z *work, 
#  int *lwork, 
#  d *rwork, 
#  bint *bwork, 
#  int *info
# )
# bind to the real space variant of the function
addr2 = get_cython_function_address('scipy.linalg.cython_lapack', 'cgees')
functype = ctypes.CFUNCTYPE(None,
                            _ptr_int, # JOBVS
                            _ptr_int, # SORT
                            _ptr_sngl, # SELECT
                            _ptr_int, # N
                            _ptr_sngl, # A
                            _ptr_int, # LDA
                            _ptr_int, # SDIM
                            _ptr_sngl, # W
                            _ptr_sngl, # VS
                            _ptr_int, # LDVS
                            _ptr_sngl, # WORK
                            _ptr_int, # LWORK
                            _ptr_sngl, # RWORK
                            _ptr_int, # BWORK
                            _ptr_int, # INFO
                            )
cgees_fn = functype(addr2)


@njit
def numba_zgees(x):
    JOBVS = np.array([86], np.int32) # ord('V') to get Schur vectors
    SORT = np.array([78], np.int32)  # ord('N') to not sort
    SELECT = np.empty(1, np.float64)
    _M, _N = x.shape
    N = np.array(_N, np.int32)
    A = x     # in & out
    LDA = np.array(_N, np.int32)
    SDIM = np.empty(1, np.int32) # out
    W = np.empty(_N, dtype=np.complex128) #out
    VS = np.empty((_N, _N), dtype=np.complex128) #out
    LDVS = np.array(_N, np.int32)
    WORK = np.empty((1,), dtype=np.complex128) #out
    LWORK = np.array(-1, dtype=np.int32)
    RWORK = np.empty(_N, dtype=np.float64)
    BWORK = np.empty(_N, dtype=np.int32)
    INFO = np.empty(1, dtype=np.int32)

    def check_info(info):
        if info[0] != 0:
            raise RuntimeError("INFO indicates problem")

    zgees_fn(JOBVS.ctypes,
             SORT.ctypes,
             SELECT.ctypes,
             N.ctypes,
             A.view(np.float64).ctypes,
             LDA.ctypes,
             SDIM.ctypes,
             W.view(np.float64).ctypes,
             VS.view(np.float64).ctypes,
             LDVS.ctypes,
             WORK.view(np.float64).ctypes,
             LWORK.ctypes,
             RWORK.ctypes,
             BWORK.ctypes,
             INFO.ctypes)
    check_info(INFO)
    #print("Calculated workspace size as", WORK[0])
    WS_SIZE = np.int32(WORK[0].real)
    LWORK = np.array(WS_SIZE, np.int32)
    WORK = np.empty(WS_SIZE, dtype=np.complex128)
    zgees_fn(JOBVS.ctypes,
             SORT.ctypes,
             SELECT.ctypes,
             N.ctypes,
             A.view(np.float64).ctypes,
             LDA.ctypes,
             SDIM.ctypes,
             W.view(np.float64).ctypes,
             VS.view(np.float64).ctypes,
             LDVS.ctypes,
             WORK.view(np.float64).ctypes,
             LWORK.ctypes,
             RWORK.ctypes,
             BWORK.ctypes,
             INFO.ctypes)
    check_info(INFO)
    return A, VS.T # .T the schur vectors, Fortran vs C order


@njit
def numba_cgees(x):
    JOBVS = np.array([86], np.int32) # ord('V') to get Schur vectors
    SORT = np.array([78], np.int32)  # ord('N') to not sort
    SELECT = np.empty(1, np.float32)
    _M, _N = x.shape
    N = np.array(_N, np.int32)
    A = x     # in & out
    LDA = np.array(_N, np.int32)
    SDIM = np.empty(1, np.int32) # out
    W = np.empty(_N, dtype=np.complex64) #out
    VS = np.empty((_N, _N), dtype=np.complex64) #out
    LDVS = np.array(_N, np.int32)
    WORK = np.empty((1,), dtype=np.complex64) #out
    LWORK = np.array(-1, dtype=np.int32)
    RWORK = np.empty(_N, dtype=np.float32)
    BWORK = np.empty(_N, dtype=np.int32)
    INFO = np.empty(1, dtype=np.int32)

    def check_info(info):
        if info[0] != 0:
            raise RuntimeError("Function return indicates problem")

    cgees_fn(JOBVS.ctypes,
             SORT.ctypes,
             SELECT.ctypes,
             N.ctypes,
             A.view(np.float32).ctypes,
             LDA.ctypes,
             SDIM.ctypes,
             W.view(np.float32).ctypes,
             VS.view(np.float32).ctypes,
             LDVS.ctypes,
             WORK.view(np.float32).ctypes,
             LWORK.ctypes,
             RWORK.ctypes,
             BWORK.ctypes,
             INFO.ctypes)
    check_info(INFO)
    #print("Calculated workspace size as", WORK[0])
    WS_SIZE = np.int32(WORK[0].real)
    LWORK = np.array(WS_SIZE, np.int32)
    WORK = np.empty(WS_SIZE, dtype=np.complex64)
    cgees_fn(JOBVS.ctypes,
             SORT.ctypes,
             SELECT.ctypes,
             N.ctypes,
             A.view(np.float32).ctypes,
             LDA.ctypes,
             SDIM.ctypes,
             W.view(np.float32).ctypes,
             VS.view(np.float32).ctypes,
             LDVS.ctypes,
             WORK.view(np.float32).ctypes,
             LWORK.ctypes,
             RWORK.ctypes,
             BWORK.ctypes,
             INFO.ctypes)
    check_info(INFO)
    return A, VS.T # .T the schur vectors, Fortran vs C order

@njit
def expm2(a):
    # We have to convert the array to fortran order and return it in c order which is the default in numpy
    mat_in = np.asfortranarray(np.copy(a).astype(complex128))
    t, vs = numba_zgees(mat_in)
    d = np.diag(t)
    d = np.exp(d)
    #print(d.dtype)
    exp = (vs @ (np.diag(d)).astype(complex128) @ vs.conj().T)
    return np.ascontiguousarray(exp)

@njit
def expm3(a):
    # We have to convert the array to fortran order and return it in c order which is the default in numpy
    mat_in = np.asfortranarray(np.copy(a).astype(complex64))
    t, vs = numba_cgees(mat_in)
    d = np.diag(t)
    d = np.exp(d)
    #print(d.dtype)
    exp = (vs @ (np.diag(d)).astype(complex64) @ vs.conj().T)
    return np.ascontiguousarray(exp)

def expm(a):
    use_qutip = False
    if isinstance(a, qu.Qobj):
        qutip_a = a
        a = a.data.todense()
        use_qutip = True
    exp = expm2(a)
    if use_qutip:
        return qu.Qobj(exp, dims=qutip_a.dims)
    else:
        return exp

@njit
def resize(a, new_size):
    new = np.zeros(new_size)
    idx = 0
    while True:
        newidx = idx + len(a)
        if newidx > new_size:
            new[idx:] = a[:new_size-newidx]
            break
        new[idx:newidx] = a
        idx = newidx
    return new
