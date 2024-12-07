import simos
import pytest
import qutip as qu
import sympy as sp
import numpy as np

@pytest.fixture
def NVops():
    Se = qu.qeye(3)
    Ie = qu.qeye(2)
    Sx = qu.jmat(1,'x')
    Sy = qu.jmat(1,'y')
    Sz = qu.jmat(1,'z')
    Sx = qu.tensor(Sx,Ie)
    Sy = qu.tensor(Sy,Ie)
    Sz = qu.tensor(Sz,Ie)
    Ix = qu.jmat(1/2,'x')
    Iy = qu.jmat(1/2,'y')
    Iz = qu.jmat(1/2,'z')
    Ix = qu.tensor(Se,Ix)
    Iy = qu.tensor(Se,Iy)
    Iz = qu.tensor(Se,Iz)

    NVm = qu.basis(3,0)
    NV0 = qu.basis(3,1)
    NVp = qu.basis(3,2)
    Im = qu.basis(2,0)
    Ip = qu.basis(2,1)
    NVm = qu.tensor(NVm*NVm.dag(),Ie)
    NV0 = qu.tensor(NV0*NV0.dag(),Ie)
    NVp = qu.tensor(NVp*NVp.dag(),Ie)
    Im = qu.tensor(Se,Im*Im.dag())
    Ip = qu.tensor(Se,Ip*Ip.dag())

    identitity = qu.tensor(Se,Ie)

    return Sx, Sy, Sz, Ix, Iy, Iz, NVm, NV0, NVp, Im, Ip, identitity

def convert_to_numpy(*args):
    return [simos.backends.numpy.QNdarray(arg.full(),dims=arg.dims) for arg in args]

def convert_to_sympy(*args):
    # ensure that we have rationals
    def _handle_floats(arg_np):
        arg_signs = (-(arg_np < 0).astype(int))*2 + 1
        arg_signs = sp.Matrix(arg_signs)
        arg_prep = (2*arg_np)**2
        arg_prep = arg_prep.astype(int)
        arg_sympy = sp.Matrix(arg_prep)
        arg_sympy = arg_sympy.applyfunc(lambda x: sp.sqrt(x))/2
        arg_sympy = sp.matrices.dense.matrix_multiply_elementwise(arg_sympy,arg_signs)
        return arg_sympy
        
    new_args = []
    for arg in args:
        arg_real = np.real(arg.full())
        arg_imag = np.imag(arg.full())
        arg_sympy = _handle_floats(arg_real) + sp.I*_handle_floats(arg_imag)       
        new_args.append(simos.backends.sympy.QSMatrix(arg_sympy,dims=arg.dims))
    return new_args

def convert_to_sparse(*args):
    return [simos.backends.sparse.Qcsc_matrix(arg.full(),dims=arg.dims) for arg in args]

@pytest.mark.parametrize("method", ['qutip','numpy','sympy','sparse'])
class TestSpinsystem:    

    @pytest.mark.parametrize("val", np.arange(0.5,10,0.5))
    def test_single_spin(self, method, val):
        systemdef = []
        systemdef.append({'name':'S', 'val':val})
        s = simos.create_system(systemdef,method=method)
        assert s
        return None
    
    def test_nv_spinsystem(self, NVops,method):
        systemdef = []
        systemdef.append({'name':'S', 'val':1,'type':'NV'})
        systemdef.append({'name':'I', 'val':1/2})
        s = simos.create_system(systemdef,method=method)

        Sx, Sy, Sz, Ix, Iy, Iz, NVm, NV0, NVp, Im, Ip, identitity = NVops

        if method == 'numpy' or method == 'numba':
            Sx, Sy, Sz, Ix, Iy, Iz, NVm, NV0, NVp, Im, Ip, identitity = convert_to_numpy(Sx, Sy, Sz, Ix, Iy, Iz, NVm, NV0, NVp, Im, Ip, identitity)
        elif method == 'sympy':
            #print(Sx)
            Sx, Sy, Sz, Ix, Iy, Iz, NVm, NV0, NVp, Im, Ip, identitity = convert_to_sympy(Sx, Sy, Sz, Ix, Iy, Iz, NVm, NV0, NVp, Im, Ip, identitity)
            #print(Sx)
        elif method == 'sparse':
            Sx, Sy, Sz, Ix, Iy, Iz, NVm, NV0, NVp, Im, Ip, identitity = convert_to_sparse(Sx, Sy, Sz, Ix, Iy, Iz, NVm, NV0, NVp, Im, Ip, identitity)
        
        def _all(test):
            if method == 'qutip':
                return test
            elif method == 'numpy' or method == 'numba':
                return test.all()
            elif method == 'sympy':
                return test
            elif method == 'sparse':
                return test.toarray().all()
            
        # Check identity
        # assert _all(s.id == identitity)    

        # Check spinoperators
        assert _all((s.Sx == Sx))
        assert _all((s.Sy == Sy))
        assert _all((s.Sz == Sz))
        assert _all((s.Ix == Ix))
        assert _all((s.Iy == Iy))
        assert _all((s.Iz == Iz))

        # Check projectors
        assert _all((s.Sp[-1] == NVm))
        assert _all((s.Sp[0] == NV0))
        assert _all((s.Sp[1] == NVp))
        assert _all((s.Sp[-1.0] == NVm))
        assert _all((s.Sp[0.0] == NV0))
        assert _all((s.Sp[1.0] == NVp))
        assert _all((s.Ip[-1/2] == Im))
        assert _all((s.Ip[1/2] == Ip))
        assert _all((s.Ip[-0.5] == Im))
        assert _all((s.Ip[0.5] == Ip))

        return None
    
    def test_nv_level_model(self,method):
        S = {'name':'S', 'val':1}
        GS = {'name':'GS', 'val':0}
        ES = {'name':'ES','val':0}
        SS = {'name':'SS','val':0}
        s = simos.create_system(([S,(GS,ES)],SS),method=method)
        assert s
        return None