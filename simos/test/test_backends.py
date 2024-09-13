import pytest
import simos
import numpy as np
import qutip as qu

def create(method,data,dims=None):
    if method == 'qutip':
        return qu.Qobj(data,dims=None)
    elif method == 'numpy':
        return simos.backends.numpy.QNdarray(data,dims=dims)
    elif method == 'sympy':
        return simos.backends.sympy.QSMatrix(data,dims=dims)
    elif method == 'sparse':
        return simos.backends.sparse.Qcsc_matrix(data,dims=dims)
    else:
        raise ValueError("Invalid method")  

def get_data(method,obj):
    if method == 'qutip':
        return obj.full()
    elif method == 'numpy':
        return obj
    elif method == 'sympy':
        return np.array(obj).astype(np.complex128)
    elif method == 'sparse':
        return obj.toarray()

@pytest.mark.parametrize("method", ['qutip','numpy','sympy','sparse'])
class TestDataStructure:
    @pytest.mark.parametrize("Ndim", 2**np.arange(1,10))
    def test_create(self,Ndim,method):
        data = np.random.randint(-10,10,(Ndim,Ndim))
        dims = [Ndim,Ndim]
        obj = create(method,data,dims)
        assert obj.dims == dims
        assert np.allclose(get_data(method,obj),data)

    @pytest.mark.parametrize("funname", ['ket','bra', 'isket', 'isbra', 'ket2dm', 'dm2ket', 'operator_to_vector', 'vector_to_operator', 'jmat', 'identity', 'diags', 'block_diagonal', 'expect', 'tensor', 'spre', 'spost', 'liouville', 'liouvillian', 'applySuperoperator', 'differentiation_matrix', 'concatenate', 'split', 'tidyup', 'data', 'changeelement'])
    def test_common_functions(self,method,funname):
        getattr(getattr(simos.backends,method), funname)

@pytest.mark.parametrize("j1, m1, j2, m2, J, M", [
    (1, 0, 1, 0, 1, 0),
    (0.5, -0.5, 0.5, 0.5, 1, 1),
    (1, -1, 1, 1, 2, 0),
    # Add more test cases as needed
])
def test_cg(j1, m1, j2, m2, J, M):
    # Calculate the Clebsch-Gordan coefficient using the numpy backend
    cg_nb = simos.backends.numpy.cg(j1, m1, j2, m2, J, M)

    # Calculate the Clebsch-Gordan coefficient using qutip
    cg_qu = qu.clebsch(j1, j2, J, m1, m2, M)
    
    # Check that the results are the same
    assert np.isclose(cg_nb, cg_qu), f"Expected {cg_qu}, but got {cg_nb}"

def test_eigen_consistency():
    system_def = []
    system_def.append({'name':'S','val':1,'type':'NV'})
    system_def.append({'name':'I','val':1/1,'type':'19F'})


    s = simos.System(system_def,method='qutip')
    H = simos.D*s.Sz**2 + simos.ye*s.Sz

    vals,vectors = H.eigenstates()
    
    s = simos.System(system_def,method='numpy')
    H = simos.D*s.Sz**2 + simos.ye*s.Sz
    valsNP,vectorsNP = H.eigenstates()

    for i in range(len(vals)):
        assert np.allclose(vals[i],valsNP[i])
        assert np.allclose(vectors[i].full(),vectorsNP[i])
    
    s = simos.System(system_def,method='sympy')
    H = simos.D*s.Sz**2 + simos.ye*s.Sz
    valsSP,vectorsSP = H.eigenstates()

    for i in range(len(vals)):
        assert np.allclose(vals[i],complex(valsSP[i]))
        assert np.allclose(vectors[i].full(),simos.data(vectorsSP[i]).astype(complex))