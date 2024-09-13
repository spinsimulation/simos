import simos as sos
import numpy as np
import qutip as qu
import pytest

@pytest.mark.parametrize("method", ["sympy", "numpy","sparse"])
def test_operator_to_vector_qutip_consistency(method):
    test = np.array([[1,2],[3,4]])
    qrho = qu.Qobj(test)
    
    rho = sos.tidyup(test, method=method)
    rhodata = sos.data(rho).astype(float)
    assert np.allclose(rhodata, qrho.full())

    qvec = qu.operator_to_vector(qrho)
    vec = sos.operator_to_vector(rho)
    vecdata = sos.data(vec).astype(float)
    assert np.allclose(vecdata, qvec.full())

    # Check that the inverse operation is consistent
    qrho2 = qu.vector_to_operator(qvec)
    rho2 = sos.vector_to_operator(vec)
    rho2data = sos.data(rho2).astype(float)
    assert np.allclose(rho2data, qrho2.full())
    assert np.allclose(rho2data, qrho.full())


@pytest.mark.parametrize("method", ["sympy", "numpy","sparse"])
def test_spre_qutip_consistency(method):
    test = np.array([[1,2],[3,4]])
    qrho = qu.Qobj(test)
    qs = qu.spre(qrho)

    rho = sos.tidyup(test, method=method)
    s = sos.spre(rho)
    sdata = sos.data(s).astype(float)
    assert np.allclose(sdata, qs.full())

@pytest.mark.parametrize("method", ["sympy", "numpy","sparse"])
def test_spost_qutip_consistency(method):
    test = np.array([[1,2],[3,4]])
    qrho = qu.Qobj(test)
    qs = qu.spost(qrho)

    rho = sos.tidyup(test, method=method)
    s = sos.spost(rho)
    sdata = sos.data(s).astype(float)
    assert np.allclose(sdata, qs.full())

@pytest.mark.parametrize("method", ["sympy", "numpy","sparse"])
def test_sprepost_qutip_consistency(method):
    test = np.array([[1,2],[3,4]])
    qrho = qu.Qobj(test)
    qs = qu.spre(qrho) - qu.spost(qrho)

    sprepost = getattr(sos.backends,method).sprepost
    rho = sos.tidyup(test, method=method)
    s = sprepost(rho)
    sdata = sos.data(s).astype(complex)
    assert np.allclose(sdata, qs.full())

@pytest.mark.parametrize("method", ["sympy", "numpy","sparse"])
def test_lindblad_qutip_consistency(method):
    test = np.array([[1,2],[3,4]])
    qrho = qu.Qobj(test)
    qL = qu.lindblad_dissipator(qrho)

    rho = sos.tidyup(test, method=method)
    L = sos.lindbladian(rho)
    Ldata = sos.data(L).astype(complex)
    assert np.allclose(Ldata, qL.full())

@pytest.mark.parametrize("method", ["sympy", "numpy","sparse"])
def test_liouvillian_qutip_consistency(method):
    test = np.array([[1,2],[3,4]])
    qrho = qu.Qobj(test)
    qL0 = qu.liouvillian(None,[qrho])
    qL1 = qu.liouvillian(qrho,[])
    qL2 = qu.liouvillian(qrho, [qrho])
    
    rho = sos.tidyup(test, method=method)
    L0 = sos.liouvillian(None,[rho])
    L1 = sos.liouvillian(rho,[])
    L2 = sos.liouvillian(rho, [rho])

    L0data = sos.data(L0).astype(complex)
    L1data = sos.data(L1).astype(complex)
    L2data = sos.data(L2).astype(complex)
    assert np.allclose(L0data, qL0.full())
    assert np.allclose(L1data, qL1.full())
    assert np.allclose(L2data, qL2.full())
