import simos as sos
import pytest
import numpy as np
import sympy as sp


@pytest.mark.parametrize("method", ['qutip','numpy','sympy','sparse'])
def test_construction_one_spin(method):
    systemarray = [{'name':'A','val':1/2}]
    s = sos.System(systemarray,method=method)
    assert hasattr(s,'Aid')
    assert hasattr(s,'Ax')
    assert hasattr(s,'Ay')
    assert hasattr(s,'Az')
    assert hasattr(s,'Aplus')
    assert hasattr(s,'Aminus')
    assert hasattr(s,'Ap')
    s.Ap[-0.5]
    s.Ap[0.5]


def check_all_spinops_there(sys,name,val):
    assert hasattr(sys,name+'id')
    assert hasattr(sys,name+'x')
    assert hasattr(sys,name+'y')
    assert hasattr(sys,name+'z')
    assert hasattr(sys,name+'plus')
    assert hasattr(sys,name+'minus')
    assert hasattr(sys,name+'p')
    for m in np.arange(-val,val+1):
        getattr(sys,name+'p')[m]


@pytest.mark.parametrize("method", ['qutip','numpy','sympy','sparse'])
def test_two_spins(method):
    # Construction
    S = {'val':1/2,'name':'S'}
    I = {'val':1/2,'name':'I'}
    s = sos.System([S,I], method=method)
    check_all_spinops_there(s,'S',1/2)
    check_all_spinops_there(s,'I',1/2)

    # Coupled product states
    s.add_ghostspin('C',['S','I'])
    check_all_spinops_there(s,'C_3',1)
    s.C_1p[0]
    assert hasattr(s, 'C_1id')
    # Make sure that the coupled product states are correctly defined
    def _all(test):
        if method == 'qutip':
            return test
        elif method == 'numpy' or method == 'numba':
            return test.all()
        elif method == 'sympy':
            return test
        elif method == 'sparse':
            return test.toarray().all()
        
    # Partial trace
    op = s.Sx + s.Ix
    op_red, op_red2, r = sos.subsystem(s, op, ['S'], keep = True)

    # Check that transformation matrices exist
    assert hasattr(s,'toC')
    assert hasattr(s,'fromC')


@pytest.mark.parametrize("method", ['qutip','numpy','sympy','sparse'])
def test_two_electronic_levels(method):
    # Construction
    GS = {'val':0, 'name':'GS'}
    ES = {'val':0, 'name':'ES'}
    a = sos.System((GS,ES), method=method)
    assert hasattr(a,'GSid')
    assert hasattr(a,'ESid') 

    T = np.array([[0.5,0.5],[0.5,-0.5]])
    a.add_basis(T,'A',['g','e'])
    assert hasattr(a,'A_gid')
    assert hasattr(a,'A_eid')
    op = a.A_gid + a.A_eid
    op_red, op_red2, r = sos.subsystem(a, op, ['ES'], keep = True)


@pytest.mark.parametrize("method", ['qutip','numpy','sympy','sparse'])
def test_state(method):
    S = {'val':1/2,'name':'S'}
    I = {'val':1/2,'name':'I'}
    s = sos.System([S,I], method=method)

    k0 = sos.state(s, 'S[0.5],I[-0.5]')
    r0 = sos.ket2dm(k0)
    rv0 = sos.dm2vec(r0)

    def _all(test):
        if method == 'qutip':
            return test
        elif method == 'numpy' or method == 'numba':
            return test.all()
        elif method == 'sympy':
            return test
        elif method == 'sparse':
            return test.toarray().all()
    
    # lets co back
    r0recon = sos.vec2dm(rv0)
    assert _all(r0recon == r0)
    k0recon = sos.dm2ket(r0)
    print(k0recon)
    print(k0)
    #todo fix qutip
    if method != 'qutip':
        assert _all(k0recon == k0)

@pytest.mark.parametrize("method", ['qutip','numpy','sympy','sparse'])
def test_rates(method):
    GS = {'val':0, 'name':'GS'}
    ES = {'val':0, 'name':'ES'}
    a = sos.System((GS,ES), method=method)
    rates = {}
    kex = 5e6
    kdec = 10e6
    rates["GS -> ES"] = kex
    rates["ES -> GS"] = kdec
    cops = sos.transition_operators(a, rates)

    S = {'val':1/2,'name':'S','T1':1e-3}
    s = sos.System([S], method=method)
    cops = sos.relaxation_operators(s)

    H0 = s.Sz
    k0 = sos.state(s, 'S[0.5]')
    r0 = sos.ket2dm(k0)
    r1 = sos.evol(H0,1,r0,c_ops=cops)

def test_paper_listing1():
    import sympy as sp
    import simos as sos
    r,th,phi,y1,y2,t = sp.symbols('r,theta,phi,gamma_1,gamma_2,t', real=True,positive=True)
    system_def = []
    system_def.append({'name':'S','val':1/2})
    system_def.append({'name':'I','val':1/2})
    s = sos.System(system_def,'sympy')
    H0 = sos.dipolar_coupling(s,'S','I',y1,y2,r,th,phi,approx='secular')
    Hsimple = sos.symbolic_replace_all(H0,1,r,th).simplify()
    psi0 = sos.state(s,'S[-0.5],I[-0.5]')
    psi = sos.rot(s.Sx,sp.pi/2,psi0)
    psi = sos.evol(Hsimple,t/2,psi)
    psi = sos.rot(s.Sy,sp.pi,psi)
    psi = sos.rot(s.Iy,sp.pi,psi)
    psi = sos.evol(Hsimple,t/2,psi)
    psi = sos.rot(s.Sx,-sp.pi/2,psi)
    m = sos.expect(s.Sz,psi)
    mexpected = sp.cos((3*t*sp.sin(th)**2)/(8*sp.pi*r**3) - t/(4*sp.pi*r**3))
    test = sp.simplify(sp.re(m).simplify()[0] - mexpected/2)
    assert test == 0


def test_paper_listing1_num():
    import sympy as sp
    import simos as sos
    #r,th,phi,y1,y2,t = sp.symbols('r,theta,phi,gamma_1,gamma_2,t', real=True,positive=True)
    r, th, phi = 1, 1, 1
    y1, y2 = 10e9,10e9
    system_def = []
    system_def.append({'name':'S','val':1/2})
    system_def.append({'name':'I','val':1/2})
    s = sos.System(system_def,'numpy')

    H0 = sos.dipolar_coupling(s,'S','I',y1,y2,r,th,phi,approx='secular')
    psi0 = sos.state(s,'S[-0.5],I[-0.5]')
    store = []
    dt = 0.01
    for i in range(100):
        psi = sos.rot(s.Sx,np.pi/2,psi0)
        psi = sos.evol(H0,i*dt/2,psi)
        psi = sos.rot(s.Sy,np.pi,psi)
        psi = sos.rot(s.Iy,np.pi,psi)
        psi = sos.evol(H0,i*dt/2,psi)
        psi = sos.rot(s.Sy,np.pi,psi)
        psi = sos.rot(s.Sx,-np.pi/2,psi0)
        m = sos.expect(s.Sz,psi)
        store.append(m)

