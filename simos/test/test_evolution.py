import simos
import pytest
import qutip as qu
import sympy as sp
import numpy as np

@pytest.mark.parametrize("method", ['qutip','numpy','sympy','sparse'])
class TestEvolutionRates:
    #@pytest.mark.parametrize("val", np.arange(0.5,10,0.5))
    def test_nuclear_decay_with_fit(self,method):
        GS = {'val': 0 , 'name':'GS'}
        ES = {'val': 0 , 'name':'ES'}
        nuc = simos.System((GS, ES), method = method)

        rates = {"ES->GS": 2}
        rates = simos.tidyup_ratedict(nuc, rates)
        c_ops = simos.transition_operators(nuc, rates)
        dt = 1
        N = 100

        U = simos.evol(0*nuc.id, dt, c_ops = c_ops)
        obs  = []
        rho = nuc.ESid.copy()
        for idx_N in range(N):
            meas = simos.expect(nuc.ESid , rho)
            obs.append(np.real(meas))
            assert np.isclose(float(np.imag(meas)),float(0), rtol = 1e-6)
            rho = simos.applySuperoperator(U, rho)
        
        # curvefit
        from scipy.optimize import curve_fit
        def func(t, a, b):
            return a * np.exp(-b * t)
        
        popt, pcov = curve_fit(func, np.arange(N)*dt, obs, p0 = [2, 0.1])
        assert np.isclose(popt[1], rates["ES->GS"], rtol = 1e-6)
        return None
    
def test_analytically_correct():
    GS = {'val': 0 , 'name':'GS'}
    ES = {'val': 0 , 'name':'ES'}
    nuc = simos.System((GS, ES), method = 'sympy')

    r = sp.Symbol('r', real = True, positive = True)
    t = sp.Symbol('t', real = True, positive = True)
    rates = {"ES->GS": r}
    rates = simos.tidyup_ratedict(nuc, rates)
    c_ops = simos.transition_operators(nuc, rates)

    U = simos.evol(0*nuc.id, t, c_ops = c_ops)
    rho = simos.applySuperoperator(U, nuc.ESid)
    meas = simos.expect(nuc.ESid , rho)
    assert meas == sp.exp(-r*t)
    
    return None

########################
# Test a rotating frame Rabi
########################
class TestRotatingFrame:
    @pytest.mark.parametrize("method", ['qutip','numpy','sparse'])
    def test_rabi_numerical(self,method):
        system_def = []
        system_def.append({'name':'I','val':1/2})
        s = simos.System(system_def, method = method)

        Hx = s.Ix/10

        rho0 = simos.pol_spin(1,method=method)
        dt = 1
        N = 100

        rho = rho0.copy()
        store = []

        for i in range(N):    
            meas = simos.expect(rho0, rho)
            assert np.isclose(float(np.imag(meas)),float(0))
            store.append(np.real(meas))
            rho = simos.evol(Hx, dt, rho)

        # curvefit
        from scipy.optimize import curve_fit
        def func(x, f):
            return (np.cos(2*np.pi*f*x)+1)/2

        startval = 0.01
        popt, pcov = curve_fit(func, np.arange(N)*dt, store, p0 = [startval],method='lm')
        assert np.isclose(popt[0], 0.1/(2*np.pi), rtol = 1e-9)
        return None
    
    def test_rabi_analytical(self):
        system_def = []
        system_def.append({'name':'I','val':1/2})
        s = simos.System(system_def, method = 'sympy')

        f = sp.Symbol('f', real = True, positive = True)
        Hx = 2*sp.pi*f*s.Ix

        rho0 = simos.pol_spin(1,method='sympy')
        t = sp.Symbol('t', real = True, positive = True)
        rho = simos.evol(Hx, t, rho0)
        meas = simos.expect(rho0, rho)
        meas = sp.simplify(meas)
        solution_that_should_be = (sp.cos(2*sp.pi*f*t)+1)/2
        assert sp.simplify(meas - solution_that_should_be) == 0
    

    @pytest.mark.parametrize("method", ['qutip','numpy','sparse','sympy'])
    def test_rot_pi_pulse(self,method):
        system_def = []
        system_def.append({'name':'S','val':1/2})
        s = simos.System(system_def, method = method)

        rho0 = simos.pol_spin(1,method=method)

        def _all(test):
            if method == 'qutip':
                return test
            elif method == 'numpy' or method == 'numba':
                return test.all()
            elif method == 'sympy':
                return test
            elif method == 'sparse':
                return test.toarray().all()
    
        if method == 'sympy':
            pi = sp.pi
        else:
            pi = np.pi
        
        rho = simos.rot(s.Sx,pi,rho0)
        assert _all((rho == simos.pol_spin(-1,method=method)))

        rho = simos.rot(s.Sy,pi,rho0)
        assert _all((rho == simos.pol_spin(-1,method=method)))

        rho = simos.rot(s.Sy,pi,rho0)
        assert _all((rho == simos.pol_spin(-1,method=method)))

    def test_pihalf_pulse(self):
        # Test pi/2 pulse
        # A pi/2 pulse around the x-axis should produce |0> + |1>
        # A pi/2 pulse around the y-axis should produce |0> + i|1>
        #rho0 = simos.pol_spin(1,method=method)
        rho0 = qu.basis(2,0)
        system_def = []
        system_def.append({'name':'S','val':1/2})
        s = simos.System(system_def, method = "qutip")
        
        # y pulse
        wanted = np.zeros(2,dtype = complex)
        wanted[0] = 1/np.sqrt(2)
        wanted[1] = 1/np.sqrt(2)
        rho = simos.rot(s.Sy,np.pi/2,rho0)
        out = rho.full()
        out = out.flatten()
        print(out)
        print(wanted)
        assert np.allclose(out,wanted)

        # now a x pulse
        wanted[0] = 1/np.sqrt(2)
        wanted[1] = -1j/np.sqrt(2)
        rho = simos.rot(s.Sx,np.pi/2,rho0)
        out = rho.full()
        out = out.flatten()
        print(out)
        print(wanted)
        assert np.allclose(out,wanted)
    
    @pytest.mark.parametrize("method", ['qutip','numpy'])
    def test_consistency(self,method):
        system_def = []
        system_def.append({'name':'S', 'val':1/2, 'pol':1})
        s = simos.System(system_def,method=method)
        #rho0 = simos.gen_rho0(s)

        Uevol = simos.evol(s.Sz,1)*simos.evol(s.Sz+s.Sx,1)
        Uprop = simos.prop(s.Sz,1,carr1=np.array([1,0]),H1=s.Sx,magnus=False)
        def _all(test):
            if method == 'qutip':
                return test
            elif method == 'numpy' or method == 'numba':
                return test.all()
            elif method == 'sympy':
                return test
            elif method == 'sparse':
                return test.toarray().all()
        print(Uevol)
        print(Uprop)
        assert _all((Uevol == Uprop))


