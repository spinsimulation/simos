import pytest
import numpy as np
import simos

def _all_close(test1,test2,method):
    if method == 'qutip':
        if hasattr(test1,'full'):
            test1 = test1.full()
        if hasattr(test2,'full'):
            test2 = test2.full()
        return np.allclose(test1,test2)
    elif method == 'numpy' or method == 'numba':
        return np.allclose(test1,test2)
    elif method == 'sympy':
        #return print((test1).evalf(),(test2).evalf())
        test1A = test1.evalf()
        test2A = test2.evalf()
        test1A = np.array(test1A).astype(np.complex128)
        test2A = np.array(test2A).astype(np.complex128)
        return np.allclose(test1A,test2A)

    elif method == 'sparse':
        if hasattr(test1,'toarray'):
            test1 = test1.toarray()
        if hasattr(test2,'toarray'):
            test2 = test2.toarray()
        return np.allclose(test1,test2)

@pytest.mark.parametrize("m", [1,2])
class TestDifferentiationMatrix:
    @pytest.mark.parametrize("N", np.arange(1,11))  
    def test_symbolic_consistency(self,m,N):
        num_matrix = simos.backends.numpy.differentiation_matrix(N, m)
        sym_matrix = simos.backends.sympy.differentiation_matrix(N, m).evalf()
        sym_matrix = np.array(sym_matrix).astype(np.float64)
        assert np.allclose(num_matrix, sym_matrix)

    @pytest.mark.parametrize("N", np.arange(1,11))
    def test_numeric_consistency(self,m,N):
        num_matrix = simos.backends.numpy.differentiation_matrix(N, m,method='optimal')
        num_matrix2 = simos.backends.numpy.differentiation_matrix(N, m,method='fft')
        assert np.allclose(num_matrix, num_matrix2)

    def test_manual(self,m):
        if m ==1:
            N = 4
            num_matrix = simos.backends.numpy.differentiation_matrix(N, m)
            manual = np.array([[0,1,0,-1],[-1,0,1,0],[0,-1,0,1],[1,0,-1,0]])/2
            assert np.allclose(num_matrix, manual)
        elif m == 2:
            N = 4
            num_matrix = simos.backends.numpy.differentiation_matrix(N, m)
            manual = np.array([[-3,2,-1,2],[2,-3,2,-1],[-1,2,-3,2],[2,-1,2,-3]])/2
            assert np.allclose(num_matrix, manual)
        else:
            raise ValueError("m must be 1 or 2")
        
        if m == 1:
            N = 3
            num_matrix = simos.backends.numpy.differentiation_matrix(N, m)
            manual = np.array([[0,1,-1],[-1,0,1],[1,-1,0]])*np.sqrt(3)/3
            assert np.allclose(num_matrix, manual)
        elif m == 2:
            N = 3
            num_matrix = simos.backends.numpy.differentiation_matrix(N, m)
            manual = np.array([[-2,1,1],[1,-2,1],[1,1,-2]])/3
            assert np.allclose(num_matrix, manual)
        else:
            raise ValueError("m must be 1 or 2")

@pytest.mark.parametrize("key, values, weights, dynamics", [
    ("a", [0, 1, 2], [1, 1, 1], {1: 1}),
    ("b", [3, 4], [0.5, 0.5], {0: 1}),
])
class TestStochasticLiouvilleParameters:
    def test_add_parameters(self, key, values, weights, dynamics):
        params = simos.StochasticLiouvilleParameters()
        params.add(key, values, weights, dynamics)
        assert key in params
        assert params[key].values == values
        assert params[key].weights == weights
        assert params[key].dynamics == dynamics

    def test_getitem_setitem(self, key, values, weights, dynamics):
        params = simos.StochasticLiouvilleParameters()
        param = simos.focker._StochasticLiouvilleParameter(values, weights, dynamics)
        params[key] = param
        assert params[key].values == values
        assert params[key].weights == weights
        assert params[key].dynamics == dynamics

    def test_separable_property(self, key, values, weights, dynamics):
        params = simos.StochasticLiouvilleParameters()
        params.add(key, values, weights, dynamics)
        assert params.separable == (len(dynamics) == 1 and 0 in dynamics)

    def test_dof_property(self, key, values, weights, dynamics):
        params = simos.StochasticLiouvilleParameters()
        params.add(key, values, weights, dynamics)
        assert params.dof == len(values)

    def test_tensor_values(self, key, values, weights, dynamics):
        params = simos.StochasticLiouvilleParameters()
        params.add(key, values, weights, dynamics)
        tensor_values = list(params.tensor_values())
        expected_values = [{key: v} for v in values]
        assert tensor_values == expected_values

    def test_tensor_weights(self, key, values, weights, dynamics):
        params = simos.StochasticLiouvilleParameters()
        params.add(key, values, weights, dynamics)
        tensor_weights = list(params.tensor_weights())
        expected_weights = [w / sum(weights) for w in weights]
        assert tensor_weights == expected_weights

@pytest.mark.parametrize("method", ["qutip","sparse","numpy","sympy"])
class TestStochachsticEvolution:
    def test_consistent_hopping_around_same_Hamiltonian(self,method):
        system_def = [{'name':'S','val':1/2}]
        s = simos.System(system_def,method=method)
        params = simos.StochasticLiouvilleParameters()
        params['a'].values = [0,1,2]
        params['a'].dynamics = {1:1.0}
        params["a"].weights = [1,1,1]

        def Hfun(a):
            if method == 'sympy':
                return (s.Sz+s.Sx).evalf()
            else:
                return s.Sz+s.Sx
        
        out = simos.stochastic_evol(Hfun,params,1,simos.pol_spin(1,method))
        out_reference = simos.evol(Hfun(0),1,simos.pol_spin(1,method))

        assert _all_close(out,out_reference,method)
    
    def test_consistent_powder_averaging(self,method):
        system_def = [{'name':'S','val':1/2}]
        s = simos.System(system_def,method=method)
        params = simos.StochasticLiouvilleParameters()
        params['a'].values = [0,1,2]
        params['a'].dynamics = {0:1}
        params["a"].weights = [1,1,1]

        def Hfun(a):
            return s.Sz+a*s.Sx
        
        out = simos.stochastic_evol(Hfun,params,1,simos.pol_spin(1,method))
        
        sot = []
        for i in range(3):
            Hi = Hfun(params['a'].values[i])
            sot.append(simos.evol(Hi,1,simos.pol_spin(1,method)))            

        if method == 'sympy':
            out_reference = (sot[0]+sot[1]+sot[2])/3
        else:
            out_reference = (sum(sot)/3).unit()

        assert _all_close(out,out_reference,method)
    
    def test_consistent_powder_averaging_with_mixed_state(self,method):
        system_def = [{'name':'S','val':1/2}]
        s = simos.System(system_def,method=method)
        params = simos.StochasticLiouvilleParameters()
        params['a'].values = [0,1,2]
        params['a'].dynamics = {1:0}
        params["a"].weights = [1,1,1]

        def Hfun(a):
            return s.Sz+a*s.Sx
        
        out,expectations = simos.stochastic_evol(Hfun,params,1,simos.pol_spin(0.7,method),e_ops = [s.Sz])
        
        sot = []
        for i in range(3):
            Hi = Hfun(params['a'].values[i])
            sot.append(simos.evol(Hi,1,simos.pol_spin(0.7,method)))            

        #out_reference = (sum(sot)/3).unit()
        #out_reference = (sot[0]+sot[1]+sot[2])/3
        if method == 'sympy':
            out_reference = (sot[0]+sot[1]+sot[2])/3
        else:
            out_reference = (sum(sot)/3).unit()

        assert _all_close(out,out_reference,method)
    
