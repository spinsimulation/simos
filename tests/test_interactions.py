import simos
import sympy as sp
import numpy as np
import pytest


#@pytest.mark.parametrize("method", ['qutip','numpy','sympy','sparse'])
class TestDipolar:
    #@pytest.mark.parametrize("val", np.arange(0.5,10,0.5))
    def test_consitency(self):
        g1 = sp.Symbol("g1")
        g2 = sp.Symbol("g2")
        x = sp.Symbol("x")
        y = sp.Symbol("y")
        z = sp.Symbol("z")

        D = simos.dipolar_spatial(g1, g2, [x, y, z], mode = "cart", case = "matrix")
        Dreal = simos.dipolar_spatial(simos.ye, simos.ye, [[1e-9, 0.5e-9, 0.2e-9]], mode = "cart", case = "matrix")

        Dsubs = D[1,1].subs({x:1e-9, y: 0.5e-9, z : 0.2e-9, g1 : simos.ye, g2 : simos.ye, sp.physics.quantum.constants.hbar : simos.hbar, sp.Symbol("mu0") : simos.mu_0, sp.pi: np.pi})
        assert Dreal[1,1] == Dsubs

        return None
    


class TestSpherical:
    def test_spherical1(self):
        testi = np.array([[1,2,-3], [40,5,6], [7,8,9]])
        spheri = simos.mat2spher(testi)
        testi2 = simos.spher2mat(spheri)
        def _all(test):
            return test.all()  
        assert _all((testi-testi2 < 1e-16))
        return None 


@pytest.mark.parametrize("method", ['qutip','numpy','sparse'])
class TestAnisotropicCouplingandInteractionHamiltonian:
    def test_spherical2(self,method):
        A = {'val': 1/2, 'name':'A'}
        B = {'val': 1/2 , 'name':'B'}
        # Construct the system. 
        s = simos.System([A,B], method = method)
        S = simos.spherbasis(s)
        testi = np.array([[1,2,-3], [40,5,6], [7,8,9]])
        test = simos.AnisotropicCoupling(mat = testi, frame = "lab", labonly = False)
        
        H1 = simos.interaction_hamiltonian(s, "A", test, spin2 = "B", frame = "lab")
        r0 = S["AB"][0][0][0]*test.spher("lab")[0][0] 
        r1 = S["AB"][1][1][0]*test.spher("lab")[1][0]   - S["AB"][1][1][-1]*test.spher("lab")[1][-1] - S["AB"][1][1][1]*test.spher("lab")[1][1]
        r2 = S["AB"][2][2][0]*test.spher("lab")[2][0]   - S["AB"][2][2][-1]*test.spher("lab")[2][-1] - S["AB"][2][2][1]*test.spher("lab")[2][1]  + S["AB"][2][2][-2]*test.spher("lab")[2][-2] + S["AB"][2][2][2]*test.spher("lab")[2][2] 
        H2 = r0+r1+r2


        assert (simos.data(H1) - simos.data(H2) < 1e-6).all()
        return None 
    
    def test_params(self,method):
        A = {'val': 1/2, 'name':'A'}
        B = {'val': 1/2 , 'name':'B'}
        # Construct the system. 
        s = simos.System([A,B], method = method)
        S = simos.spherbasis(s)
        testi = np.array([[1,2,3], [40,5,6], [7,8,9]])
        test = simos.AnisotropicCoupling(mat = testi, frame = "lab", labonly = False)
        params = test.params()

        test2 = simos.AnisotropicCoupling(iso = params["iso"], delta= params["delta"], eta = params["eta"], axz = params["axz"], ayz = params["ayz"], axy = params["axy"], euler = [params["alpha"], params["beta"], params["gamma"]])

            
        assert (np.abs(simos.data(test.mat("lab")) - simos.data(test2.mat("lab"))) < 1e-6).all()
        assert (simos.data(test.mat("pas")) - simos.data(test2.mat("pas")) < 1e-6).all()

        return None 