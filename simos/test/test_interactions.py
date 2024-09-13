import simos
import sympy as sp
import numpy as np



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