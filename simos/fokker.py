from . import backends
from itertools import product
import numpy as _np
from .propagation import  evol
from .qmatrixmethods import tidyup
#from coherent import dipolar_coupling

######################################################
############# Basic Fokker Planck           ##########
######################################################

class _StochasticLiouvilleParameter:
    """Single Parameter for the stochastic Liouville equation.
    """
    def __init__(self, values=None, weights=None, dynamics={}):
        """Initializes the parameter.

        Using the dynamics dictionary, the parameter can be time-dependent. The keys of the dictionary are the orders of the time-derivatives and the values are the coefficients of the corresponding derivative.

        :param values: Values of the parameter. If None, the parameter is considered to be a constant.
        :type values: list, optional
        :param weights: Weights of the values. If None, the parameter is considered to be uniform.
        :type weights: list, optional
        :param dynamics: Dynamics of the parameter. If None, the parameter is considered to be static.
        :type dynamics: dict, optional
        """
        self.values = values
        self.weights = weights
        self.dynamics = dynamics
    
    def __str__(self):
        return f"Values: {self.values}, Weights: {self.weights}, Dynamics: {self.dynamics}"
    
    def __repr__(self):
        return f"Values: {self.values}, Weights: {self.weights}, Dynamics: {self.dynamics}"

class StochasticLiouvilleParameters(dict):
    """Stochastic Liouville Parameter set."""
    # Example definition:
    # params = simos.StochasticLiouvilleParameters()
    # params['a'].values = [0,1,2]
    # params['a'].dynamics = {1:1}
    # params["a"].weights = [1,1,1]

    def __init__(self,*args):
        """Initializes the parameter set.
        
        params = simos.StochasticLiouvilleParameters()
        params['a'].values = [0,1,2]
        params['a'].dynamics = {1:1}
        params["a"].weights = [1,1,1]
        """
        super().__init__()
        
    def __getitem__(self, key):
        if key not in self:
            super().__setitem__(key, _StochasticLiouvilleParameter())
        return super().__getitem__(key)
    
    def __setitem__(self, key, value):
        if key not in self:
            super().__setitem__(key, _StochasticLiouvilleParameter())
        super().__setitem__(key, value)

    def add(self, key, values, weights=None, dynamics={}):
        self[key].values = values
        self[key].weights = weights
        self[key].dynamics = dynamics

    # properties
    @property
    def separable(self):
        """Returns True if the parameter set is separable, i.e., the dynamics are separable. This is the case if all dynamics are None or only have a dynamics key set to 0.
        """
        # all dynamics are None or only have a key 0 and no other keys
        is_separable = True
        for key in self:
            if self[key].dynamics is not None:
                if len(self[key].dynamics.keys()) != 1:
                    is_separable = False
                    return is_separable
                elif 0 not in self[key].dynamics:
                    is_separable = False
                    return is_separable
        return is_separable
        
    def tensor_values(self):
        """returns an iterable that gives a dictionary of all possible values with their names
        
        Example:
        For elements self["a"].values = [1,2] and self["b"].values = [3,4]
        the output will be {"a":1, "b":3}, {"a":1, "b":4}, {"a":2, "b":3}, {"a":2, "b":4}
        """

        # get all keys
        keys = list(self.keys())
        # get all values
        values = [self[key].values for key in keys]

        # get all possible combinations
        for combination in product(*values):
            yield {keys[i]:combination[i] for i in range(len(keys))}
        
    @property
    def dof(self):
        """Number of degrees of freedom of the parameter set."""
        return _np.prod([len(self[key].values) for key in self])

    def tensor_mixer(self,method,boundary='periodic'):
        """Returns the mixer for the parameter set. 
        
        Determined by the dynamics key of the corresponding parameters.
        """
        differentiation_matrix = getattr(backends, method).differentiation_matrix
        tensor = getattr(backends, method).tensor
        tidyup = getattr(backends, method).tidyup
        individual_mixers = []
        ones = []
        total_size = 1
        mixer_dims = []
        for key in self:
            N = len(self[key].values)
            mixer_dims.append(N)
            total_size *= N
            M = _np.zeros((N,N)) 
            M = tidyup(M)
            if self[key].dynamics is not None:
                for order in self[key].dynamics:
                    M += differentiation_matrix(N,order,boundary=boundary)*self[key].dynamics[order]
            individual_mixers.append(M)
            one = _np.ones((N,N))
            ones.append(tidyup(one))
        # Now we build the mixer
        one_mixer = _np.zeros((total_size,total_size))
        #mixer_dims = [ones[0].dims[0]]*len(individual_mixers)
        #mixer_dims = [mixer_dims,mixer_dims]
        
        #mixer_dims = [mixer_dims,mixer_dims]
        one_mixer = tidyup(one_mixer)
        
        for i in range(len(individual_mixers)):
            curr = ones.copy()
            curr[i] = individual_mixers[i]
            ti = tensor(curr)
            ti.dims = [[ti.shape[0]],[ti.shape[1]]]
            one_mixer += ti
        return one_mixer
    
    def tensor_weights(self, renorm=True,symbolic=False):
        """Returns an iterable that gives the weights for all possible combinations of the parameter set.
        """
        if symbolic and renorm:
            raise NotImplementedError("Symbolic renormalization is not implemented yet.")
        weights = []
        for key in self:
            if self[key].weights is None:
                wi = _np.ones(len(self[key].values))
                if renorm:
                    wi = [wii/sum(wi) for wii in wi]
                weights.append(wi)
            else:
                wi = self[key].weights
                if renorm:
                    wi = [wii/sum(wi) for wii in wi]
                weights.append(wi)

        for combination in product(*weights):
            # multiply
            yield _np.prod(combination)

    
def stochastic_evol(Hfun, params: StochasticLiouvilleParameters, t, *rho, c_ops=[], e_ops=[], boundary_conditions='periodic', space='hilbert', wallclock='global',method=None):
    """Propagates a density matrix or a state vector under a stochastic Liouville equation.
    
    
    :param Hfun: Hamiltonian function.
    :param params: Parameters set spanning the Fokker-Planck space.
    :param t: Time at which the evolution is calculated.
    :param rho: Density matrix or state vector.
    :param c_ops: Collapse operators.
    :param e_ops: Expectation values.
    :param boundary_conditions: Boundary conditions for the differentiation matrix. Currently only 'periodic' is implemented.
    :param space: Space in which the calculation is performed. Either 'hilbert' or 'liouville'.
    :param wallclock: Wallclock time for the calculation to advance.
    :param method: Method to use for the calculation. If None, the method is chosen automatically based on rho.

    """
    if c_ops is None:
        c_ops = []
    if e_ops is None:
        e_ops = []
    if method is None:
        method = backends.get_backend(rho)
    tensor = getattr(backends, method).tensor
    identity = getattr(backends, method).identity
    liouvillian = getattr(backends, method).liouvillian
    if method == "sympy":
        I  = backends.get_calcmethod('I', 'symbolic')
    else:
        I  = backends.get_calcmethod('I', 'numeric')
    dm2ket = getattr(backends, method).dm2ket
    split = getattr(backends, method).split
    expect = getattr(backends, method).expect
    ket2dm = getattr(backends, method).ket2dm
    vec2dm = getattr(backends, method).vec2dm
    dm2vec = getattr(backends, method).dm2vec
    concatenate = getattr(backends, method).concatenate
    isket = getattr(backends, method).isket
    block_diagonal = getattr(backends, method).block_diagonal

    if len(rho) == 1:
        rho = rho[0]
        rho_is_pure = isket(rho)
        dims  = rho.dims 
    elif len(rho) == 0:
        rho = None
    else:
        raise ValueError("Too many arguments for rho.")

    
    if params.separable:
        # Create the Liouvillians for each parameter separately
        rhos = []
        Us = []
        for pi in params.tensor_values():
            Hi = Hfun(**pi)
            if rho is None:
                Ui = evol(Hi,t,c_ops=c_ops,wallclock=wallclock)
                Us.append(Ui)
            else:
                ri = evol(Hi,t,rho,c_ops=c_ops,wallclock=wallclock)
                rhos.append(ri)
        if rho is None:
            return block_diagonal(Us)
    else:
        if rho is not None:
            # Check if rho is pure
            try:
                rho = dm2ket(rho)
                rho_is_pure = True
            except ValueError:
                rho = dm2vec(rho)
            # We have to go to the liouville space in case there is relaxation or the state is mixed
            if(c_ops == [] and rho_is_pure):
                space = 'hilbert'
            else:
                space = 'liouville'

        # Create the full Liouvillian list
        L_list = []
        if space == 'hilbert':
            L_list = [Hfun(**pi) for pi in params.tensor_values()]
        elif space == 'liouville':
            L_list = [I*liouvillian(Hfun(**pi), c_ops) for pi in params.tensor_values()]
        else:
            raise ValueError("Space must be either 'hilbert' or 'liouville'.")
        

        # Now we build the Liouvillian / Hamiltonian
        # This is block-diagonal with the Liouvillians in L_list as blocks
        L_size = L_list[0].shape[0]
        L = block_diagonal(L_list)

        # Now we build the mixer
        mix = params.tensor_mixer(method,boundary=boundary_conditions)

        M = tensor([mix,identity(L_size)])
        D = L + I*M

        dof = params.dof

        if rho is None:
            U = evol(D,t,wallclock=wallclock)
            return U
        else:
            rhoF = concatenate([rho]*dof,dims = [[dof,rho.shape[0]],[1,1]])

            # Propagate
            rhoF = evol(D,t,rhoF)

            # Split rhoF
            rhos = split(rhoF, dof,is_hilbert=(space == 'hilbert'))
            # Convert to density matrices
            if rho_is_pure:
                rhos = [ket2dm(ri) for ri in rhos]
            #else:
            #    rhos = [vec2dm(ri) for ri in rhos]

    
    # Calculate expectation values
    e_values = []
    for e_op in e_ops:
        e_values.append([expect(e_op,ri) for ri in rhos])

    weights = params.tensor_weights()

    rhos = [wi*ri for wi,ri in zip(weights,rhos)]
    # Average rho (normalized by the weights already)
    if method == 'sympy':
        rhoAVG = rhos[0]
        for ri in rhos[1:]:
            rhoAVG += ri
        #rhoAVG = rhoAVG.unit()
    else:
        rhoAVG = sum(rhos).unit()
    
    rhoAVG = tidyup(rhoAVG, dims = dims)
    if e_ops == []:
        return rhoAVG
    else:
        return rhoAVG, e_values
