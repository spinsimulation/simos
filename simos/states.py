from . import backends
import scipy as _sc
import numpy as _np
from .propagation import rot
from .qmatrixmethods import tidyup, ket2dm
from .constants import kB

###########################################################
# Functions to generate intitial states
###########################################################


def state(system, statename : str):
    """ Generates an initial state (pure state).

    :param System system: An instance of the System class.
    :param str statename: A string specifying the initial state.

    :returns: The initial state (ket state vector).
    """
    operator = system.id
    sub_names = []
    names = statename.split(",") 
    for name in names:
        # Remove white spaces
        name = name.replace(" ", "")
        # Extract spin spin sublevel, if provided.  
        if '[' in name and  ']' in name: 
            value = float(name.split('[')[1].split(']')[0])
            name = name.split('[')[0]   
            operator = operator * getattr(system, name+"p")[value]         
        else:
            value = None
            operator = operator * getattr(system, name+"id")
        # Verify that members only occur once.
        if name  in sub_names:
            raise ValueError("State has an invalid format. Level '" + name + "' may only occur once.")
        sub_names.append(name)
    # Return ket.
    tidyup_fun = getattr(getattr(backends, system.method), "tidyup")
    out = tidyup_fun(operator.diag(), dims  = [system.dims[0], [1]])
    return out.unit()


def state_product(*states, method = "qutip"):
    """ Returns the tensor product of a list of quantum objects (states or density matrices).
     
    If the list contains at least one density matrix, the return object will be a density matrix (state vectors will be automatically converted).
    All states are normalized before tensor product is being executed.
    
    :param *states: List of quantum objects, supports state vectors and density matrices. 

    :Keyword Arguments:
        * *method* (''str'') -- Backend.

    :returns: Tensor product of all states.
    """
    isket_fun =  getattr(getattr(backends, method), 'isket')
    ket2dm_fun = getattr(getattr(backends, method), 'ket2dm')
    tensor_fun = getattr(getattr(backends, method), 'tensor')
    is_there_a_dm = False
    for state in states:
        if not isket_fun(state):
            is_there_a_dm = True

    if is_there_a_dm :
        dms = []
        for state in states:
            if isket_fun(state):
                dms.append(ket2dm_fun(state).unit())
            else:
                dms.append(state.unit())
        return tensor_fun(dms)
    else:
        states_normed = []
        for state in states:
            states_normed.append(state.unit())
    return tensor_fun(states_normed)


def pol_spin(polarization,method='qutip'):
    """Density matrix for a (partially) polarized spin 1/2.

    :params polarization: Polarization -1 .. 0 .. 1. Negative polarization values correspond to a polarization in the opposite spin state than positive polarization values.
    :returns: Density matrix for the (partially) polarized spin 1/2.
    """
    tidyup_fun =  getattr(getattr(backends,method), 'tidyup')
    return tidyup_fun([[(1+polarization),0],[0,(1-polarization)]])/2



def thermal_state(H, T = 300,  unit = "Kelvin"):
    """ Density matrix for a thermal state at a temperature T under a given  Hamiltonian H.

    :param H: System Hamiltonian.
    :param T: Temperature.

    :Keyword Arguments:
        * *unit* (''str'') -- Specifies temperature unit. Can be Kelvin or Degree Celcius.    

    :returns: Thermal density matrix. 
    """

    # Tidy up the input. 
    H = tidyup(H)
    # Handle temperature units.
    if unit == "Kelvin":
        T = T
    elif unit == "Celsius":
        T = T + 273.15
    else:
        raise ValueError("Invalid unit specified for temperature.")
    # Get the eigenvectos and values for the Hamiltonian.
    vals, vecs = H.eigenstates()
    Ezero = _np.min(vals)
    # Mutliply Boltzmann-weighted vectors.
    out = 0*H
    for val, vec  in zip(vals, vecs):
        bfac = _np.exp(val-Ezero)/(kB*T)
        out = out + bfac*ket2dm(vec)
    return out.unit()


