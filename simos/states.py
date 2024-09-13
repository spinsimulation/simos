from . import backends
import scipy as _sc
from .propagation import rot


###########################################################
# Functions to generate intitial states
###########################################################


def state(system, statename : str, mode = "dm"):
    """ Generates an initial state (pure state).

    :param System system: An instance of the System class.
    :param str statename: A string specifying the initial state.

    :Keyword Arguments:
        * *mode* (''str'') -- Can be dm or ket. If dm, the state is returned as a density matrix, if ket a state vector is returned.

    :returns: The initial state.
    
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
    # Return ket or dm, depending on mode.
    if mode == "dm":
        return operator.unit()
    elif mode == "ket":
        tidyup_fun = getattr(getattr(backends, system.method), "tidyup")
        out = tidyup_fun(operator.diag(), dims  = [system.dims[0], [1]])
        return out.unit()
    else:
        raise ValueError("Invalid mode. Mode can be dm or ket.")


def rotate_state(system, dm, *args):
    """  Transforms a density matrix specified in a molecular frame of reference in to the laboratory frame. 
    The orientation of the molecular frame relative to the labframe can be specified in various formalisms.
    
    :param system: An instance of the System class.
    :param dm: The state (ket or dm) that will be rotated.
    :param *args: For qutip, numpy and sparse backends this can be a scipy.spatial.transform.Rotation object. For sympy this must, for all other backends it can, be the euler angles (zxy) of the rotation.  

    :returns: The density matrix after the coordinate transformation.
    
    """
    # If a ket is provided, make a dm.
    isket_fun = getattr(getattr(backends, system.method), 'isket')
    ket2dm_fun = getattr(getattr(backends, system.method), 'ket2dm')
    if isket_fun(dm):
        dm = ket2dm_fun(dm)
    # Connstruct the global x, y and z operators for spin rotations.
    xop = 0*system.id
    yop = 0*system.id
    zop = 0*system.id
    for level in system.system:
        if level["val"] > 0:
            xop += getattr(system, level["name"]+"x")
            yop += getattr(system, level["name"]+"y")
            zop += getattr(system, level["name"]+"z")
    # Determine the rotation matrix from the input.
    if len(args) ==1:
        R = args
        if isinstance(R, _sc.spatial.transform.Rotation): 
            R = args[0]
            euler = R.as_euler("zxy")
        else:
            euler = args[0]
    elif len(args) ==3:
        euler = [args[0], args[1], args[2]]
    else:
        raise ValueError("Invalid input for Rotation.")
    # Rotate the state.
    dm = rot(zop, euler[0], dm)
    dm = rot(xop, euler[1], dm)
    dm = rot(yop, euler[2], dm)
    return dm


def state_product(*states : list, method = "qutip"):
    """ Returns the tensor product of a list of quantum objects (states or density matrices).
     
    If the list contains at least one density matrix, the return object will be a density matrix (state vectors will be automatically converted).
    All states are normalized before tensor product is being executed.
    
    :param list *states: List of quantum objects, supports state vectors and density matrices. The routine will automatically noram
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
    """Density matrix for a paritally polarized spin 1/2.

    :params polarization: Polarization -1 .. 0 .. 1. Negative polarization values correspond to a polarization in the opposite spin state than positive polarization values.
    :returns: Density matrix for the (partially) polarized spin 1/2.
    """
    tidyup_fun =  getattr(getattr(backends,method), 'tidyup')
    return tidyup_fun([[(1+polarization),0],[0,(1-polarization)]])/2


def gen_rho0(spinsystem, NV_name='S'):
    rho_arr = []
    method = spinsystem.method
    for i in range(len(spinsystem.system)):
        spin = spinsystem.system[i]
        spin_type = spin['type']
        if spin_type in ['NV','NV-','NV+']:
            op0 = getattr(spinsystem, NV_name + 'p')[0]
            op0 = op0.ptrace(i)/2
            rho_arr.append(op0.unit())
        else:
            if 'pol' in spinsystem.system[i].keys():
                rho_arr.append(pol_spin(spinsystem.system[i]['pol']))
            else:
                 rho_arr.append(pol_spin(1).unit())
    rho0 = getattr(getattr(backends,method), 'tensor')(rho_arr).unit()
    return getattr(getattr(backends,method), 'tidyup')(rho0)
