import numpy as _np
import math as _ma
import qutip as _qu
import sympy as _sp
from .constants import *
from .trivial import f2w, w2f, cart2spher, spher2cart
from . import backends
import copy
from itertools import combinations
from .constants import gyromagnetic_ratios

#############################################################
# Functions to simplify spin-spin and spin-field interactions
#############################################################

###########################################################
# Spin-Field Interactions
###########################################################

def field_interaction(spin_system, spin_name, y, *args, mode = "cart"):
    """ Returns the Hamiltonian for a generic spin-field ineraction.

    :param System spin_system: An instance of the System class.
    :param str spin_name: Name of the spin that is interacting with the field. 
    :param y: Coupling strength of the interaction, can be a scalar or a 3x3 matrix.
    :param *args: Magnetic/electric field, specified as a single argument [x,y,z] or three separate arguments x,y,z.

    :Keyword Arguments:
        * *mode* (``str``) --
          Can be 'cart' or 'spher', specifies whether the field is provided in cartesian or spherical coordinates.
    
    :returns: Hamiltonian of the spin-field interaction, a 3x3 matrix.
    """
    # Extract the field from *args.
    if len(args) == 3:
        field = _np.array([args[0], args[1], args[2]]).transpose()
    elif len(args) == 1:
        field = _np.array(args[0])
    else:
        raise ValueError("Wrong number of coordinates provided")
    if mode == "spher":
        field = spher2cart(field)
    # Construct the Hamiltonian of the interaction. 
    if len(_np.shape(y)) == 0 :
        field = y*_np.array(field)
    elif _np.shape(y) == (3,3): 
        field  = _np.matmul(_np.array(y), _np.array(field)) 
    else:
        raise ValueError("Gyromagnetic ratio has to be a scalar or a 3x3 matrix.")
    H = spin_system.id * 0
    H +=  field[0]*getattr(spin_system, spin_name + 'x')
    H +=  field[1]*getattr(spin_system, spin_name + 'y') 
    H +=  field[2]*getattr(spin_system, spin_name + 'z')
    return H

def auto_zeeman_interaction(spin_system, *args , mode = "cart", rotating_frame= []):
    """ Returns the combined Zeeman Hamiltonian for all spins in a system using the tabulated isotropic gyromagnetic
    ratios of their type (i.e. does not not consider any chemical shielding!). 

    :param System spin_system: An instance of the System class.
    :param *args: Magnetic field, specified as a single argument [x,y,z] or three separate arguments x,y,z. 

    :Keyword Arguments:
        * *mode* (``str``) --
          Can be 'cart' or 'spher', specifies whether the field is provided in cartesian or spherical coordinates.
        * *rotating_frame* (``list``) --
          A list of all spin types for which a rotating frame approximation is used (i.e. their Zeeman interactions are not considered). 
    
    :returns: The Zeeman Hamiltonian, a 3x3 matrix.

    """
    H = spin_system.id * 0
    for spin in spin_system.system:
        if spin["val"] >0: # only spins are processed, levels ignored 
            spin_name = spin['name']
            spin_type = spin["type"]
            if spin_type in rotating_frame: # check if rot_frame applies
                y = 0
            else:
                try:
                    y = getattr(gyromagnetic_ratios,'y' + spin_type)
                except:
                    raise ValueError("Unknown spin type. Type is " + str(spin_type) )
            H += field_interaction(spin_system, spin_name, y, *args, mode = mode)
    return H


###########################################################
# Spin-spin interaction with generic coupling tensor
###########################################################

def spinspin_coupling(spin_system, Name1, Name2, A,  approx = None): 
    """ Returns the Hamiltonian for a generic spin-spin interaction.

    :param System spin_system: An instance of the System class.
    :param str Name1: Name of the first spin of the pair. 
    :param str Name2: Name of the second spin of the pair. 
    :param A: Coupling strength of the interaction, can be a scalar or a 3x3 matrix.

    :Keyword Arguments:
        * *approx* (``str``) --
          Can be 'none', 'full', 'secular'. If 'None' or 'full', the interaction is not truncated. If 'secular' the interaction is truncated to the secular component.
    
    :returns: Hamiltonian of the spin-spin interaction, a 3x3 matrix.
    """
    # Enable case insensitive approximation.
    approx = approx.lower()
    # Fetch spin vectors.
    Sx = getattr(spin_system, Name1+'x')
    Sy = getattr(spin_system, Name1+'y')
    Sz = getattr(spin_system, Name1+'z')
    Ix = getattr(spin_system, Name2+'x')
    Iy = getattr(spin_system, Name2+'y')
    Iz = getattr(spin_system, Name2+'z') 
    # Test validity of A and assemble coupling Hamlitonian.
    A  = _np.array(A)
    if len(_np.shape(A)) == 0 :# y is a scalar
        if approx in  ['full', None]:
            return  A*(Sx*Ix + Sy*Iy + Sz*Iz)
        elif approx in ["secular"]:
            return  A*(Sz*Iz)
        else:
            raise ValueError('Invalid approximation provided.')
    elif _np.shape(A) == (3,3):  # y is a tensor
        H = 0*spin_system.id
        if approx in ["full", None]:
            H +=  A[0,0]*Sx*Ix +  A[0,1]*Sx*Iy + A[0,2]*Sx*Iz
            H +=  A[1,0]*Sy*Ix +  A[1,1]*Sy*Iy + A[1,2]*Sy*Iz
            H +=  A[2,0]*Sz*Ix +  A[2,1]*Sz*Iy + A[2,2]*Sz*Iz
        elif approx in ["secular"]:
            H +=  A[2,2]*Sz*Iz            
        else:
            raise ValueError('Invalid approximation provided.')            
    else:
        raise ValueError("A has to be a scalar or a 3x3 matrix.")


###########################################################
# Purely dipolar spin-spin interaction
###########################################################

def dipolar_spatial(y1, y2, *args, mode = 'spher', case = 'matrix'):
    """ Returns a dipolar coupling tensor in angular frequencies. 

        This is a vecorized routine, multiple dipolar couplings are calculated if a list of coordinates is provided.

    :param y1: Gyromagnetic ratio of the first spin, a scalar. 
    :param y2: Gyromagnetic ratio of the second spin, a scalar. 
    :param *args: Distance vector of the spin pair in the laboratory frame of reference. Can be specified as a single argument or three separate arguments. If multiple vectors are provided, a list of dipolar couplings will be returned for all of these.

    :Keyword Arguments:
        * *mode* (``str``) --
          Can be 'cart' or 'spher', specifies whether the vectors are provided in cartesian or spherical coordinates.
        * *case* (``str``) --
          Can be 'matrix' or 'alphabet'. If 'matrix', the function returns the coupling tensor as a 3x3 matrix. If 'alphabet', the dipolar alphabet is returned as a dictionary.

    :returns: The dipolar coupling tensor, either a 3x3 matrix or a dictionary with keys 'A','B','C','D','E','F'.
    """
    # Argument structure handling.
    if len(args) == 3:
        coords = _np.array([args[0], args[1], args[2]]).transpose()
    elif len(args) == 1:
        coords = _np.array(args[0])
    else:
        raise ValueError("Wrong number of coordinates provided")
    # Symbolic required?
    if all(backends.get_backend(x) != 'sympy' for x in [y1, y2, list(args)]):
        calcmode = "numeric"
    else:
        calcmode = "symbolic"
    sqrt_fun =  backends.get_calcmethod("sqrt", calcmode)
    sin_fun =  backends.get_calcmethod("sin", calcmode)
    cos_fun =  backends.get_calcmethod("cos", calcmode)
    pow_fun = backends.get_calcmethod("pow", calcmode)
    exp_fun = backends.get_calcmethod("exp", calcmode)
    array_fun = backends.get_calcmethod("array", calcmode)
    multelem_fun = backends.get_calcmethod("multiply", calcmode)
    hbar = backends.get_calcmethod("hbar", calcmode)
    mu0 = backends.get_calcmethod("mu0", calcmode)
    pi  = backends.get_calcmethod("pi", calcmode)
    I = backends.get_calcmethod("I", calcmode)
    D1 = array_fun(_np.identity(3, dtype = int))
    coords = array_fun(coords)
    pre = mu0*hbar*y1*y2/(4*pi)
    # If matrix desired, calculate dipolar coupling matrix D = pre/r^3 *(D1 - 3*D2)
    if case == 'matrix':
        if mode == 'spher':
            coords = spher2cart(coords)
        if len(coords.shape) == 1: 
            coords = array_fun([coords])
        pre =  pre*pow_fun(pow_fun(sqrt_fun(pow_fun(coords[:,0],2) + pow_fun(coords[:,1],2) + pow_fun(coords[:,2],2)), 3), -1)
        normalizer = pow_fun(sqrt_fun(pow_fun(coords[:,0],2) + pow_fun(coords[:,1],2) + pow_fun(coords[:,2], 2)),-1)
        coords[:,0] = multelem_fun(coords[:,0],normalizer)
        coords[:,1] = multelem_fun(coords[:,1],normalizer)
        coords[:,2] = multelem_fun(coords[:,2],normalizer)
        D2 = _np.einsum('ai, aj -> aij', coords, coords)
        Dipolar = (D1-3*D2)
        Dipolar =  _np.einsum('i, iab -> iab', pre, Dipolar)
        if _np.shape(Dipolar)[0] == 1: # de-pack 
            return array_fun(Dipolar[0])
        else:
            return array_fun(Dipolar)
    # If dipolar alphabet desired, get dipolar alphabet.
    elif case == 'alphabet': 
        if mode == 'cart':
            coords = cart2spher(coords)
        if len(coords.shape) == 1:
            r = array_fun([coords[0]])
            theta = array_fun([coords[1]])
            phi = array_fun([coords[2]])        
        else: 
            r = array_fun(coords[:,0])
            theta = array_fun(coords[:, 1])
            phi = array_fun(coords[:,2])
        Alphabet = {}
        pre = pre*pow_fun(pow_fun(r, 3), -1)
        Alphabet['a'] =     multelem_fun(pre, array_fun(_np.ones(len(r), dtype = int )-3*pow_fun(cos_fun(theta),2)))
        Alphabet['b'] = -1* multelem_fun(pre, array_fun(_np.ones(len(r), dtype = int )-3*pow_fun(cos_fun(theta),2)))/4
        Alphabet['c'] = -3* multelem_fun(pre, multelem_fun(sin_fun(theta), multelem_fun(cos_fun(theta), exp_fun(-I*phi))))/2
        Alphabet['d'] = -3* multelem_fun(pre, multelem_fun(sin_fun(theta), multelem_fun(cos_fun(theta), exp_fun(I*phi))))/2
        Alphabet['e'] = -3* multelem_fun(pre, multelem_fun(pow_fun(sin_fun(theta),2), exp_fun(-2*I*phi)))/4
        Alphabet['f'] = -3* multelem_fun(pre, multelem_fun(pow_fun(sin_fun(theta),2), exp_fun(2*I*phi)))/4
        if _np.shape(Alphabet["a"])[0] == 1: # de-pack 
            for key in Alphabet.keys():
                Alphabet[key] = Alphabet[key][0]
            return Alphabet
        else:
            return Alphabet
    else:
        raise ValueError('Invalid case provided. Case can only be "matrix or "alphabet".')
    

def dipolar_coupling(spinsystem, Name1, Name2, y1, y2, *args, mode = 'spher', approx = 'Full'):  
    """ Returns a dipolar coupling Hamiltonian.
   
    :param System spinsystem: An instance of the System class.
    :param str Name1: Name of the first spin.
    :param str Name2: Name of the second spin.
    :param y1: Gyromagnetic ratio of the first spin, a scalar. 
    :param y2: Gyromagnetic ratio of the second spin, a scalar. 
    :param *args: Distance vector of the coupled spins in the laboratory frame of reference. Can be specified as a single argument or three separate arguments. 

    :Keyword Arguments:
        * *mode* (``str``) --
          Can be 'cart' or 'spher', specifies whether the vectors are provided in cartesian or spherical coordinates.
        * *approx* (``str``) --
          Can be 'none', 'full','secular', 'pseudosecular', 'hfi'.  If 'None' or 'full', the interaction is not truncated. If 'secular' (pseudosecular) the interaction is truncated to the secular (pseudosecular) component. 

    :returns: The dipolar coupling Hamiltonian, a 3x3 matrix, in angular frequencies. 
     
    """   
    # Enable case insensitive approximation. 
    approx = approx.lower()
    if approx == 'pseudosecular':
        approx = 'ab'
    # Get necessary spin operators for all cases.
    Sz = getattr(spinsystem, Name1 + 'z')
    Iz = getattr(spinsystem, Name2 + 'z')
    Sx = getattr(spinsystem, Name1 + 'x')
    Ix = getattr(spinsystem, Name2 + 'x')
    Sy = getattr(spinsystem, Name1 + 'y')
    Iy = getattr(spinsystem, Name2 + 'y')
    # Calculate.
    if approx in ['full', 'secular','pseudosecular', 'hfi', 'none', None]:
        Dipolar = dipolar_spatial(y1, y2, *args, mode = mode, case = 'matrix')
        if approx in ['full', 'none', None]:
            Hdip =   Dipolar[0,0]*Sx*Ix +  Dipolar[0,1]*Sx*Iy + Dipolar[0,2]*Sx*Iz
            Hdip +=  Dipolar[1,0]*Sy*Ix +  Dipolar[1,1]*Sy*Iy + Dipolar[1,2]*Sy*Iz
            Hdip +=  Dipolar[2,0]*Sz*Ix +  Dipolar[2,1]*Sz*Iy + Dipolar[2,2]*Sz*Iz
        elif approx == 'secular':
            Hdip = Dipolar[2,2]*Sz*Iz
        elif approx == 'pseudosecular':
            Hdip =  Dipolar[0,0]*Sx*Ix + Dipolar[1,1]*Sy*Iy + Dipolar[2,2]*Sz*Iz
        elif approx == 'hfi':
            Hdip = Dipolar[2,0]*Sz*Ix +  Dipolar[2,1]*Sz*Iy + Dipolar[2,2]*Sz*Iz
    else:
        # additional operators 
        Splus = getattr(spinsystem, Name1 + 'plus')
        Iplus = getattr(spinsystem, Name2 + 'plus')
        Sminus = getattr(spinsystem, Name1 + 'minus')
        Iminus = getattr(spinsystem, Name2 + 'minus')
        Dipolar = dipolar_spatial(y1, y2, *args, mode = mode, case = 'alphabet')
        Hdip = 0*spinsystem.id
        for letter in approx:
            if letter == 'a':
                    Hdip += Dipolar[letter]*Sz*Iz
            elif letter == 'b':
                    Hdip += Dipolar[letter]*((Splus*Iminus + Sminus*Iplus))
            elif letter == 'c':
                    Hdip += Dipolar[letter]*((Splus*Iz + Sz*Iplus))
            elif letter == 'd':
                    Hdip += Dipolar[letter]*((Sminus*Iz + Sz*Iminus))
            elif letter == 'e':
                    Hdip += Dipolar[letter]*(Splus*Iplus)
            elif letter == 'f':
                    Hdip += Dipolar[letter]*(Sminus*Iminus)
            elif letter in [',', ';', ' ']:
                next
            else:
                raise ValueError('Invalid approximation provided.')
                 
    return Hdip

