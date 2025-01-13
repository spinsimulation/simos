import numpy as _np
import math as _ma
import qutip as _qu
import sympy as _sp
from .constants import *
from .trivial import cart2spher, spher2cart
from . import backends
import copy
from itertools import combinations
from .constants import gyromagnetic_ratios
from typing import Union
from .spherical import *
from .qmatrixmethods import eigenstates
import warnings

#######################################
# Spatial Part of Coherent Interactions
#######################################

class AnisotropicCoupling():
    """ A class for representing the spatial component of (anisotropic) coherent interactions.

    Attributes
    ----------
    
    R : 
       Scipy rotation object / sympy quaternion that specifies rotation from principal axes system to laboratory frame of reference.
    _mat_lab: 
        Matrix  representation in the laboratory frame of reference.
    _mat_pas: 
        Matrix representation in the principal axes frame of reference. 
    _spher_lab: dict
        Spherical tensor representation in the laboratory frame of reference.
    _spher_pas: dict
        Spherical tensor representation in the principal axes frame of reference.        
    
    """
    # Constructor. 
    def __init__(self, labonly = True, euler = [0,0,0], **kwargs):
        """ Constructs an instance of the AnisotropicCoupling class.

        :Keyword Arguments:
        * *labonly* (''bool'') -- If True, only the LAB representation is constructed.
        * *euler* -- Euler angles for PAS to LAB rotation, zyz convention.
        * *mat* -- Matrix representation of interaction.
        * *spher* -- Spherical tensor representation of interaction.
        * *iso* -- Isotropic value of interaction.
        * *delta* -- Anisotropy of interaction.
        * *eta* -- Asymmetry of interaction.
        * *axy* -- xy matrix element in PAS.
        * *axz* -- xz matrix element in PAS.
        * *ayz* -- yz matrix element in PAS.
        """
        # Input validity check.
        if (   ("mat" in kwargs and any(i in kwargs for i in ["spher", "iso", "delta", "eta", "axy", "axz", "ayz"]))
            or ("spher" in kwargs and any(i in kwargs for i in ["mat", "iso", "delta", "eta", "axy", "axz", "ayz"]))):
            raise ValueError("Input overdefines the interaction")
        if (("delta" in kwargs and not "eta" in kwargs) or ("eta" in kwargs and not "delta" in kwargs)) :
            raise ValueError("Definition of rank 2 is incomplete")
        if  _np.shape(euler) != (3,):
            raise ValueError("Wrong euler angle format.")
        # Probe for symbolic vs numeric. Load required math-methods. 
        if all(backends.get_backend(kwargs.get(x)) != 'sympy' for x in kwargs) and all(backends.get_backend(e) != 'sympy' for e in euler):
            self.calcmode = "numeric"
        else:
            self.calcmode = "symbolic"
        frommatrix_fun = backends.get_calcmethod("frommatrix", self.calcmode)
        fromeuler_fun = backends.get_calcmethod("fromeuler", self.calcmode)
        toeuler_fun = backends.get_calcmethod("toeuler", self.calcmode)
        rotate_fun  = backends.get_calcmethod("rotate", self.calcmode) 
        eigenstates_fun =  backends.get_calcmethod("eigenstates", self.calcmode) 
        pow_fun = backends.get_calcmethod("pow", self.calcmode) 
        sqrt_fun = backends.get_calcmethod("sqrt", self.calcmode)
        I = backends.get_calcmethod("I", self.calcmode)
        # Initialise interaction from matrix or spherical tensor.
        if ("mat" in kwargs or "spher" in kwargs):
            if "mat" in kwargs:
                mat_tmp = kwargs.get("mat")
                spher_tmp = mat2spher(mat_tmp)
            elif "spher" in kwargs:
                spher_tmp = kwargs.get("spher")
                mat_tmp = spher2mat(spher_tmp)
            if self.calcmode == "symbolic" and labonly == False:
                warnings.warn("Symbolic mode does not allow for automated PAS determination, assuming your input has been specified in PAS.")
                self._mat_pas = mat_tmp
                self._spher_pas = spher_tmp
                self.R = fromeuler_fun("zyz", euler)
                self._mat_lab   = rotate_fun(self._mat_pas, self.R,  False)
                self._spher_lab = mat2spher(self._mat_lab)
            elif labonly:
                self._mat_lab = mat_tmp
                self._spher_lab = spher_tmp
            else:
                # Probe for pas vs lab input. 
                # First prepare the spherical tensor and extract 
                # symmetric part of matrix (A2) for diagonalization. 
                sym_mat = spher2mat(spher_tmp, selrank = [2])
                _ , T_topas = eigenstates_fun(sym_mat)
                R_tmp =  frommatrix_fun(T_topas.transpose())  
                euler_tmp = toeuler_fun(R_tmp, "zyz")
                # If all euler angles are approx. 0 the interaction was provided
                # in the PAS already. Otherwise lab rep. was given.
                if all(_np.abs(euler_tmp[i]) < 1e-6 for i in range(3)):
                    self.R = fromeuler_fun("zyz", euler)
                    self._mat_pas = mat_tmp
                    self._spher_pas = spher_tmp
                    self._mat_lab   = rotate_fun(self._mat_pas, self.R,  False)
                    self._spher_lab = mat2spher(self._mat_lab)
                else:
                    self.R = R_tmp
                    self._mat_lab  = mat_tmp
                    self._spher_lab = spher_tmp
                    self._mat_pas = rotate_fun(self._mat_lab, self.R,  True)
                    self._spher_pas = mat2spher(self._mat_pas)  
        # Initialise interaction from parameters.
        else:
            self.R = fromeuler_fun("zyz", euler)
            self._spher_pas = {0 : { 0 : 0}, 1 : {-1 : 0, 0: 0, 1 :0}, 2 : {-2 :0, -1 : 0, 0: 0, 1:0, 2:0}}  
            # isotropic, rank 0, component
            if "iso" in kwargs:
                self._spher_pas[0][0] += -1* sqrt_fun(3)* kwargs.get("iso")
            # anisotropic, rank 2, component
            if "delta" in kwargs:
                self._spher_pas[2][0] +=  sqrt_fun(3*pow_fun(2, -1)) * kwargs.get("delta") 
                self._spher_pas[2][2] += -pow_fun(2, -1) * kwargs.get("delta") * kwargs.get("eta")
                self._spher_pas[2][-2]+= -pow_fun(2,-1)  * kwargs.get("delta") * kwargs.get("eta")    
            if "axy" in kwargs:
                self._spher_pas[1][0] +=  -I*sqrt_fun(2)*kwargs.get("axy")
            if "axz" in kwargs:
                self._spher_pas[1][1] += 1*kwargs.get("axz")
                self._spher_pas[1][-1] += 1*kwargs.get("axz")
            if "ayz" in kwargs:
                self._spher_pas[1][1] +=   -1*I*kwargs.get("ayz")
                self._spher_pas[1][-1] +=   +1*I*kwargs.get("ayz")
            self._mat_pas  = spher2mat(self._spher_pas)
            self._mat_lab = rotate_fun(self._mat_pas, self.R)
            self._spher_lab = mat2spher(self._mat_lab)
    # Getters.
    def mat(self, *args):
        if len(args) == 0:
            frame = "lab"
        else:
            frame = args[0]
        if frame == "lab":
            return self._mat_lab
        elif frame =="pas":
            return self._mat_pas
    def spher(self, *args):
        if len(args) == 0:
            frame = "lab"
        else:
            frame = args[0]
        if frame == "lab":
            return self._spher_lab
        elif frame =="pas":
            return self._spher_pas
    def params(self):
        params = {}
        trace_fun = backends.get_calcmethod("trace", self.calcmode)
        pow_fun  = backends.get_calcmethod("pow", self.calcmode)
        toeuler_fun = backends.get_calcmethod("toeuler", self.calcmode)
        params["iso"] =  pow_fun(3,-1)*trace_fun(self._mat_pas)
        params["delta"] = self._mat_pas[2,2] - params["iso"]
        params["eta"] =  (self._mat_pas[1,1] - self._mat_pas[0,0])*pow_fun(params["delta"],-1)
        params["axy"] = self._mat_pas[0,1]
        params["axz"] = self._mat_pas[0,2]
        params["ayz"] = self._mat_pas[1,2]
        euler = toeuler_fun(self.R, "zyz") 
        params["alpha"]= euler[0]
        params["beta"] = euler[1]
        params["gamma"] = euler[2]
        return params
    def alphabet(self, *args):
        if len(args) == 0:
            frame = "lab"
        else:
            frame = args[0]
        if frame == "lab":
            alph = {}
            alph["a"] = self._mat_lab[2,2]
            alph["b"] = 1/4*(self._mat_lab[0,0] + self._mat_lab[1,1])
            alph["c"] = 1/2* self._spher_lab[2][1]
            alph["d"] = -1/2*self._spher_lab[2][-1]
            alph["e"] = 1/2*self._spher_lab[2][2]
            alph["f"] = 1/2* self._spher_lab[2][-2]
        elif frame == "pas":
            alph = {}
            alph["a"] = self._mat_pas[2,2]
            alph["b"] = 1/4*(self._mat_pas[0,0] + self._mat_pas[1,1])
            alph["c"] = 1/2*self._spher_pas[2][1]
            alph["d"] = -1/2*self._spher_pas[2][-1]
            alph["e"] = 1/2*self._spher_pas[2][2]
            alph["f"] = 1/2*self._spher_pas[2][-2]
        return alph


def dipolar_spatial(y1, y2, *args, mode = 'spher', case = 'matrix'):
    """ Returns a dipolar coupling tensor in angular frequencies. This is a vecorized routine, multiple dipolar couplings are calculated if a list of coordinates is provided.

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
    # Symbolic required?
    if all(backends.get_backend(x) != 'sympy' for x in [y1, y2, list(args)]):
        calcmode = "numeric"
    else:
        calcmode = "symbolic"
    # Define Functions for symbolic/numeric case.
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
    # Argument structure handling.
    if len(args) == 3:
        coords = array_fun([args[0], args[1], args[2]]) 
    elif len(args) == 1:
        coords = array_fun(args[0])
    else:
        raise ValueError("Wrong number of coordinates provided")
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



########################################
##### Interaction Hamiltonian  #########
########################################


def interaction_hamiltonian(system, spin1, A, approx = "None", frame = "lab", mode = "cart", **kwargs):
    """ Hamiltonian of a coherent spin-spin or spin-field interaction. 

    :param System system: An instance of the System class.
    :param str spin1: Name of the (first) spin of the interaction.
    :param A: Spatial part of the interaction. Can be an instance of the AnisotropicCoupling class, a scalar or a 3x3 matrix. If a 3x3 matrix is provided,
              it is assumed that the latter specifies the interaction in the laboratory frame of reference.


    :Keyword Arguments:
        * *approx* (`` str`` ) --
          Specifies the approximation. Can be "None", "Secular", "Pseudosecular", "Hfi" or any letter "a"-"f" and combinations thereoff for spin-spin interactions and  "None" or "Secular" for spin-field interactions.
        * *frame* (``str``) --
          Can be 'pas' or 'lab', specifies whether the Hamiltonian is returned in the principal axis system of the interaction or the laboratory frame of reference.
        * *spin2* (`` str`` ) --     
          Specifies the name of the 2nd spin in spin-spin interactions.
        * *field* --     
          Specifies the magnetic field for spin-field interactions.
        * *mode* (``str``) --
          Can be 'cart' or 'spher', specifies whether the field is provided in cartesian or spherical coordinates.

    :returns: Hamiltonian of the spin-spin or spin-field interaction, a 3x3 matrix.
    """
    # Enable case insensitive approximation.
    approx = approx.lower()
    # Process spatial information.
    if not isinstance(A, AnisotropicCoupling) and len(_np.shape(A)) == 0:
        # Isotropic value was provided; generate accordingly.
        A = AnisotropicCoupling(iso = A, labonly = (frame == "lab"))
    elif not isinstance(A, AnisotropicCoupling) and _np.shape(A) == (3,3):
        # Matrix was provided,
        A = AnisotropicCoupling(mat = A, frame = "lab", labonly = (frame == "lab"))
    else:
        if not isinstance(A, AnisotropicCoupling):
            raise ValueError("Invalid format for spatial part of interaction.")
    # Fetch spin vectors of first spin.
    Sx = getattr(system, spin1+'x')
    Sy = getattr(system, spin1+'y')
    Sz = getattr(system, spin1+'z')
    # Prepare storage.
    H = 0*system.id
    # If spin-spin interaction. 
    if "spin2" in kwargs:
        spin2 = kwargs.get("spin2")
        Ix = getattr(system, spin2+'x')
        Iy = getattr(system, spin2+'y')
        Iz = getattr(system, spin2+'z')
        # Standard approximation case.
        if approx in ['full', 'secular','pseudosecular', 'hfi', 'none', None]:
            if isinstance(A, AnisotropicCoupling):
                A = A.mat(frame)
            if approx in ['full', 'none', None]:
                H =   A[0,0]*Sx*Ix +  A[0,1]*Sx*Iy + A[0,2]*Sx*Iz
                H +=  A[1,0]*Sy*Ix +  A[1,1]*Sy*Iy + A[1,2]*Sy*Iz
                H +=  A[2,0]*Sz*Ix +  A[2,1]*Sz*Iy + A[2,2]*Sz*Iz
            elif approx == 'secular':
                    H =  A[2,2]*Sz*Iz
            elif approx == 'pseudosecular':
                    H =  A[0,0]*Sx*Ix + A[1,1]*Sy*Iy + A[2,2]*Sz*Iz
            elif approx == 'hfi':
                    H =  A[2,0]*Sz*Ix +  A[2,1]*Sz*Iy + A[2,2]*Sz*Iz
        # Alphabet-like approximation cases.
        else:
            A  = A.alphabet(frame)
            Splus = getattr(system, spin1 + 'plus')
            Iplus = getattr(system, spin2 + 'plus')
            Sminus = getattr(system, spin1 + 'minus')
            Iminus = getattr(system, spin2 + 'minus')
            for letter in approx:
                if  letter == 'a':
                    H += A[letter]*Sz*Iz
                elif letter == 'b':
                    H += A[letter]*((Splus*Iminus + Sminus*Iplus))
                elif letter == 'c':
                    H += A[letter]*((Splus*Iz + Sz*Iplus))
                elif letter == 'd':
                    H += A[letter]*((Sminus*Iz + Sz*Iminus))
                elif letter == 'e':
                    H += A[letter]*(Splus*Iplus)
                elif letter == 'f':
                    H += A[letter]*(Sminus*Iminus)
                elif letter in [',', ';', ' ']:
                    next
                else:
                    raise ValueError('Invalid approximation provided.')    
    # If spin-field interaction. 
    else:
        # Symbolic vs numeric.
        if system.method == "sympy":
            calcmode = "symbolic"
        else:
            calcmode = "numeric"
        # Prep. the field input.
        array_fun = backends.get_calcmethod("array", calcmode)
        field = array_fun( kwargs.get("field"))
        if mode == "spher":
            field  = spher2cart(field)
        if _np.shape(field) == (3,1):
            field = field.transpose()
        # Calculate the Hamiltonian.
        A = A.mat(frame)
        matmul_fun = backends.get_calcmethod("matmul", calcmode)
        field  = matmul_fun(A,  field) 
        if  approx in ['full', 'none', None]:
            H += field[0]*getattr(system, spin1 + 'x')
            H +=  field[1]*getattr(system, spin1 + 'y') 
            H +=  field[2]*getattr(system, spin1 + 'z')
        elif approx in ["secular"]:
            H +=  field[2]*getattr(system, spin1 + 'z')     
        else:
            raise ValueError('Invalid approximation provided.')     
    return H 



def zeeman_interaction(system, spin, y, *args, mode = "cart", **kwargs):
    """ Zeeman Hamiltonian of a single spin.

    :param System system: An instance of the System class.
    :param str spin: Name of the selected spin.
    :param y: Spatial part of the interaction, can be a scalar , a 3x3 matrix, an instance of the AnisotropicCoupling class or a tuple/list with 3 (or 6) entries specifing isotropic chemical shift, span and skew (and 3 euler angles).
    :param *args: Magnetic field, specified as a single argument [x,y,z] or three separate arguments x,y,z. 

    :Keyword Arguments:
        * *mode* (``str``) --
          Can be 'cart' or 'spher', specifies whether the field is provided in cartesian or spherical coordinates.
        * *approx* (`` str`` ) --
          Specifies the approximation. Can be "None" or "Secular".

    :returns: Hamiltonian of the Zeeman interaction, a 3x3 matrix.
    """
    # Extract the field from *args.
    if len(args) == 3:
        field = [args[0], args[1], args[2]]
    elif len(args) == 1:
        field = args[0]
    else:
        raise ValueError("Wrong number of coordinates provided")
    if mode == "spher":
        field = spher2cart(field)
    # Construct the Zeeman Hamiltonian from Anisotropic coupling, scalar or matrix.
    if isinstance(y, AnisotropicCoupling) or len(_np.shape(y)) == 0 or _np.shape(y) ==(3,3):
        return interaction_hamiltonian(system, spin, y, field= field,  mode = mode, frame = "lab")
    # Construct the Zeeman Hamiltonian from isotropic shift, span, skew (and euler angles).
    else:
        if system.method == "sympy":
            calcmode = "symbolic"
        else:
            calcmode = "numeric"
        pow_fun = backends.get_calcmethod("pow", calcmode)
        omega = y[1]
        kappa = y[2]
        delta = pow_fun(6, -1) * (omega*kappa + 3*omega)
        eta = - pow_fun(2,-1)* (omega*kappa - omega) *pow_fun(delta, -1)
        if len(y) == 6:
            euler = [y[3], y[4], y[5]]
        else:
            euler = [0,0,0]
        A = AnisotropicCoupling(iso = y[0], delta = delta,  eta = eta, euler = euler, frame = "pas")
        return interaction_hamiltonian(system, spin, A, field= field, mode = mode, frame = "lab")


def auto_zeeman_interaction(spin_system, *args , mode = "cart", rotating_frame= []):
    """ Returns the combined Zeeman Hamiltonian for all spins in a system using the tabulated isotropic gyromagnetic
    ratios of their type (i.e. does not not consider any chemical shielding or anisotropy!). 

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
            H += zeeman_interaction(spin_system, spin_name, y, *args, mode = mode)
    return H


def dipolar_coupling(system, spin1, spin2, y1, y2, *args, mode = 'spher', **kwargs):  
    """ Returns a dipolar coupling Hamiltonian.
   
    :param System system: An instance of the System class.
    :param str spin1: Name of the first spin.
    :param str spin2: Name of the second spin.
    :param y1: Gyromagnetic ratio of the first spin, a scalar. 
    :param y2: Gyromagnetic ratio of the second spin, a scalar. 
    :param *args: Distance vector of the coupled spins in the laboratory frame of reference. Can be specified as a single argument or three separate arguments. 

    :Keyword Arguments:
        * *mode* (``str``) --
          Can be 'cart' or 'spher', specifies whether the vectors are provided in cartesian or spherical coordinates.
        * *approx* (``str``) --
          Can be 'none', 'full','secular', 'pseudosecular', 'hfi' or a string with all desired letters of the dipolar alphabet.  If 'None' or 'full', the interaction is not truncated. If 'secular' (pseudosecular) the interaction is truncated to the secular (pseudosecular) component. 
        
    :returns: The dipolar coupling Hamiltonian, a 3x3 matrix, in angular frequencies. 
     
    """   
    D = AnisotropicCoupling(mat = dipolar_spatial(y1, y2, *args, mode = mode, case = 'matrix'), labonly = True)
    return interaction_hamiltonian(system, spin1, D, spin2 = spin2,  frame = "lab", **kwargs)

def zfs_interaction(system, spin, ZFS, **kwargs):
    """ Returns the zero-field-splitting (ZFS) Hamiltonian.

    :param System system: An instance of the System class.
    :param str spin: Name of the spin for which the interaction is defined.
    :param ZFS: Spatial part of the ZFS interaction, can be a 3x3 matrix or an instance of the AnisotropicCoupling class or a tuple/list with 2 (or 5) entries which specify D and E (and 3 euler angles).

    :returns: The zero-field-splitting Hamiltonian, a 3x3 matrix. 

    """ 
    # define ZFS Hamiltonian from matrix or predefined Anistropic Coupling instance 
    if isinstance(ZFS, AnisotropicCoupling) or _np.shape(ZFS) ==(3,3):
        return interaction_hamiltonian(system, spin, ZFS, spin2 = spin, frame = "lab")
    # define ZFS Hamiltonian from D (and E) parameters, and, possibly, euler angles
    else:
        if system.method == "sympy":
            calcmode = "symbolic"
        else:
            calcmode = "numeric"
        array_fun = backends.get_calcmethod("array", calcmode)
        pow_fun = backends.get_calcmethod("pow", calcmode)
        ZFS_mat = array_fun(_np.zeros((3,3), dtype = int))
        if len(_np.shape(ZFS)) == 0:
            D = ZFS
            E = 0
            euler = [0,0, 0]
        else:
            D = ZFS[0] 
            if len(ZFS) > 1 :
                E = ZFS[1]
            else:
                E = 0
            if len(ZFS) == 5:
                euler = [ZFS[2], ZFS[3], ZFS[4]]
            else:
                euler = [0,0,0]
        ZFS_mat[0,0] = pow_fun(-3, -1)*D + E
        ZFS_mat[1,1] = pow_fun(-3, -1)*D - E
        ZFS_mat[2,2] = 2*pow_fun(3, -1)*D
        print(ZFS_mat)
        A = AnisotropicCoupling(mat = ZFS_mat, euler = euler, frame = "pas")
        return interaction_hamiltonian(system, spin, A, spin2 = spin, frame = "lab", **kwargs)
    


def quad_interaction(system, spin, QUAD, **kwargs):
    """ Returns the quadrupole coupling Hamiltonian.
   
    :param System system: An instance of the System class.
    :param str spin: Name of the spin for which the interaction is defined.
    :param QUAD: Spatial part of the quadrupole interaction, can be a 3x3 matrix or an instance of the AnisotropicCoupling class or a tuple/list with 2 (or 5) entries which specify anisotropy and assymetry of the interaction  (and 3 euler angles).
        
    :returns: The quadrupole Hamiltonian, a 3x3 matrix.
    """ 
    # define quadrupole Hamiltonian from matrix or predefined Anistropic Coupling instance 
    if isinstance(QUAD, AnisotropicCoupling) or _np.shape(QUAD) ==(3,3):
        return interaction_hamiltonian(system, spin, QUAD, spin2 = spin, frame = "lab")
    # define ZFS Hamiltonian from CQ and anisotropy parameters, and, possibly, euler angles
    else:
        delta = QUAD[0]
        eta = QUAD[1]
        if len(QUAD) == 5:
            euler = [QUAD[2], QUAD[3], QUAD[4]]
        else:
            euler = [0,0,0]
        A = AnisotropicCoupling(delta = delta, eta = eta, euler = euler)
        return interaction_hamiltonian(system, spin, A, spin2 = spin, frame = "lab", **kwargs)



