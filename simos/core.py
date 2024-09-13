from importlib.resources import open_binary
import numpy as _np
from . import backends
from typing import Union
from .constants import *
from .trivial import *
import functools
import operator
import warnings
import importlib.util

class WallClock:
    def __init__(self,start_time=0.0):
        self.time = start_time
    
    def inc(self, time):
        self.time += time

    def reset(self):
        self.time = 0.0
    
    def phase(self,freq,is_angular=False):
        if is_angular:
            return 2*_np.pi*freq*self.time
        else:
            return freq*self.time

######
# Global WallClock
######
globalclock = WallClock()


############################################################
# System Class
############################################################

class System(): 
    """ A class for representing composite quantum systems of spins and electronic levels. Stores all operators of the system.

    Attributes
    ----------
    
    method : str
        Backend that was used for system construction.
    system : list
        A list of dictionaries for all members of the system.
    subsystems : list
        A list of all separable subsystems of the Hilbert space, i.e. subsystems that were combined with a tensor product. Each subsystem is stored as a tuple holding the members' names.
    dim : int
        Hilbert space size of the system.
    dims : list
        Substructure of the Hilbert space (i.e. list of Hilbert space sizes of subsystems that were combined with a tensor product, in analogy to qutips dims attribute)
    id
        Identity operator of the system.
    
    """
    # Constructor 
    def __init__(self, system_array : Union[list, tuple], method = 'qutip'):
        """ Constructs an instance of the system class.

        :param list, tuple system_array: A series of dictionaries combined via potentially nested lists and tuples specifying all members of the system (spins or electronic levels) and how to construct their combined Hilbert space.
        
        :Keyword Arguments:
        * *method* (''str'') -- Specifies the backend.
        """
        # Ensure validity of the names.
        # This is a minimal check. Obvious issues (e.g. chosing same name for multiple levels) do not raise warnings or errors. 
        names = [i["name"] for i in totalflatten(system_array)]
        for name in names:
            if name.isalnum() == False:
                raise ValueError("Member names must only contain alphanumerical letters.")
        # Build the system.
        self.method = method 
        ops, dim, independent_names = _build(system_array, method = method) 
        self.subsystems = independent_names
        self.dim = dim 
        id_fun = getattr(getattr(backends,method), 'identity')
        self.id = id_fun(self.dim, dims = list(ops.items())[0][1].dims)
        for name, op in ops.items():
            setattr(self,name,op)
        self.system = totalflatten(system_array)
        self.dims = self.id.dims

    # Get, set or delete level properties after construction
    def get_property(self, level_name : str, property_name : str):
        """ Returns a property of a spin or an electronic level of the system.

        :param str level_name: Name of the spin or electronic level.
        :param str property_name: Name of the property.
        :returns: The property.
        """
        for level in self.system:
            if level["name"] == level_name:
                return level[property_name]
            
    def set_property(self, level_name :str , property_name : str , value):
        """ Sets a property of a spin or an electronic level of the system to a desired value.

        :param str level_name: Name of the spin or electronic level.
        :param str property_name: Name of the property.
        :param value: New value of the property.
        """
        if property_name == "val" or property_name == "name":
            raise ValueError("Properties name and val cannot be changed after construction.")
        else:
            index = self.get_index(self, level_name)
            self.system[index][property_name] = value 

    def del_property(self, level_name : str, property_name : str):
        """Deletes a property of a spin or an electronic level of the system.

        :param str level_name: Name of the spin or electronic level.
        :param str property_name: Name of the property.
        """    
        if property_name == "val" or property_name == "name":
            raise ValueError("Properties name and val cannot be changed after construction.")
        else:
            index = self.get_index(self, level_name)
            del self.system[index][property_name]

    # Get level positions in system list or Hilbert space substructure
    def get_index(self, level_name : str): 
        """Returns the position of a spin or an electronic level in the system-list.

        :param str level_name: Name of the spin or electronic level.
        :returns: Position index.
        """
        for ind_level, level in enumerate(self.system):
            if level["name"] == level_name:
                return ind_level
        return None     
    
    def get_subsystem_index(self, subsystem_name : str):
        """Returns the position of the subsystem that contains the spin or electronic level in the complete Hilbert space of the system.

        :param str level_name: Name of the spin or electronic level.
        :returns: Position index.
        """
        if subsystem_name in self.subsystems:
            return self.subsystems.index(subsystem_name)
        else:
            for ss in self.subsystems:
                if isinstance(ss, tuple):
                    if subsystem_name in ss:
                        return self.subsystems.index(ss)
                    
    # Basis Transformations
    def add_ghostspin(self, newname : str, spinnames : list, returnnames = False):
        """Couples N spins of the systems with spin values j1, j2 ... to N spins with values |j1-j2-...|, |j1-j2-...+1|, ..., |j1+j2+...| and stores the spin operators of the new, 
        coupled spins as attributes of the system class. 

        :param str newname: Will be used as a prefix for the spin operators of the new spins.
        :param list spinnames: List of the names of all spins that will be coupled.

        :Keyword Arguments:
        * *returnnames* (``bool``) -- If True, a list with the names of all new spins is returned. 
        """
        origname = newname
        newname =   newname + "_"
        # Extract projectors of the selected spins and verify that they are not ghostspins themselves.
        ps = []
        for spinname in spinnames:
            ps.append(getattr(self, spinname+"p"))
            tm = [_np.max(ps[-1][i].diag()) for i in ps[-1].keys()]
            if _np.min(tm) <  (1-1e-9): 
                raise ValueError("Ghostspin cannot be built from another ghostspin. Use native spins of system instead.") # no ghostspins allowed
        # Perform the coupling of the first two spins.
        newspin, U1dummy, U1 = _couple_spins(ps[0], ps[1], method = self.method, T = None)
        spin1mult = []
        spin1name = []
        spin1 = []
        for sm in newspin.keys():
            spin1mult.append(sm)
            spin1name.append(str(int(2*sm+1)))
            spin1.append(newspin[sm])
        # Couple the remaining spins.
        for i_spin in range(2, len(spinnames)):
            spin2 = ps[i_spin]
            tmp_spin1mult = []
            tmp_spin1name = []
            tmp_spin1 = []
            for i_sm, sm in enumerate(spin1mult):
                newspin, Umatrixdummy, Umatrix = _couple_spins(spin1[i_sm]["p"], spin2, method = self.method, T = U1)  
                for tmp_sm in newspin.keys():
                    tmp_spin1mult.append(tmp_sm)
                    tmp_spin1name.append(spin1name[i_sm]+"_"+str(int(2*tmp_sm+1)))
                    tmp_spin1.append(newspin[tmp_sm])
                if i_sm == 0:
                    tmp_U1 = Umatrix
                else:
                    tmp_U1 = _np.matmul(Umatrix, tmp_U1)
            spin1mult = tmp_spin1mult
            spin1name = tmp_spin1name
            spin1 = tmp_spin1
            U1 = _np.matmul(tmp_U1, U1)
        # Set the attributes in the systems.
        # Spin operators:
        for i_sm, sm in enumerate(spin1mult):
            newspin = spin1[i_sm]
            for opkey in newspin.keys(): # loop over operator keys
                setattr(self, newname+spin1name[i_sm]+opkey, newspin[opkey]) # set projector 
        # Transformation matrices:
        setattr(self, "to"+origname, U1)
        setattr(self, "from"+origname, _np.transpose(U1))
        # If required, return the names of the new spins
        if returnnames:
            return [newname+i for i in spin1name]

    def add_basis(self, T : _np.ndarray , basis_name : str, state_names : Union[list, tuple, _np.ndarray]):
        """Constructs a complete alternate basis for the system using a user-defined transformation matrix. The identity operators of the new basis states and 
        transformation matrices for conversion between the original and the new basis are stored as attributes of the system. 

        :param numpy.ndarray T: Transformation matrix, describing the transformation from the native system basis (Zeeman basis) to the new basis. 
        :param str basis_name: Name of the basis. 
        :param list, tuple, numpy.ndarray state_names: Series of names for the states of the new basis. A name must be provided for each individual state. 
        
        """
        # Verify validity of the input
        if len(state_names) != self.dim:
            raise ValueError("Number of basis names does not match number of states.")
        elif _np.shape(T)[0] !=_np.shape(T)[1]:
            raise ValueError("Transformation matrix must be diagonal.")
        elif _np.shape(T)[0] != self.dim:
            raise ValueError("Invalid dimension of transformation matrix. Must match the number of basis states.")
        else:
            # Construct and set the id operators of the new basis states in the original system basis. 
            tidyup_fun = getattr(getattr(backends, self.method), 'tidyup')
            ket2dm_fun = getattr(getattr(backends, self.method), 'ket2dm')
            for idx, name in enumerate(state_names):
                    levelket = tidyup_fun(T[idx, :])
                    leveldm = ket2dm_fun(levelket)
                    setattr(self, basis_name + "_" + name +"id", leveldm)
            # Also store the transformation matrices to and from the new basis. 
            setattr(self, "to"+basis_name, T)
            setattr(self, "from"+basis_name, _np.transpose(T))


# System variant for numba
if importlib.util.find_spec('numba') != None:
    import numba
    def create_system(system_array, method = 'qutip'):
        if method != 'numba':
            return System(system_array, method = method)
        else:
            obj = System(system_array, method = 'numpy')
            all_vars = vars(obj)
            class SystemNumba():
                def __init__(self):
                    self.method = 'numba'
            type_dict = {'method':numba.types.unicode_type}
            for key, value in all_vars.items():
                if isinstance(value, _np.ndarray):
                    if len(value.shape) == 1:
                        type_dict[key] = numba.complex128[:]
                    elif len(value.shape) == 2:
                        type_dict[key] = numba.complex128[:,:]
                    elif len(value.shape) == 3:
                        type_dict[key] = numba.complex128[:,:,:]
            
            nbc = numba.experimental.jitclass(type_dict)(SystemNumba)
            obj2 = nbc()
            for key, value in all_vars.items():
                if key in type_dict.keys():
                    setattr(obj2, key, value)
            return obj2



############################################################
# Subsystems  (system-aware partial trace)
############################################################

def subsystem(system, op, selection, keep = True):
    """ Extracts a subsystem from a system. If possible, the extraction is done with a partial trace. Otherwise,
    the subsystem is extracted with a projection operation. 
    
    :param System system: Instance of the System-class.
    :param op: Operator of the system.
    :param dict, list, str selection: Specifies the subsystem. Members can be specified with their names or their entire definition, mutiple members can be provided in a list. 

    :Keyword Arguments:
    * *keep* (``bool``) -- If True, selection specifies the system remaining after partial trace. If False, selection specifies the subsystem that is removed.
    
    :returns:
        - The operator in the subsytem.
        - :py:class:`dict` A dictionary with instructions on how to project operators of this subsystem back into the full system. 
    """

    # Process the "selection" input.
    if isinstance(selection, (str, dict)):
        selection = [selection]
    selection = totalflatten(selection)
    for idx, item in enumerate(selection):
        if isinstance(item, dict):
            selection[idx] = item["name"]
    # Prepare the dictionary that will hold the "reverse" instruction.
    reverse = {}
    # Perform a clean partial trace, if possible. 
    # Probe whether the subsystem is separable.
    ind_selection = []
    for ind_subsystem, subsystem in enumerate(system.subsystems):
        if all([i in selection for i in subsystem]):
            ind_selection.append(ind_subsystem)
    # If yes, perform a clean partial trace.
    if len(ind_selection) > 0: 
        # If keep is set to False, "invert" the selection.
        if keep == False:
            all_selection = list(_np.arange(0, len(system.subsystems)))
            ind_selection = [ elem for elem in all_selection if elem not in ind_selection]
        # Take the partial trace.
        op = op.ptrace(ind_selection)
        # Write the revert-recipe.
        reverse["tensor_ind"] = ind_selection
        reverse["tensor_dims"] = system.dims[0]
        return op, reverse 
    # If the subsystem is not separable, extract anyways.
    else:
        warnings.warn("The requested subsystem is not separable. The extraction might be an unphysical operation.")
        # If keep is set to False, "invert" the selection.
        if keep is False:
            selection = [i["name"] for i in system.system if i not in selection]
        tidyup_fun = getattr(getattr(backends, system.method), 'tidyup')
        # Find basis states in the current Hilbert space that belong to our selection.
        projector = 0*system.id
        for member in selection:
            projector += getattr(system, member+"id")
        pos =_np.where(_np.abs(projector.diag()) > 1e-9)[0]
        reverse["direct_pre"] = _np.min(pos)
        reverse["direct_post"] = (system.dim-1) - _np.max(pos)        
        # Determine the substructure of this Hilbert subspace.
        full_members = []
        part_members = []
        sel_members = []
        for member in system.system:                
            id = getattr(system, member["name"] + "id")
            idx = _np.where(_np.abs(id.diag()) > 1e-9)[0]
            if (_np.max(idx)>=_np.max(pos)) and  (_np.min(idx) <=_np.min(pos)):
                full_members.append(member)
                if member["name"] in selection:
                    sel_members.append(member)
            elif (_np.max(idx)>=_np.min(pos) and _np.max(idx)<=_np.max(pos)) or (_np.min(idx)>=_np.min(pos) and _np.min(idx)<=_np.max(pos)):
                part_members.append(member)
        # If possible, remove separable subsystems that were not selected. 
        if part_members == []:
            dim1d = [int(2*i["val"]+1) for i in full_members]
            dims = [dim1d, dim1d]
            out =  tidyup_fun(op[_np.min(pos):_np.max(pos)+1, _np.min(pos):_np.max(pos)+1], dims = dims)
            if len(sel_members) is not len(full_members):
                idx = [full_members.index(i) for i in sel_members]
                norm = _np.prod(dim1d)/_np.prod(_np.array(dim1d)[idx])
                out = (1/norm)*out.ptrace(idx)
                reverse["tensor_ind"] = idx
                reverse["tensor_dims"] = dims[0]
            return out, reverse 
        else:
            dims = [[len(pos)], [len(pos)]]
            out =  tidyup_fun(op[_np.min(pos):_np.max(pos)+1, _np.min(pos):_np.max(pos)+1], dims = dims)
        return out, reverse
    

def reverse_subsystem(system, op, reverse_instruction):
    """ Projects an operator formulated in the Hilbert space of a subsystem back into the Hilbert space of the full system.
    This function is intended as a counterpart to the subsystem method and difficult to use as a standalone method.
    
    :param System system: Instance of the System-class.
    :param op: Operator of the subsystem.
    :param dict reverse_instruction: A dictionary with instructions on how to project operators of this subsystem back into the full system. Such a dictionary is automatically generated by the subsystem method.

    :returns: The operator in the Hilbert space of the full system.  
    """    

    identity_fun = getattr(getattr(backends,system.method), 'identity')
    tensor_fun =  getattr(getattr(backends,system.method), 'tensor')
    # Revert partial trace. 
    if "tensor_ind" in reverse_instruction.keys():
        idx = reverse_instruction["tensor_ind"]
        dims =  reverse_instruction["tensor_dims"]
        if isinstance(idx, int):
            op = _project_up(op, dims, idx, method= system.method)
        else: 
            all_ops = _np.empty(len(dims), dtype = object)
            for ind_o in range(len(all_ops)):
                if ind_o in idx:
                    all_ops[ind_o] = op.ptrace(idx.index(ind_o))
                else:
                    all_ops[ind_o] = identity_fun(dims[ind_o])
            op = tensor_fun(all_ops)
    # Fill Hilbert.
    if "direct_pre" in reverse_instruction.keys():
        pre =  reverse_instruction["direct_pre"]
        post =  reverse_instruction["direct_post"]
        op = _fill_hilbert(op, pre, post, method = system.method)
    return op 


#####################################################################
# Functions for Construction of Operators for Spins and Spinsystems 
#####################################################################

def spinops(spindef, method="qutip", prefix=""):
    """Returns single-spin operators. 
     
    :param dict spindef: Definition of the spin for which operators should be build. Must contain a key 'val' that specifies the spin value. Can further contain a key 'type' to specify the spin type and a key 'name' for the spin name. If spin type is 'NV', 'NV-'or 'NV+' a series of additional operators are constructed.

    :Keyword Arguments:
    * *method* (``str``) -- Backend. Default is qutip.
    * *prefix* (``str``) -- A prefix that is added to the keys of the return dictionary.

    :returns: Dictionary that holds the single spin operators, keys are 'id','x','y','z','+','-','p' for 'default' spin type with spin > 0 and 'id' if spin = 0.
    """
    val = spindef['val']  # spin value
    N = int(2 * val + 1)  # Hilbert space size for individual spin
    if spindef.get('type') is None: 
        spindef['type'] = 'default'
    ops = {}
    # Pull backend-specific methods
    jmat_fun =  getattr(getattr(backends,method), 'jmat')
    diag_fun =  getattr(getattr(backends,method), 'diags')
    identity_fun = getattr(getattr(backends,method), 'identity')
    # Construct spin operators
    ops[prefix+'id'] = identity_fun(int(2 * val + 1))
    if val == 0:
        return ops
    ops[prefix+'x'] = jmat_fun(val, 'x')
    ops[prefix+'y'] = jmat_fun(val, 'y')
    ops[prefix+'z'] = jmat_fun(val, 'z')
    ops[prefix+'plus'] = jmat_fun(val, '+')
    ops[prefix+'minus'] = jmat_fun(val, '-')
    # Projectors
    projectors = {}
    for n in range(N):
        m_S = -val + n
        op_proj = _np.zeros(N, dtype=int)
        op_proj[n] = 1
        projectors[m_S] = diag_fun(op_proj,0)
    ops[prefix+'p'] = projectors
    # Special operators for specific spin types
    if spindef['type'] in ['NV', 'NV-', 'NV+']:
        if spindef['val'] != 1:
            raise ValueError('NVs must be spin=1 spins')
        ops[prefix+"op_x_red"] = 1 / 2 * (diag_fun([0, 1], 1) + diag_fun([0, 1], -1))
        ops[prefix+"op_y_red"] = 1 / 2j * (diag_fun([0, 1], 1) + diag_fun([0, -1], -1))
        ops[prefix+"op_z_red"] = diag_fun([0, 0, -1], 0)
        ops[prefix+"op_x_red2"] = 1 / 2 * (diag_fun([1, 0], 1) + diag_fun([1, 0], -1))
        ops[prefix+"op_y_red2"] = 1 / 2j * (diag_fun([1, 0], 1) + diag_fun([-1, 0], -1))
        ops[prefix+"op_z_red2"] = diag_fun([1, 0, 0], 0)
    return ops

def _project_up(operator, dimensions : list, idx : int, method='qutip'):
    """Projects an operator to a higher dimensional Hilbert space by performing a tensor product with a series of identity matrices.
        If the operator is a dictionary, list or tuple of operators, the function returns a dictionary, list or tuple of projected operators, respectively.
    
    :param operator: Operator that will be projected.
    :param list dimensions: Hilbert space dimension of the identity matrices used for the tensor product.
    :param int idx: Position of the operator in the series of tensor products.

    :Keyword Arguments:
    * *method* (``str``) -- Backend.

    :returns: Operator in a higher dimensional Hilbert space.
    """
    if isinstance(operator, dict):
        return {k:_project_up(operator[k],dimensions,idx, method=method) for k in operator.keys()}
    elif isinstance(operator, list):
        return [_project_up(op,dimensions,idx, method=method) for op in operator]
    elif isinstance(operator, tuple):
        return tuple([_project_up(op,dimensions,idx, method=method) for op in operator])
    else:
        identity_fun = getattr(getattr(backends,method), 'identity')
        tensor_fun = getattr(getattr(backends,method), 'tensor')
        identities = [identity_fun(e) for e in dimensions]   
        return tensor_fun(identities[:idx]+[operator]+identities[idx+1:])

def _fill_hilbert(operator, pre_fill:int, post_fill:int, method='qutip'):
    """ Projects an operator to a higher dimensional Hilbert space by performing a direct sum with empty operators.
        If the operator is a dictionary, list or tuple, the function returns a dictionary, list or tuple of projected operators, respectively.

    :param operator: 
    :param int pre_fill: Hilbert space dimension added before the operator.
    :param int post_fill: Hilbert space dimension added after the operator.

    :Keyword Arguments:
    * *method* (``str``) -- Backend.

    :returns: Operator in a higher dimensional Hilbert space.
    """
    if isinstance(operator, dict):
        return {k:_fill_hilbert(operator[k],pre_fill,post_fill, method=method) for k in operator.keys()}
    elif isinstance(operator, list):
        return [_fill_hilbert(op,pre_fill,post_fill, method=method) for op in operator]
    elif isinstance(operator, tuple):
        return tuple([_fill_hilbert(op,pre_fill,post_fill, method=method) for op in operator])
    else:
        directsum_fun = getattr(getattr(backends, method), 'directsum')
        return directsum_fun(operator, pre_fill, post_fill)

def _build(si : Union[dict, list, tuple], method='qutip'):
    """ Builds a system from an instruction in a recursive manner.
     
    :param dict, list, tuple si: Recipe for system construction. 

    :Keyword Arguments:
    * *method* (``str``) -- Backend.
    
    :returns: An instance of the System class.
    """
    if isinstance(si, dict):
        dimensions = int(si['val']*2+1)
        operators  = spinops(si, method=method, prefix=si['name'])
        independent_names = [si['name']]
    elif isinstance(si, (list, tuple)):
        # Subsystems need to be built
        operators = []
        dimensions = [] 
        si = flatten(si)
        independent_names = []
        for i,li in enumerate(si):
            curr_ops, curr_dims, names_to_add = _build(li, method=method)
            dimensions.append(curr_dims)
            operators.append(curr_ops)
            independent_names.append(names_to_add)
        # Sub_operators need to be combined with the tensor product
        if isinstance(si, list):
            for idx in range(len(operators)):
                operators[idx] = _project_up(operators[idx],dimensions,idx, method=method)
            dimensions = functools.reduce(operator.mul, dimensions)    
            # Make independent names a list
            independent_names = flatten(independent_names)
        # Sub_operators need to be combined with the direct sum       
        elif isinstance(si, tuple):
            for idx in range(len(operators)):
                pre_fill = int(_np.sum(dimensions[:idx]))
                post_fill = int(_np.sum(dimensions[idx+1:]))
                operators[idx] = _fill_hilbert(operators[idx],pre_fill,post_fill, method=method)
            dimensions = functools.reduce(operator.add, dimensions)
            # Make independent names a tuple
            independent_names = tuple(totalflatten(independent_names))
    else:
        raise ValueError("Unknown input type.")
    return fuse_dictionaries(operators), dimensions, independent_names


def _couple_spins(p1 : dict, p2 : dict, method = "qutip", T = None):
    # Load backend-specific functions.
    if method == 'sympy':
        cg_fun = backends.get_calcmethod('cg',  'symbolic')
    else:
        cg_fun = backends.get_calcmethod('cg',  'numeric')
    replace_fun = getattr(getattr(backends, method), 'changeelement')
    data_fun = getattr(getattr(backends, method ), 'data' )
    tidyup_fun = getattr(getattr(backends, method ), 'tidyup' )
    # Extract the spin values of spin 1 and spin 2.
    p1 = p1.copy()
    p2 = p2.copy()
    j1 = [i for i in p1.keys()][-1]
    j2 = [i for i in p2.keys()][-1]
    # Initialise the transformtion matrix as a np.ndarray.
    U = 0*data_fun(p1[j1])
    # Get the indices for all m_s - basis states of the selected spins. 
    # This also allows to catch an error if the spins are of different subsystems. 
    # Also store the "missing" positions, i.e. all other basis states. 
    pos = []
    for m1 in p1.keys():
        if isinstance(T, _np.ndarray):
            p1[m1] = p1[m1].transform(T)
        for m2 in p2.keys():
            if isinstance(T, _np.ndarray):
                p2[m2] = p2[m2].transform(T)
            try:
                pos.append(_np.where(abs((p1[m1]*p2[m2]).diag()) > 1e-9)[0])
            except Exception as error:
                print("Spins of different subsystems cannot be coupled.")
    pos = _np.array(pos)
    pos = pos.astype(int)
    mispos = _np.where( _np.in1d( _np.arange(U.shape[0]), pos.flatten()) == False)[0] 
    # Build the coupled spins.
    Jvals = _np.unique(_np.array([j1+j2, _np.abs(j1-j2)]))
    S = dict.fromkeys(Jvals)
    i = 0
    for J in Jvals:
        Jpos = []
        Mvals = _np.arange(-J, J+1, 1)
        # Exctract the spin operators for this multiplicity.
        operators  = spinops({"val": J}, method= method, prefix="") 
        S[J] = dict.fromkeys(operators.keys())
        S[J]["p"] = dict.fromkeys(Mvals)
        for M in Mvals:
            operator = 0*p1[j1].diag()
            # Build the coupled ket using Clebsch Gordan coefficients.
            for m1 in _np.arange(-j1, j1+1, 1):
                for m2 in _np.arange(-j2, j2+1, 1):
                    operator +=  cg_fun(j1, m1,j2, m2, J, M) * (p1[m1]*p2[m2]).diag()
            # Build density matrices from ket.
            # This has to be done separate for all duplicates of the states in the Hilbert space.
            S[J]["p"][M] =  0*p1[j1] 
            for j in range(pos.shape[1]):
                indices = pos[ : , j]
                subop = 0*p1[j1].diag()
                for indice in indices:
                    subop[indice] = operator[indice]
                    U[pos[i,j], indice] = operator[indice]
                subdm = tidyup_fun(_np.outer(_np.conj(subop), subop))
                subdm.dims = p1[j1].dims
                S[J]["p"][M] +=  subdm
            if isinstance(T, _np.ndarray):
                # Perform backtransformation to the original zeeman basis.
                S[J]["p"][M] = S[J]["p"][M].transform(_np.transpose(T)) 
            Jpos.append(i) # save positions of this J value
            i += 1
        # Also set the other spin operators.
        # Take care this is now eigenbasis of the coupled state so different backtransform will be performed.
        for opkey in operators.keys(): # loop over spin ops
            if opkey != "p":
                op = operators[opkey] # get operator 
                S[J][opkey] = 0*p1[j1] # initialise   
                for j in range(pos.shape[1]): # loop again over multiples of copies 
                    for ind_J1, J1 in enumerate(Jpos):
                        index1 = pos[J1, j]
                        for ind_J2, J2 in enumerate(Jpos): 
                            index2 = pos[J2, j]
                            S[J][opkey] = replace_fun(S[J][opkey], index1, index2, op[ind_J1, ind_J2])
    for i in mispos:
        U[i,i] = 1
    Upure = U
    if isinstance(T, _np.ndarray):
        U = _np.matmul(U, T) # first transform with T, then transform with U 
    for J in Jvals:
        for opkey in S[J].keys():
            if opkey != "p":
                S[J][opkey] = S[J][opkey].transform(_np.transpose(U))
    U[_np.abs(U) < 1e-10] = 0   # tidyup transformation matrix
    return S, U, Upure















