import numpy as _np
from . import backends
import warnings

###########################################################
# Automatic spin-Relaxation
###########################################################

def relaxation_operators(spin_system,temperature='infinity'):
    """Returns the collapse operators for a spin system for which T1 and T2 times have been supplied in the system dictionary.

    :param spin_system: The spin system for which the collapse operators are to be calculated.
    
    :Keyword Arguments:
        * *temperature* ('infinity') -- The temperature at which the system is to be simulated. The default is 'infinity'.
    
    :returns list: A list of collapse operators.
    """
    # doi: 10.1002/mrc.1242
    # doi:10.1016/j.chemphys.2011.03.014 (old)
    if temperature != 'infinity':
        raise NotImplementedError("Temperature dependence not implemented yet. Only infinite temperature limit is supported.")
    if temperature == 'infinity':
        d = 0
    all_cops = []
    if spin_system.method == "sympy":
        kind = "symbolic"
    else:
        kind = "numeric"
    sqrt = backends.get_calcmethod("sqrt", kind)
    for spin in spin_system.system:
        spin_name = spin['name']
        T1_defined = 'T1' in spin.keys()
        T2_defined = 'T2' in spin.keys()
        if T1_defined:
            T1 = spin['T1']
        if T2_defined:
            T2 = spin['T2']
 
        #projectors = getattr(spin_system, spin_name + 'p')
        if T2_defined:
            if T1_defined:
                rate = 2/T2 - 1/T1
            else:
                rate = 2/T2
            rate = sqrt(rate)
            #for level in projectors.values():
            #    #print(level)
            #    all_cops.append(rate*level)
            zop = getattr(spin_system, spin_name + 'z')
            all_cops.append(rate*zop)
        
        # Problematic for S>1 => Fix later
        if T1_defined:
            raiseing = getattr(spin_system, spin_name + 'plus')
            lowering = getattr(spin_system, spin_name + 'minus')
            rate_up = (1+d)/(4*T1/2)
            rate_down = (1-d)/(4*T1/2)
            all_cops.append(sqrt(rate_up)*raiseing)
            all_cops.append(sqrt(rate_down)*lowering)
            #if spin['val'] == 1/2:
            #    rate = sqrt(2/(spin['T1']))
            #elif spin['val'] == 1:
            #    #print("HEY")
            #    rate = sqrt(2/(3*spin['T1']))
            #else:
            #    raise(ValueError('Only Spins with S<1 supported.'))
            #all_level_keys = projectors.keys()
            #print(rate)
            #Sx = getattr(spin_system, spin_name + 'x')
            #for lev in all_level_keys:
            #    PP = projectors[lev]
            #    all_cops.append(rate*Sx*PP)
    return all_cops


###########################################################
# Collapse Operators for Level Transitions (laser excitation,
#, decay processes, chemical reaction) and all related
# functions
###########################################################


def tidyup_ratedict(spinsystem, rates):
    """Tidies up a rate dictionary, i.e. ensures that the keys have a correct format for usage with the transition_operators method.
    
    :param System spinsystem: The spin system for which the collapse operators are to be calculated. (future use)
    :param dict rates: A rate dictionary for level transitions. The keys of the dictionary describe the levels between which transitions occur, the values are the rates of the transitions. e.g. {'A->B': 1e6, 'B->A': 1e6} describes a transition from level A to level B with a rate of 1e6 s^-1 and vice versa.
    
    :returns dict: A tidy rate dictionary. 
    """
    tidyrates = {}
    # Method to add a rate to a dictionary.
    def _add_rate(name1, name2, rate, key):
        new_names = []
        for name_num in [name1, name2]: # further process source and sink descriptors
            sub_names = []
            full_names = []
            names = name_num.split(",") 
            for name in names:
                name = name.replace(" ", "") #get rid off white spaces 
                if '[' in name and  ']' in name: # if spin sublevel provided, extract 
                    value = float(name.split('[')[1].split(']')[0])
                    name = name.split('[')[0]            
                else:
                    value = None
                if name  in sub_names: #check if this level has already appeared once 
                    raise ValueError("Rate '" + key + "' has an invalid format. Level '" + name + "' may only occur once.")
                sub_names.append(name)
                if value is not None:
                    full_names.append(name+"["+str(value)+"]")   
                else:
                    full_names.append(name)
            full_names = sorted(full_names)
            full_names = "".join([i+"," for i in full_names[0:-1]] + [full_names[-1]])
            new_names.append(full_names)
        new_key = new_names[0]+"->" +new_names[1]
        if new_key in tidyrates.keys():
            raise ValueError("Rate '" + new_key + "' cannot be specified more than once.")
        else:
            tidyrates[new_key] = rate
    # Process all keys of the dictionary one by one. 
    for key, rate in rates.items(): 
            # Split the process in source and sink.
            parts_two = key.split("<->")
            if len(parts_two) == 2:
                _add_rate(parts_two[0], parts_two[1], rate, key)
                _add_rate(parts_two[1], parts_two[0], rate, key)
            else:
                part_forward = key.split("->")
                if len(part_forward) == 2:
                    _add_rate(part_forward[0], part_forward[1], rate, key)
                else:
                    part_backward = key.split("<-")
                    if len(part_backward) == 2:
                        _add_rate(part_backward[1], part_backward[0], rate, key)
                    else:
                        raise ValueError("Rate '" + key + "' has an invalid format.")  
    return tidyrates

def transition_operators(spinsystem, rates):
    """Returns the collapse operators for incoherent (optical) transitions described in a rate dictionary.

    :param System spinsystem: The spin system for which the collapse operators are to be calculated.
    :param dict rates: A rate dictionary for level transitions. The keys of the dictionary describe the levels between which transitions occur, the values are the rates of the transitions. e.g. {'A->B': 1e6, 'B->A': 1e6} describes a transition from level A to level B with a rate of 1e6 s^-1 and vice versa.

    :returns list: A list of collapse operators. 

    """
    # Extract all backend specific methods.
    if spinsystem.method == "sympy":
        sqrt_fun = backends.get_calcmethod("sqrt", "symbolic")
    else:
        sqrt_fun = backends.get_calcmethod("sqrt", "numeric")
    ket_fun = getattr(getattr(backends,spinsystem.method), 'ket')
    bra_fun = getattr(getattr(backends,spinsystem.method), 'bra')
    # Process all entries of a tidy rates dictionary:
    rates = tidyup_ratedict(spinsystem, rates)
    all_cops = [] 
    transform  = False
    for key, rate in rates.items():  
            # Split the rate description in source and sink.
            part_forward = key.split("->")
            name1 = part_forward[0]
            name2 = part_forward[1]
            # Further process the descriptions of source and sink levels.
            coll_states = []
            # Probe whether transition occurs in alternate basis.
            # Prohibit transition between states of different basis sets.
            if "_" in name1:
                basisname = name1.split("_")[0]
                if "_" not in name2 or name2.split("_")[0] != basisname:
                   warnings.warn("Source and sink of a transition are not specified in the same basis. Basis transformation will be performed for both, which may be erroneous.")
                #else:
                transform = True
                Tto = getattr(spinsystem, "to"+basisname)
                Tfrom = getattr(spinsystem,"from"+basisname) 
            for name_num in [name1, name2]:
                operator = spinsystem.id                  
                names = name_num.split(",") 
                for name in names:
                    if '[' in name: # if spin sublevel provided, extract 
                        value = float(name.split('[')[1].split(']')[0])
                        name = name.split('[')[0]
                        operator = operator * getattr(spinsystem, name+"p")[value]
                    else: # if no spin sublevel provided, take level-id 
                        operator = operator * getattr(spinsystem, name+"id")
                if transform:
                    operator = operator.transform(Tto)
                coll_states.append((_np.where(abs((operator).diag()) > 1e-9)[0]))
            # Construct the collapse operators.
            for source in coll_states[0]:
                bra = bra_fun(spinsystem.dim, int(source))
                for sink in coll_states[1]: 
                    ket = ket_fun(spinsystem.dim, int(sink)) 
                    c_ops = sqrt_fun(rate)*(ket*bra)
                    if transform:
                        c_ops = c_ops.transform(Tfrom)          
                    all_cops.append(c_ops)
    return all_cops 