import numpy as _np
from . import backends
from .constants import eps_0

###########################################################
# TRIVIAL CONVERSION FUNCTIONS, sympy compatible
###########################################################

def f2w(freq):
	"""Takes a frequency and returns the corresponding angular frequency. 
	
	:param freq: A single frequency or a list, tuple or a numpy array of frequencies.
	:returns: The corresponding angular frequency.
	""" 
	# Symbolic required?
	if backends.get_backend(freq) == 'sympy': # Sympy 
		mode = "symbolic"
	else:
		mode = "numeric"    
	pi = backends.get_calcmethod("pi", mode)
	array_fun = backends.get_calcmethod("array", mode)
	try: # if multiple freqs provided 
		len(freq)
		return array_fun([2*pi*fi for fi in freq])
	except TypeError: # if single freq is provided 
		return 2*pi*freq
	
def w2f(freq):
	"""Takes an angular frequency and returns the corresponding frequency. 

	:param freq: A single angular frequency or a list, tuple or numpy array of angular frequencies.
	:returns: The corresponding frequency.
	"""
	# Symbolic required?
	if backends.get_backend(freq) == 'sympy': # Sympy 
		mode = "symbolic"
	else:
		mode = "numeric"    
	pi = backends.get_calcmethod("pi", mode)
	array_fun = backends.get_calcmethod("array", mode)
	try: # if multiple freqs provided 
		len(freq)
		return array_fun([fi/(2*pi) for fi in freq])
	except TypeError: # if single freq is provided 
		return freq/(2*pi)

def rad2deg(angle):
	""" Takes an angle in radians and returns it in degrees. 
	
	:param angle: A single angle or a list, tuple or numpy array of angles in radians.
	:returns: The angle in degrees. 
	"""
	# Symbolic required?
	if backends.get_backend(angle) == 'sympy': # Sympy 
		mode = "symbolic"
	else:
		mode = "numeric"    
	pi = backends.get_calcmethod("pi", mode)
	array_fun = backends.get_calcmethod("array", mode)
	# Transform
	try:
		len(angle)
		return array_fun([ai/pi*180 for ai in angle])
	except TypeError:
		return angle/pi*180

def deg2rad(angle):
	"""Takes an angle in degrees and returns it in radians. 
	
	:param angle: A single angle or a list, tuple or numpy array of angles in degrees.
	:returns: The angle in radians. 
	"""

	# Symbolic required?
	if backends.get_backend(angle) == 'sympy': # Sympy 
		mode = "symbolic"
	else:
		mode = "numeric"    
	pi = backends.get_calcmethod("pi", mode)
	array_fun = backends.get_calcmethod("array", mode)
	# Transform
	try:
		len(angle)
		return array_fun([ai/180*pi for ai in angle])
	except TypeError:
		return angle/180*pi

def cart2spher(*args,rad=True):
	""" Transforms cartesian into spherical coordinates.
		
	:param *args: A single or multiple sets of cartesian x, y, z  coordinates. Coordinates can be provided as a single argument, i.e. [x, y, z] or  [[x1,y1,z1], [x2,y2,z2]] or  as three separate arguments, i.e.  x, y, z or  [x1, x2], [y1, y2], [z1, z2] 
		
	:returns: The corresponding spherical coordinates.
		
	"""
	# Symbolic required?
	if all(backends.get_backend(a) != 'sympy' for a in args):
		mode = "numeric"
	else:
		mode = "symbolic"  
	sqrt_fun =   backends.get_calcmethod("sqrt", mode)
	pow_fun = backends.get_calcmethod("pow", mode)
	arctan2_fun = backends.get_calcmethod("arctan2", mode)
	array_fun = backends.get_calcmethod("array", mode)
	# Handle input structure
	if len(args) == 3: # if three arguments are provided, first is x, second y, third z
		x = array_fun(args[0])
		y = array_fun(args[1])
		z = array_fun(args[2])
	elif len(args) == 1: # if single argument is provided 
		arg = array_fun(args[0])
		if len(_np.shape(args[0])) == 1: # if only one set (x,y,z) ist given
			x = arg[0]
			y = arg[1]
			z = arg[2]
		elif len(_np.shape(args[0])) == 2: # if multiple sets (x,y,z) are given (list of lists)
			x = arg[:, 0]
			y = arg[:, 1]
			z = arg[:, 2]
	else:
		raise ValueError("Wrong number of coordinates provided")
	try: 
		len(x)
	except TypeError:
		x = array_fun([x])
		y = array_fun([y])
		z = array_fun([z])
	# Perform actual transformation
	r = sqrt_fun(pow_fun(x,2)+pow_fun(y,2)+pow_fun(z,2))
	theta = arctan2_fun(sqrt_fun(pow_fun(x,2)+pow_fun(y,2)),z)
	phi = arctan2_fun(y,x)
	if not rad:
		theta = rad2deg(theta)
		phi = rad2deg(phi)
	if len(r) > 1:
		return array_fun([r, theta, phi]).transpose()
	else:
		return array_fun([r, theta, phi]).transpose()[0]

def spher2cart(*args, rad=True):
	""" Transforms spherical into cartesian coordinates. 
	
	:param *args: A single or multiple sets of spherical coordinates r, theta, phi. Coordinates can be provided as a single argument, i.e.  [r, theta, phi] or  [[r1,theta1,phi1], [r2,theta2,phi2]] or as three separate arguments, i.e.  r, theta, phi or  [r1, r2], [theta1, theta2], [phi1 , phi2].
	:returns: The corresponding cartesian coordinates.
	"""
	# Symbolic required?
	if all(backends.get_backend(a) != 'sympy' for a in args):
		mode = "numeric"
	else:
		mode = "symbolic"  
	sin_fun = backends.get_calcmethod("sin", mode)
	cos_fun = backends.get_calcmethod("cos", mode)
	multelem_fun = backends.get_calcmethod("multiply", mode)
	array_fun = backends.get_calcmethod("array", mode) 
    # Handle input structure
	if len(args) == 3:
		r   = array_fun(args[0])
		th  = array_fun(args[1])
		phi = array_fun(args[2])
	elif len(args) == 1:
		arg = array_fun(args[0])
		if len(_np.shape(arg)) == 1:
			r = arg[0]
			th = arg[1]
			phi = arg[2]
		elif len(_np.shape(arg)) == 2:
			r = arg[:, 0]
			th = arg[:, 1]
			phi = arg[:, 2]
	else:
		raise ValueError("Wrong number of coordinates provided")
	try:
		len(r)
	except TypeError:
		r = array_fun([r])
		th = array_fun([th])
		phi = array_fun([phi])	    
    # Perform actual transformation
	if not rad:
		th = rad2deg(th)
		phi = rad2deg(phi)
	x = multelem_fun(r, multelem_fun(sin_fun(th),cos_fun(phi)))
	y = multelem_fun(r, multelem_fun(sin_fun(th),sin_fun(phi)))
	z = multelem_fun(r, cos_fun(th))
	if len(x) > 1:
		return array_fun([x,y,z]).transpose()
	else:
		return array_fun([x,y,z]).transpose()[0]


###########################################################
# TRIVIAL FUNCTIONS TO HANDLE LISTS;DICTS; etc. 
###########################################################

def flatten(S): 
	""" Flattens a nested list of lists, i.e. makes a single list out of it.
	
	:param list S: A nested list of lists.
	:returns: The flattened list.
	"""
	lookfor = type(S)
	if not isinstance(S, lookfor):
		if lookfor == list:
			return [S]
		elif lookfor == tuple:
			return (S)
	if len(S) == 0:
		return S
	if isinstance(S[0], lookfor):
		return flatten(S[0]) + flatten(S[1:])
	return S[:1] + flatten(S[1:])

def totalflatten(S): 
	""" Flattens nestest lists and tuples into a single, flat list.
	
	:param list S: A nested list or tuples of lists and/or tuples.
	:returns: The flattened list.
	"""
	if len(S) == 0:
		return S
	if isinstance(S[0], tuple) or isinstance(S[0], list):
		return list(totalflatten(S[0])) + list(totalflatten(S[1:]))
	return list(S[:1]) + list(totalflatten(S[1:]))

def fuse_dictionaries(*args):
	""" Fuses dictionaries into a single dictionary.

	:param *args: A list of dictionaries or a series of dictionary arguments:

	:returns: A single, fused dictionary.
	"""	
	if len(args) == 1:
		arg = args[0]
	else:
		arg = [*args]
	
	if isinstance(arg, dict):
		return arg
	else:
		tmp = arg[0].copy()
		for i in range(1,len(arg)):
			tmp.update(arg[i])
		return tmp
    

###########################################################
# Handle labframe simulations 
###########################################################


def mafm(T_in,f): 
    """Takes a given time and rounds the value to the multiples of the period of a given frequency

    :param T_in : Time, scalar or array
    :param f: Frequency (cycles) to take for the rounding
    """
    d = 1/f
    M = T_in/d
    M_round = _np.round(M)
    if M_round.shape == ():
        return (M_round*d)
    else:
        M_start = M_round[0]
        M_stop = M_round[-1]
        step = _np.ceil(_np.abs((M_stop-M_start)/len(T_in)))*_np.sign(M_stop-M_start)
        M_arr = _np.arange(M_start,M_stop,step)
        return (M_arr*d)    


def t2p(total_time,f):
    """Takes a time and returns the corresponding phase at a given frequency
	"""
    return (2*_np.pi*f*total_time) % (2*_np.pi)



######################################################
# Calculate Electric Fields of Charges 
######################################################

def efield(*args, eps_rel = 1, mode='cart'):
    """ Calculates the electric field of individual charges as well as total electric field in [V/m]
    Vectorized.
    """
    # Argument structure handling
    if len(args) == 4:
        coords = _np.array([args[0], args[1], args[2]]).transpose()
        q = args[3]
    elif len(args) == 2:
        coords = _np.array(args[0])
        q = args[1]
    else:
        raise ValueError("Wrong number of coordinates provided")
    if not isinstance(q, _np.ndarray):
         q = _np.array(q)
    # Spher2cart if necessary
    if mode =='spher':
        coords = spher2cart(coords)
    if len(coords.shape) == 1: # pack if single vector provided 
        coords = _np.array([coords])
    # Calculate Efield
    normalizer = _np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)
    pre = q/(4*_np.pi*eps_0*eps_rel)
    E =  _np.einsum('a, ab -> ab', pre/normalizer**3, coords)
    Etot = _np.sum(E, axis = 0)
    return E, Etot



################################
# Mimic Shot Noise
################################


def _noise_gen(xi):
    return _np.random.poisson(xi) #+ xi
add_shot_noise = _np.vectorize(_noise_gen)


####################
# Rotate 3x3 Matrix 
####################


def rotate_matrix(mat, *args, inverse = False, **kwargs):
    """ Rotates a matrix by applying a "sandwich product" with a rotation matrix. The rotation matrix can be specified using a scipy Rotation object, a sympy Quaternion object or three euler angles.

	Do not use this routine for quantum objects. Use .transform() or rotate_operator() instead.

	:param mat: Input matrix.
	:param *args: Rotation descripton: 
	 * Scipy Rotation object (scipy.spatial.transform.Rotation)
	 * 3 euler angles. If euler angles are provided, the keyword argument convention=... must be given (e.g. convention="zyx").
	 * sympy Quaternion object (sympy.quaternion.Quaternion)

	:Keyword Arguments:
        * *inverse* (''bool'') -- If True, rotation is inverted.
		* *convention* (''bool'') -- e.g. "zyz" or "zxy", only required if euler angles are provided, axes convention. 
	
	:Example:
	>>> rotate_matrix([[0,0,0],[1,0,0],[0,0,0]], np.pi/2, 0, 0, convention = "zxy")
	>>> R = Rotation.from_euler("zxy", [np.pi/2, 0, 0])
	>>> rotate_matrix([[0,0,0],[1,0,0],[0,0,0]], R)

	:returns: The rotated matrix. 
    """

    # Symbolic required?
    if backends.get_backend(mat) == 'sympy' or any(backends.get_backend(x) == "sympy" for x in args):
        calcmethod = "symbolic"
    else:
        calcmethod = "numeric"
    array_fun = backends.get_calcmethod("array", calcmethod)
	# Correct mat datatype.
    if hasattr(mat, "dims"):
        raise ValueError("This is not a quantum object method. Input should be a standard matrix / vector.")
    mat = array_fun(mat)
    # If euler angles provided, construct rotation object.
    if len(args) == 3 or (len(args) == 1 and _np.shape(args[0]) != ()):
        if "convention" not in kwargs:
            raise ValueError("Please specify axes ordering for euler angles via the convention keyword argument.")
        else:
            fromeuler_fun = backends.get_calcmethod("fromeuler", calcmethod)
            if len(args) == 3:
                R = fromeuler_fun(kwargs.get("convention"), [args[0], args[1], args[2]])
            else:
                R = fromeuler_fun(kwargs.get("convention"), args[0])
    # Single argument can  be a rotation object.
    elif len(args) == 1:
        R = args[0]
        if backends.get_backend(R) != "sympy" and calcmethod == "symbolic":
            raise ValueError("Scipy rotation input cannot handle symbolic matrix input. Use euler angles or Sympy Quaternion instead or provide matrix as a non-symbolic object.")
    else:
        raise ValueError("Invalid input for Rotation.")	
    # Perform the actual rotation. 
    rotate_fun = backends.get_calcmethod("rotate",  calcmethod)
    return rotate_fun(mat, R, inverse)
