import numpy as _np


###########################################################
# COMMON FUNCTIONS (Used e.g. for fitting)
###########################################################
def sine(x,A,f,phi):
	"""Return a sinus function (for use with lmfit)

	Parameters:
	x   : time axis
	A   : amplitude
	f   : Frequency (will be converted to angular freq.)
	phi : Phase shift in radians
	"""
	return A*_np.sin(2*_np.pi*f*x+phi)


def lor(x,amplitude,center,fwhm):
	"""Returns the Lorentzian function

	Parameter:
	x           : abscissa
	amplitude   : amplitude (max value of the Lorentzian)
	center      : resonance position
	fwhm        : full-width half-maximum width
	Note: The area under the Lorenzian curve is given by pi*amplitude*fwhm/2
	"""
	w = fwhm/2
	return amplitude/(1+(x-center)**2/w**2)


def dlor(x,amplitude,center,fwhm=None,pp_width=None):
	"""Returns the derivative of a Lorentzian

	Parameter:
	x           : abscissa
	amplitude   : amplitude
	center      : resonance position
	
	provide only 1 of two possible width arguments:   
	fwhm        : full-width half-maximum width
	pp_width    : Peak-to-Peak width
	
	Note that both width are related as fwhm=sqrt(3)*pp_width
	"""
	if (fwhm is None) and (pp_width is None):
		raise ValueError('You must either provide fwhm or pp_width')
	elif (fwhm is not None) and (pp_width is not None):
		raise ValueError('You can provide only one width parameter')
	elif (fwhm is None) and (pp_width is not None):
		w = _np.sqrt(3)*pp_width/2
	elif (fwhm is not None) and (pp_width is None):
		w = fwhm/2
	return -((2*amplitude*(-center + x))/(w**2*(1 + (-center + x)**2/w**2)**2))


def lordisp(x,amplitude,center,fwhm):
	"""Returns the Lorentzian dispersion function

	Parameter:
	x           : abscissa
	amplitude   : amplitude (max value of the Lorentzian)
	center      : resonance position
	fwhm        : full-width half-maximum width
	Note: The area under the Lorenzian curve is given by pi*amplitude*fwhm/2
	"""
	w = fwhm/2
	return amplitude/(1+(x-center)**2/w**2)*(x-center)/w

def dlordisp(x,amplitude,center,fwhm=None,pp_width=None):
	"""Returns the derivative of a Lorentzian dispersion function

	Parameter:
	x           : abscissa
	amplitude   : amplitude
	center      : resonance position
	
	provide only 1 of two possible width arguments:   
	fwhm        : full-width half-maximum width
	pp_width    : Peak-to-Peak width
	
	Note that both width are related as fwhm=sqrt(3)*pp_width
	"""
	if (fwhm is None) and (pp_width is None):
		raise ValueError('You must either provide fwhm or pp_width')
	elif (fwhm is not None) and (pp_width is not None):
		raise ValueError('You can provide only one width parameter')
	elif (fwhm is None) and (pp_width is not None):
		w = _np.sqrt(3)*pp_width/2
	elif (fwhm is not None) and (pp_width is None):
		w = fwhm/2
	return -((2*amplitude*(-center + x)**2)/(w**3*(1 + (-center + x)**2/w**2)**2)) + amplitude/(w*(1 + (-center + x)**2/w**2))

def langevin(x):
	"""Langevin function
	L(x) = coth(x) - 1/x, is zero at x=0.
	Used e.g. in the conext for Magnetic Particle Imaging (MPI). Non-linear response function of a SPION.
	"""
	eps = _np.finfo(x.dtype).eps
	return _np.piecewise(x,[_np.abs(x)<=eps,_np.abs(x)>eps],[0, lambda i: _np.cosh(i)/_np.sinh(i)-1/i])


def apodize_knee(length,roll_off):
    y1 = _np.ones(roll_off)
    y2 = _np.linspace(1,0,length-roll_off)
    return _np.concatenate((y1,y2))


def apodize_sine_bell(length,roll_off=None):
    # https://doi.org/10.1016/j.aca.2008.01.030
    x = _np.linspace(0,_np.pi,length)
    if roll_off is None:
        damping = 0.5
    else:
        damping=roll_off/length
    damping = 1/_np.tan(damping*_np.pi)
    y1 = _np.sin(x)*_np.exp(-x*damping)
    y1 = y1/_np.max(y1)
    y1[-1] = 0
    return y1


def apodize_cos2(length,roll_off,rise=0):
    x1 = _np.linspace(-_np.pi/2,0,rise)
    y1 = _np.cos(x1)**2
    y2 = _np.ones(roll_off-rise)
    x3 = _np.linspace(0,_np.pi/2,length-roll_off)
    y3 = _np.cos(x3)**2
    return _np.concatenate((y1,y2,y3))