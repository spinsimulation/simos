import numpy as np
import collections
import scipy.signal as scisig

def ffthelper(y,t,K=0,dc_block=False,scale='linear',real_input=False, power_spec=True):
    """Fourier transform the signal

    Returns [frequency axis, power-spectrum]

    Parameters:
    y : Signal which should be Fourier transformed
    t : Corresponding time axis or sampling interval

    Optional paramaters:
    K          : Undersampling factor. By default 0
    dc_block   : Removes the average of the y trace before Fourier transforming it
    scale      : By default 'linear'. Possible values: 'linear', 'dB'
    real_input : By default False. If set to true, a spectrum only with positive frequencies is returned.
    power_spec : Default True. It False, the complex fft is returned.

    For the dB scale, the spectrum is cliped at a lower bound of -200 dB (to prevent possible divergence).
    """
    y = np.array(y)
    if dc_block:
        y = y - np.average(y)
    fs = 0
    N  = len(y)
    try:
        len(t)
        fs = 1/(t[1]-t[0])
    except TypeError:
        fs = 1/t
    
    padded = np.append(y, np.zeros(K*len(y)))
    N = N+K*N
    if real_input:
        fax = np.fft.rfftfreq(N, 1/fs)
        spec = np.fft.rfft(padded)
        if power_spec:
            spec = np.abs(spec)**2
    else: 
        fax = fs*np.arange(0,N)/N-fs/2
        spec = np.fft.fftshift(np.fft.fft(padded))
        if power_spec:
            spec = np.abs(spec)**2

    if scale == 'dB':
        spec = np.clip(spec,10**(-20),None)
        spec = 10*np.log10(spec)

    return fax, spec


def _butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = scisig.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def lowpass(y, dt, cutoff, order=5):
    """Lowpass filter the input data with a Butterworth filter. Returns the filtered data.
    
    Parameter:
    y      : data input
    dt     : sampling time
    cutoff : -3dB frequency, where gain drops to 1/sqrt(2)
    order  : the order of the filter
    """
    b, a = _butter_lowpass(cutoff, 1/dt, order=order)
    y = scisig.lfilter(b, a, y)
    return y

def _butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = scisig.butter(order, normal_cutoff, btype='high', analog=True)
    return b, a

def highpass(y, dt, cutoff, order=5):
    """Highpass filter the input data with a Butterworth filter. Returns the filtered data.
    
    Parameter:
    y      : data input
    dt     : sampling time
    cutoff : -3dB frequency, where gain drops to 1/sqrt(2)
    order  : the order of the filter
    """

    b, a = _butter_highpass(cutoff, 1/dt, order=order)
    y = scisig.lfilter(b, a, y)
    return y


def us(f,tau):
    fs = 1/tau
    factor = f/fs
    #print(factor)
    factor = factor % 1
    if factor > 0.5:
        factor = 1.0 - factor
    #print(factor)
    return factor*fs



def ffthelper2d(y,t1,t2,K1=0,K2=0,power_spect=True,dc_block1=False,dc_block2=False,real_input=False):
    """
    time axis generate via tt2,tt1 = np.meshgrid(tax2,tax1)
    plot with pcolormesh(fax2,fax1,fft) (because of NMR convention)
    """
    y = np.array(y)
    
    # Zero Padding
    s = np.array(np.shape(y))
    s[0] = s[0]*(1+K1)
    s[1] = s[1]*(1+K1)
    
    # Extract sampling times if axes are given
    if isinstance(t1,collections.Iterable):
        if len(np.shape(t1)) == 1:
            t1 = t1[1]-t1[0]
        elif len(np.shape(t1)) == 2:
            t1 = t1[1,0] - t1[0,0]
        else:
            raise ValueError
    
    if isinstance(t2,collections.Iterable):
        if len(np.shape(t1)) == 1:
            t2 = t2[1]-t2[0]
        elif len(np.shape(t2)) == 2:
            t2 = t2[0,1] - t2[0,0]
        else:
            raise ValueError
    
    # Remove DC offset
    if dc_block1:
        y = y - np.mean(y,axis=0,keepdims=True)
    if dc_block2:
        y = y - np.mean(y,axis=1,keepdims=True)
        
   
    #Fourier Transform
    if real_input:
        fft = np.fft.rfft(y,n=s[1],axis=1,norm="ortho")
        fft = np.fft.fft(fft,n=s[0],axis=0,norm="ortho")
    else:
        #fft = np.fft.fft(y,n=s[1],axis=1,norm="ortho")
        #fft = np.fft.fft(fft,n=s[0],axis=0,norm="ortho")
        fft = np.fft.fftn(y,s=s,axes=(0,1),norm="ortho")

    
    fft = np.fft.fftshift(fft,axes=(0))
    if not real_input:
        fft = np.fft.fftshift(fft,axes=(1))
    
    
    # Frequency axes
    fax1 = np.fft.fftshift(np.fft.fftfreq(s[0],t1))
    
    if real_input:
        fax2 = np.fft.rfftfreq(s[1],t2)
    else:
        fax2 = np.fft.fftshift(np.fft.fftfreq(s[1],t2))
    
    fax2,fax1 = np.meshgrid(fax2,fax1)
    
    # Normalization
    normfakt = 1/np.sqrt(s[0]*s[1])

    if power_spect: 
        return fax1,fax2, np.abs(fft*normfakt)**2
    else:
        return fax1, fax2, fft
