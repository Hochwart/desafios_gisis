import numpy as np
import matplotlib.pyplot as plt
import segyio
from scipy.ndimage.filters import uniform_filter1d
from scipy.signal import butter, lfilter, hilbert

class Seismic():
    """ Class for .segy seismic trace handling and visualizing. \n
    Input: \n
    path = path for the .segy archive"""
    def __init__(self, path):
        
        self.sismica = segyio.open(path, "r", ignore_geometry=True) #Abre o arquivo com segyio

    def header(self, n):
        """ Trace header. \n
        Input: \n
        n = trace index \n
        Output: \n
        A dictionary with the desired trace header"""

        a = dict(self.sismica.header[n])
        b = list(a.keys())
        for i in range(len(b)):
            b[i] = f"{b[i]}"

        self.head = dict(zip(b, a.values()))    
            
    def single_trace(self, trace):
        """Single trace amplitudes and FFT plot. \n
        Input: \n
        trace = desired trace \n
        Output: \n
        Single trace and trace frequencies plot."""

        _, axs = plt.subplots(nrows= 1, ncols = 2, figsize=(8,16))

        points = np.arange(0, len(trace), 1)
        x = trace.T
        x = x/np.max(x)
        fig1 = axs[0].plot(x , points,'k', linewidth=0.7)
        axs[0].fill_betweenx(points, 0, x, where= x>0, color="k", interpolate=True)
        axs[0].invert_yaxis()

        freqs = np.linspace(0, 1/0.004, len(trace))
        freq = np.fft.fft(trace)
        fig2 = axs[1].plot(freqs, np.abs(freq)/np.max(np.abs(freq)))
        axs[1].set_xlim(0, 130)

        plt.show()
    
    def moving_average(self, trace, weight):
        """Moving mean filter. \n
        Input: \n
        trace  =  Desired traces \n
        weight = Weight of the moving mean. \n
        Output: \n
        An array for the moving mean of the desired traces."""

        ret = uniform_filter1d(trace, size=weight, mode = "nearest")
        self.average = ret
    
    def trace_rms(self, trace, gate, desired):
        """ Trace RMS Filter. \n
        Input: \n 
        trace = Desired Traces \n
        gate  = Desired gate for the RMS equation \n
        desired = Desired Multiplier for the RMS equation \n
        Output: \n
        An array for the AGC RMS of the desired traces."""
        self.rms = [0]*len(trace)
        
        
        for each in range(len(trace)):
            sub = np.zeros(len(trace[each]))
            parts = int(len(trace[each])/(gate/4))
            for points in range(len(trace[each])):
                
                i = points+parts
                rms_ = ((1/parts)*np.mean(trace[each][points:i]**2))**0.5
                sub[points] = (desired/rms_)
                
            self.rms[each] = trace[each]*sub
        self.rms = np.array(self.rms)
    
    def trace_integration(self, trace):
        """ Trace Integration Filter. \n
        Input: \n 
        trace = Desired Traces \n
        Output: \n
        An array for the Integral of the desired traces."""


        self.integ = [0]*len(trace)
        for each in range(len(trace)):
            self.integ[each] = np.cumsum(trace[each])
        self.integ = np.array(self.integ)

    def butter_bandpass(self, lowcut, highcut, fs, order):
        """ Band Pass Filter. \n
        Input: \n 
        trace = Desired Traces \n
        lowcut = low frequency cut \n
        highcut = high frequency cut \n
        fs = sampling frequency \n
        order = butterworth order \n
        Output: \n
        An array for the band pass filter of the desired traces."""
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        self.b, self.a = butter(order, [low, high], btype='band')
        


    def bandpass(self, trace, lowcut, highcut, fs, order):
        """ Band Pass Filter. \n
        Input: \n 
        trace = Desired Traces \n
        lowcut = low frequency cut \n
        highcut = high frequency cut \n
        fs = sampling frequency \n
        order = butterworth order \n
        Output: \n
        An array for the band pass filter of the desired traces."""

        self.filtered = [0]*len(trace)
        self.butter_bandpass(lowcut, highcut, fs, order)
        
        for each in range(len(trace)):
            self.filtered[each] = lfilter(self.b, self.a, trace[each])
        self.filtered = np.array(self.filtered)

    def rot_phase(self, trace):
        """ 90º degree Phase Rotation Filter. \n
        Input: \n 
        trace = Desired Traces \n
        Output: \n
        An array for the phase rotation filter of the desired traces."""

        self.phase_shift = [0]*len(trace)
        for each in range(len(trace)):
            xa = hilbert(trace[each])
            a = (np.imag(xa)**2+np.real(xa)**2)**0.5
            self.phase_shift[each] = -1*np.imag(xa)
        self.phase_shift = np.array(self.phase_shift)

    def TecVA(self, trace, lowcut, highcut, fs, order):
        """ TecVA Filter. \n
        Input: \n 
        trace = Desired Traces \n
        lowcut = low frequency cut \n
        highcut = high frequency cut \n
        fs = sampling frequency \n
        order = butterworth order \n
        Output: \n
        An array for the phase rotation filter of the desired traces."""

        self.trace_integration(trace)
        self.bandpass(self.integ, lowcut, highcut, fs, order)
        self.trace_rms(self.filtered, 6, 1)
        self.bandpass(self.rms, lowcut, 120, fs, order)
        self.rot_phase(self.filtered)
        self.tecva = self.phase_shift

    def show_seismic(self, trace, clip = 1e+1):
        """ Method for imaging the seismic traces.\n
        Input:  \n
        trace = Desired traces \n
        clip  = Desired clip for contrast (1e+1 Default) \n
        Output: \n
        Seismic image for the desired traces."""
        vmin, vmax = -clip, clip

        
        figsize=(10, 8)
        _, axs = plt.subplots(nrows=1, ncols=1, figsize=figsize, facecolor='w', edgecolor='k',
                            squeeze=False,
                            sharex=True)
        axs = axs.ravel()
        im = axs[0].imshow(trace.T, cmap="gray", vmin=vmin, vmax=vmax)
        plt.show()
        

path = r"Seismic\seismic.segy"

S = Seismic(path)

#Usando os filtros
#S.moving_average(S.sismica.trace.raw[:361], 7)  #Média Móvel
#S.trace_rms(S.sismica.trace.raw[::10][:2500], 6, 1)                  #RMS
#S.trace_integration(S.sismica.trace.raw[::10])
#S.bandpass(S.sismica.trace.raw[::10], 6, 40, 1/0.004, 10)
#S.rot_phase(S.sismica.trace.raw[::10])
S.TecVA(S.sismica.trace.raw[:2], 6, 40, 1/0.004, 10)
S.single_trace(S.sismica.trace.raw[0])
S.single_trace(S.tecva[0])






""" #Vendo os traços
S.single_trace(S.sismica.trace.raw[0])          
S.single_trace(S.average[0])
S.single_trace(S.AGCRMS[0]) 
S.single_trace(S.sismica.trace.raw[0]) 
S.single_trace(S.band_pass[0])"""

#Vendo a sísmica
#S.show_seismic(S.sismica.trace.raw[:361], 1e+1)
#S.show_seismic(S.average, 1e+2)
#S.show_seismic(S.sismica.trace.raw[::10], 1e+3)
#S.show_seismic(S.rms, 1e+2)
#S.show_seismic(S.integ, 1e+4)
#S.show_seismic(S.filtered, 1e+2)
#S.show_seismic(S.phase_shift, 1e+2)
#S.show_seismic(S.tecva, 1e+2)


S.sismica.close()