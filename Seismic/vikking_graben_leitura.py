import numpy as np
import matplotlib.pyplot as plt
import segyio
from scipy.ndimage.filters import uniform_filter1d

class Seismic():
    
    def __init__(self, path):
        self.sismica = segyio.open(path, "r", ignore_geometry=True)

    def header(self, n):
        a = dict(self.sismica.header[n])
        b = list(a.keys())
        for i in range(len(b)):
            b[i] = f"{b[i]}"

        self.head = dict(zip(b, a.values()))    
            
    def single_trace(self, trace):
        
        
        _, axs = plt.subplots(nrows= 1, ncols = 2, figsize=(8,16))

        points = np.arange(0, len(trace), 1)
        x = trace.T
        x = x/np.max(x)
        fig1 = axs[0].plot(x , points,'k', linewidth=0.7)
        axs[0].fill_betweenx(points, 0, x, where= x>0, color="k")
        axs[0].invert_yaxis()

        freqs = np.linspace(0, 1/0.004, len(trace))
        freq = np.fft.fft(trace)
        fig2 = axs[1].plot(freqs, np.abs(freq)/np.max(np.abs(freq)), )
        axs[1].set_xlim(0, 130)

        plt.show()
    
    def moving_average(self, trace,x , peso):
        ret = uniform_filter1d(trace[:x+1], size=peso, mode = "nearest")
        self.average = ret
    
    def agc_rms(self, trace, gate, desired):
        self.AGCRMS = [0]*len(trace[:241])
        for each in range(len(trace[:241])):
            sub = np.array_split(trace[each], len(trace[each])/(gate/4))
            for subtrace in range(len(sub)):
                rms = np.sqrt(
                              np.sum(
                                     np.square(
                                               sub[subtrace]
                                              )
                                    ) 
                                    /(len(sub[subtrace]))
                             )

                g = desired/rms
                sub[subtrace] = sub[subtrace]*g
            sub = np.concatenate(sub)
            self.AGCRMS[each] = sub
        self.AGCRMS = np.array(self.AGCRMS)

    def show_seismic(self, trace, clip = 1e+1):
        
        vmin, vmax = -clip, clip

        
        figsize=(10, 10)
        _, axs = plt.subplots(nrows=1, ncols=1, figsize=figsize, facecolor='w', edgecolor='k',
                            squeeze=False,
                            sharex=True)
        axs = axs.ravel()
        im = axs[0].imshow(trace[:241].T, cmap='gray', vmin=vmin, vmax=vmax)
        plt.show()
        
        
        

        
            
    

        

        

    



path = r"Seismic\seismic.segy"

S = Seismic(path)
""" S.single_trace(S.sismica.trace.raw[0])
S.moving_average(S.sismica.trace.raw, 1, 7)
S.single_trace(S.average[0]) """
S.agc_rms(S.sismica.trace.raw, 64, 100)
S.single_trace(S.AGCRMS[0])
S.show_seismic(S.sismica.trace.raw, 1e+2)
S.show_seismic(S.AGCRMS, 1e+2)













S.sismica.close()