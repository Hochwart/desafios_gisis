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
        
        
        fig, axs = plt.subplots(nrows= 1, ncols = 2, figsize=(8,16))

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
        ret = uniform_filter1d(trace[:x], size=peso, mode = "nearest")
        self.average = ret
            
        
            
    

        

        

    



path = r"Seismic\seismic.segy"

S = Seismic(path)
""" S.single_trace(S.sismica.trace.raw[0]) """
S.moving_average(S.sismica.trace.raw, 12, 7)
""" S.single_trace(S.average[0]) """

fig, axs = plt.subplots(nrows= 2, ncols = 2, figsize=(8,16))

points = np.arange(0, len(S.sismica.trace.raw[0]), 1)
x = S.sismica.trace.raw[0].T
x = x/np.max(x)
fig1 = axs[0][0].plot(x , points,'k', linewidth=0.7)
axs[0][0].fill_betweenx(points, 0, x, where= x>0, color="k")
axs[0][0].invert_yaxis()

freqs = np.linspace(0, 1/0.004, len(S.sismica.trace.raw[0]))
freq = np.fft.fft(S.sismica.trace.raw[0])
fig2 = axs[0][1].plot(freqs, np.abs(freq)/np.max(np.abs(freq)), )
axs[0][1].set_xlim(0, 130)

points = np.arange(0, len(S.average[0]), 1)
x = S.average[0].T
x = x/np.max(x)
fig3 = axs[1][0].plot(x , points,'k', linewidth=0.7)
axs[1][0].fill_betweenx(points, 0, x, where= x>0, color="k")
axs[1][0].invert_yaxis()

freqs = np.linspace(0, 1/0.004, len(S.average[0]))
freq = np.fft.fft(S.average[0])
fig4 = axs[1][1].plot(freqs, np.abs(freq)/np.max(np.abs(freq)), )
axs[1][1].set_xlim(0, 130)


plt.show()






S.sismica.close()