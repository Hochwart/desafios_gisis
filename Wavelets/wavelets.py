import numpy as np
import matplotlib.pyplot as plt


### Adicionar o gráfico como um método diferente
### Abrir o dt na self
### Separar a transformada para o gráfico das freq

class Wavelets:
    def __init__(self):           

        #Parametros comuns para todas as funções
        self.N      = 1000
        self.t      = np.linspace(-200/self.N, 200/self.N, self.N , endpoint="True")  
        self.freqs  = 1/(self.t[1]-self.t[0])

    def klauder(self):                                                 
        #Setando as variáveis para o cálculo


        k   = (self.f[1]-self.f[0])/self.T
        f0  = (self.f[1]+self.f[0])/2

        #Fórmula de Klauder, a transformada e normalização
        wl       = np.real(np.sin(np.pi* k * self.t * 
                          (self.T-self.t))/(np.pi*k*self.t) * 
                          np.exp(2* np.pi * 1j * f0 * self.t))
        self.fwl = np.fft.fft(wl)/np.max(wl)
        self.wl  = wl/np.max(wl)
        self.title    = f"Klauder Wavelet at {self.f[0]}Hz, {self.f[1]}Hz and {self.T}s "

    def ormsby(self):                                             
                                            


        #Fórmula de Ormsby, a transformada e normalização
        wl  = (((np.pi*self.f[3])**2)/(np.pi*self.f[3]-np.pi*self.f[2]) * 
                np.sinc(self.f[3]*self.t)**2 - ((np.pi*self.f[2])**2)/
                (np.pi*self.f[3] - np.pi*self.f[2]) * np.sinc(self.f[2]*self.t)**2) - (
                ((np.pi*self.f[1])**2)/(np.pi*self.f[1] - np.pi*self.f[0]) * 
                np.sinc(self.f[1]*self.t)**2 - 
                ((np.pi*self.f[0])**2)/(np.pi*self.f[1] - np.pi*self.f[0]) * 
                np.sinc(self.f[0]*self.t)**2)

        self.fwl = np.fft.fft(wl)/np.max(wl)
        self.wl  = wl/np.max(wl)
        self.title    = f"Ormsby Wavelet at {self.f[0]}Hz, {self.f[1]}Hz, {self.f[2]}Hz and {self.f[3]}Hz"


    def rickter(self):                                                        
        #Formula de Ricker, a transformada e normalização
        wl       = (1- 2* (np.pi**2) * (self.f**2) * (self.t**2) ) * np.exp(-1*(np.pi**2)*(self.f**2)* (self.t**2))
        self.fwl = np.fft.fft(wl)/max(wl)
        self.wl  = wl/np.max(wl) 
        self.title    = f"Rickter Wavelet at {self.f}Hz"

    def sum_sine(self):
        wl = 0
        for i in range(len(self.f[0])):
            wl += self.f[0][i] * np.sin( 2 * np.pi * self.f[1][i] * 
                                        self.t + self.f[2][i] * (np.pi/180) ) 
        self.fwl = np.fft.fft(wl)/max(wl)    
        self.wl = wl/np.max(wl)
        self.title = "Sum of Sines"


      
    
    def graph(self):

        plt.figure(figsize=(8,8))
        #Gráfico da Wavelet
        plt.subplot(2, 1, 1).set_title(f'Amplitudes for {self.title}')
        plt.plot(self.t*1000, self.wl)
        plt.grid(axis="y")
        plt.xlabel("Time (ms)")

        #Gráfico das freq
        plt.subplot(2, 1, 2).set_title(f"Spec. Density for {self.title}")
        freq = np.linspace(0.0, self.freqs, self.N)
        plt.plot(freq, np.abs(self.fwl)/(np.abs(np.max(self.fwl))))
        plt.xlim(0, 140)
        plt.xlabel("Frequency (Hz)")
        plt.grid()

        plt.tight_layout(pad=4.0)
        plt.show()


# Main
wave = Wavelets()

# Rickter
""" 
wave.f = 15
wave.rickter()
wave.graph()
 """
# Ormsby
""" 
wave.f = [10,15,40,45]
wave.ormsby()
wave.graph()
 """
# Klauder
""" 
wave.f = [15, 40]
wave.T = 5
wave.klauder()
wave.graph()
"""
# Sine Sum

wave.f = [[1, 1, 1], [45, 90, 135], [25, 10, 30]]
wave.sum_sine()
wave.graph()