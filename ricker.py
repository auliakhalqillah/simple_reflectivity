import numpy as np
import math as m  
import matplotlib.pyplot as plt 

def ricker(f,mint,maxt,stept):
    t = np.arange(mint,maxt,stept)
    # convert degrees to radians
    pii = np.pi
    r = (1-(2*(pii**2*f**2*t**2)))*np.exp(-pii**2*f**2*t**2)
    return t, r

def fourier_ricker(f,minf,maxf,stepf):
    pii = np.pi
    w = 2*pii*np.arange(minf,maxf,stepf)
    wp = 2*pii*f
    fr = ((2*w**2)/(np.sqrt(pii)*wp**3))*np.exp((-w**2/wp**2))
    return w, fr

f = 30
start = -0.5
stop = 0.5
step = (stop-start)/1000
time, rick = ricker(f,start,stop,step)
freq, fourier = fourier_ricker(f,0,80,0.5)

plt.figure(1)
plt.subplot(2,1,1)
plt.plot(time,rick)
plt.xlabel('Time(s)')
plt.ylabel('Magnitude')
plt.title('Ricker Wavelet')

plt.subplot(2,1,2)
plt.plot(freq,fourier)
plt.xlabel('\u03C9(rad/s)')
plt.ylabel('R(\u03C9)')
plt.title('Fourier Transform of Ricker Wavelet')
plt.show()