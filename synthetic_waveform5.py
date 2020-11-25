# Synthetic waveform by using reflectivity method
# Written by Aulia Khalqillah,S.Si.,M.Si

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import random 

# Function of Ricker
def ricker(f,mint,maxt,stept):
    t = np.arange(mint,maxt,stept)
    # convert degrees to radians
    pii = np.pi/180
    r = (1-(2*(pii**2*f**2*t**2)))*np.exp(-pii**2*f**2*t**2)
    return t, r

# Initial parameters
dt = 0.01
mint = 0
maxt = 1800
N = round((maxt-mint)/dt)
f_ricker = 40
nlayers = 3
nsources = 5000
phase = 45 # in degree

# Set depth by using velmod velocity model
velmod = pd.read_csv('AK135F_AVG.csv')
# velmod = pd.read_csv('PREM_1s.csv')
# velmod = pd.read_csv('velmod.csv')
Velocity = 'Vp'
th = np.zeros((velmod.shape[0],1))
for i in range(velmod.shape[0] - 1):
    th[i] = velmod['Depth'].iloc[i+1] - velmod['Depth'].iloc[i]

# Save thicknes to data frame of velmod
velmod['Thick'] =  th
print('\nVELMOD')
print(velmod.head())

# Calculate coefficient reflectivity
cr = np.zeros((nlayers-1,1)) 
for i in range(nlayers-1):
    cr[i] = (velmod['Density'].iloc[i+1] * velmod[Velocity].iloc[i+1] - velmod['Density'].iloc[i] * velmod[Velocity].iloc[i])/(velmod['Density'].iloc[i+1] * velmod[Velocity].iloc[i+1] + velmod['Density'].iloc[i] * velmod[Velocity].iloc[i])

print('\nCR')
print(cr)

# Calculate travel time between layers
td = np.zeros((nlayers,1))
d = 0 
for i in range(nlayers):
    sumdepth = d + velmod['Thick'].iloc[i]
    td[i] = (2*sumdepth)/velmod[Velocity].iloc[i]
    d = sumdepth

print('\nT0')
print(td)

# Calculate travel time that moving on x axis for each layers
xx = 0 # difference distances between source and recivier
sx = np.zeros((nsources,1))
twt = np.zeros((nlayers,nsources))
for j in range(nsources):
    t0 = 10
    for i in range(nlayers):
        t0 = t0 + td[i]
        twt[i,j] = np.sqrt(t0**2 + (xx**2/velmod[Velocity].iloc[i]**2))
    sx[j] = xx  
    random.seed(j)  
    xx = xx + (random.random())

print('\nTWT')
print(twt)
print(len(twt[0,:]))

print('\nsources')
print(len(sx))
print('\n')

# Check if the TWT exceed from maximum of time (maxt)
for j in range(nsources):
    for i in range(nlayers):
        if (twt[i,j] > maxt):
            twt[i,j] = maxt

print(len(twt[0,:]))

# Take index of RC
n = np.zeros((nlayers,nsources))
maxindex = np.zeros((nlayers,1))
for i in range(nlayers):
    for j in range(nsources):
        n[i,j] = int(twt[i,j]/dt)
    print('Max index layer-%d is %d' % (i+1,int(max(n[i,:]))))
    maxindex[i] = int(max(n[i,:]))

print('\nLENGTH OF INDEX RC')
print(len(n[0,:]))

print('\nMAX INDEX')
print(maxindex)

# Apply reflectivity to time series for each sources base don its index
RR = np.zeros((N+1,nsources))
for j in range(nsources):
    for i in range(nlayers):
        RR[int(n[i,j]),j] = cr[i-1]

# Sum the reflectivity from all sources
R = 0
for i in range(nsources):
    R = R + RR[:,i]

R = R[:-1]
print('\nR') 
print(len(R))
    
# Generate ricker wavelet
tricker, wavelet = ricker(f_ricker,-10,10,0.01)

# Apply phase to wavelet
phase = np.cos(phase*(np.pi/180))
wavelet = wavelet*phase

# Convolution between wavelet and reflectivity
trace_conv = np.convolve(wavelet,R,mode='same')

# Generate time series
final_time = np.arange(0,(len(R)*dt),dt)

plt.figure(1,figsize=(10,6))
plt.subplot(3,1,1)
plt.plot(tricker, wavelet,linewidth=0.5)
plt.grid()
plt.title('Ricker Wavelet %d Hz' % f_ricker)
plt.xlabel('Time (s)')
plt.ylabel('Magnitude')
# plt.gca().invert_yaxis()
plt.tight_layout(pad=2)

plt.subplot(3,1,2)
plt.plot(final_time, R,linewidth=0.5)
plt.grid()
plt.title('Coefficient Reflectivity')
plt.xlabel('Time (s)')
plt.ylabel('Magnitude')
# plt.gca().invert_yaxis()
plt.tight_layout(pad=2)

plt.subplot(3,1,3)
plt.plot(final_time, trace_conv,linewidth=0.5)
plt.grid()
plt.title('Synthetic Waveform')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
# plt.gca().invert_yaxis()
plt.tight_layout(pad=2)
plt.show()
