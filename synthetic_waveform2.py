import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import random 

def ricker(f,mint,maxt,stept):
    t = np.arange(mint,maxt,stept)
    # convert degrees to radians
    pii = np.pi/180
    r = (1-(2*(pii**2*f**2*t**2)))*np.exp(-pii**2*f**2*t**2)
    return t, r

def harmonic(A,f,vel,disp,mint,maxt,stept):
    t = np.arange(mint,maxt,stept)
    w = 2*(np.pi/180)*f
    k = w/vel
    # y = A*np.sin(2*(np.pi/180)*f*t)
    y = A*np.exp(-1j*(w*t) - (k*disp))
    return t,y

# Initial parameters
dt = 0.01
mint = 0
maxt = 200
N = round((maxt-mint)/dt)
v = 2500
A = 20
f = 10
f_ricker = 40
x = 0
nlayers = 10

# Generate harmonic wave of displacement
time, displacement = harmonic(A,f,v,x,mint,maxt,dt) 

# Set depth by using velmod velocity model
velmod = pd.read_csv('AK135F_AVG.csv')

# Calculate coefficient reflectivity
nsource = 10
cr = np.zeros((nlayers-1,nsource)) 
coeref = []
for i in range(nlayers-1):
    for j in range(nsource):
        cr[i,j] = (velmod['Density'].iloc[i+1] * velmod['Vp'].iloc[i+1] - velmod['Density'].iloc[i] * velmod['Vp'].iloc[i])/(velmod['Density'].iloc[i+1] * velmod['Vp'].iloc[i+1] + velmod['Density'].iloc[i] * velmod['Vp'].iloc[i])
        # coeref.append(cr[i,j])

# Calculate travel time between layers
td = np.zeros((nlayers-1,1))
d = 0 
for i in range(nlayers-1):
    sumdepth = d + velmod['Depth'].iloc[i]
    td[i] = (2*sumdepth)/velmod['Vp'].iloc[i]
    d = sumdepth

# Calculate travel time between layers
t0 = 0
xx = 0
twt = np.zeros((nlayers-1,nsource))
for j in range(nsource):
    for i in range(nlayers-1):
        twt[i,j] = np.sqrt(td[i]**2 + (xx**2/velmod['Vp'].iloc[i]**2))
        # t0 = twt[i,0]
    # xx = xx + (random.random())
    xx = xx + 10

# Calculate reflectivity
idx = []
n = np.zeros((nlayers-1,nsource))
for i in range(nlayers-1):
    for j in range(nsource):
        n[i,j] = int(twt[i,j]/dt) 

R = np.zeros((N,nsource))
for i in range(nsource):
    for k in range(nlayers-1):
        R[int(n[k,i]),i] = cr[k,i]

plt.figure(1,figsize=(10,5))
for i in range(nsource):
    plt.subplot(1,nsource,i+1)
    plt.plot(R[:,i],time,linewidth=0.5)
    plt.grid()
    plt.gca().invert_yaxis()
plt.tight_layout(pad=1)
plt.show()

# Generate ricker wavelet
tricker, wavelet = ricker(f_ricker,-10,10,0.01)

# Convolution between wavelet and reflectivity
trace_conv = np.zeros((N,nsource))
for j in range(nsource):
    trace_conv[:,j] = np.convolve(wavelet,R[:,j],mode='same')

plt.figure(2,figsize=(10,5))
for i in range(nsource):
    plt.subplot(1,nsource,i+1)
    plt.plot(trace_conv[:,i],time,linewidth=0.5)
    plt.grid()
    plt.gca().invert_yaxis()
plt.tight_layout(pad=1)
plt.show()

# # Convolution between trace_conv
# trace = trace_conv[:,0]
# for i in range(nsource-1):
#     # trace = np.convolve(trace,trace_conv[:,i+1],mode='same')
#     trace = trace * trace_conv[:,i+1]

# plt.figure(3,figsize=(3,6))
# plt.plot(trace,time,linewidth=0.5)
# plt.gca().invert_yaxis()
# plt.show()

# # Cross-correlate
# # disp = np.convolve(trace,displacement,mode='same')
# disp = displacement * trace

# plt.figure(4,figsize=(3,6))
# plt.plot(np.real(disp),time,linewidth=0.5)
# plt.gca().invert_yaxis()
# plt.show()