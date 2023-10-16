## Implements https://twiki.cern.ch/twiki/pub/LHCbPhysics/B2JPsiPhiK/LHCb-ANA-2015-042.pdf Eqs. 28 - 32

import numpy as np
import matplotlib.pyplot as plt
from mpmath import * 
from cmath import phase

M_A = 2
M_B = 2
beta_AB = 0.3

# set precision
mp.dps = 500

def mu_AB(ma,mb):
    return ma * mb / (ma + mb)

def Z(m,ma,mb,mu,beta):
    return 4 * mu / beta**2  * (ma+mb-m)

def I(x):
    return 0.5 * sqrt(pi) * (1.0 - sqrt(pi*x) * exp(x) * erfc(sqrt(x)) )

def Pi(m,ma,mb,beta):
    mu = mu_AB(ma,mb)
    return - mu*beta/(sqrt(2)*pi**2) * I(Z(m,ma,mb,mu,beta)) 


## I(Z)

Z_values = np.linspace(-10, 10, 500)
I_values = [-I(-Z) for Z in Z_values]

Re_I = [i.real for i in I_values]
Im_I = [i.imag for i in I_values]

# Plot the Re/Im
plt.figure(figsize=(11, 9))
plt.subplot(2, 2, 1)
plt.plot(Z_values, Re_I, label='-Re I(Z)', color='red')
plt.plot(Z_values, Im_I, label='-Im I(Z)', color='blue')
plt.xlabel('-Z')
plt.ylabel('Negative Real and Imaginary parts of -I(Z)')
plt.legend()
plt.grid(True)

# Plot Argand
plt.subplot(2, 2, 2)
plt.plot(Re_I, Im_I, label='-Im I(Z) vs -Re I(Z)', color='green')
plt.xlabel('-Re I(Z)')
plt.ylabel('-Im I(Z)')
plt.xlim([-1, 1])  
plt.ylim([-1, 1])
plt.grid(True)

# Calculate magnitude squared and phase
magnitude_squared = [abs(i)**2 for i in I_values]
phases = [phase(i) for i in I_values]

# Plot Magnitude Squared of I(Z)
plt.subplot(2, 2, 3)
plt.plot(Z_values, magnitude_squared, label='|I(Z)|^2', color='purple')
plt.xlabel('Z')
plt.ylabel('|I(Z)|^2')
plt.grid(True)

# Plot Phase of I(Z)
plt.subplot(2, 2, 4)
plt.plot(Z_values, phases, label='Phase of I(Z)', color='orange')
plt.xlabel('Z')
plt.ylabel('Phase (radians)')
plt.grid(True)

#plt.show()

## Pi(m)
m_values = np.linspace(3, 4.05, 500)
I_values = [Pi(m,M_A,M_B,beta_AB) for m in m_values]

Re_I = [i.real for i in I_values]
Im_I = [i.imag for i in I_values]

# Plot Re/Im
plt.figure(figsize=(11, 9))
plt.subplot(2, 2, 1)
plt.plot(m_values, Re_I, label='Re Pi(m)', color='red')
plt.plot(m_values, Im_I, label='Im Pi(m)', color='blue')
plt.xlabel('m')
plt.ylabel('Real and Imaginary parts of Pi(m)')
plt.legend()
plt.grid(True)

# Plot Argand
plt.subplot(2, 2, 2)
plt.plot(Re_I, Im_I, label='Im Pi(m) vs Pi(m)', color='green')
plt.xlabel('Re Pi(m)')
plt.ylabel('Im Pi(m)')
# plt.xlim([-1, 1])  
# plt.ylim([-1, 1])
plt.grid(True)

# Calculate magnitude squared and phase
magnitude_squared = [abs(i)**2 for i in I_values]
phases = [phase(i) for i in I_values]

# Plot Magnitude Squared 
plt.subplot(2, 2, 3)
plt.plot(m_values, magnitude_squared, label='|Pi(m)|^2', color='purple')
plt.xlabel('m')
plt.ylabel('|Pi(m)|^2')
plt.grid(True)

# Plot Phase 
plt.subplot(2, 2, 4)
plt.plot(m_values, phases, label='Phase of Pi(m)', color='orange')
plt.xlabel('m')
plt.ylabel('Phase (radians)')
plt.grid(True)

plt.show()
