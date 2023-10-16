# Implements pipi S-wave from https://arxiv.org/pdf/1102.2183.pdf
import math
import numpy as np
import matplotlib.pyplot as plt

# Masses
m_pi = 0.13957
m_K = 0.496
m_eta = 0.54751

# Constants from paper
s_M = 0.85**2
s_0 = 4 * m_K**2
z_0 = m_pi
s_l = 1.05**2
s_h = 1.42**2

# CFD parameters from paper
B_0 = 7.14
B_1 = -25.3
B_2 = -33.2
B_3 = -26.2

d_0 = 226.5
c = -81
B = 93.3
C = 48.7
D = -88.3

B_0_S2 = -79.4
B_1_S2 = -63
z2 = 0.1435
B_h2 = 32

# Helper functions
def k_abs(s,m):
    return math.sqrt(abs(s / 4 - m**2))
def omega(s,s0):
    return (math.sqrt(s) - math.sqrt(s0 - s)) / (math.sqrt(s) + math.sqrt(s0 - s))

## S0 wave
def cotDelta1(s):
    w = omega(s,s_0)
    k = k_abs(s,m_pi)
    return math.sqrt(s)/(2*k) * m_pi**2/(s-0.5*z_0**2) * ( z_0**2/(m_pi*math.sqrt(s)) + B_0 + B_1 * w + B_2 * w**2 + B_3 * w**3 )

def delta1(s):
    return  180./np.pi * ( np.pi/2 - np.arctan(cotDelta1(s)) )

diff = 0.0001**2 # 0.1 MeV
d_M = delta1(s_M)
d_M_p = (delta1(s_M+diff) - delta1(s_M-diff)) / (2*diff)

def delta2(s):
    k2 = k_abs(s,m_K)
    kM = k_abs(s_M,m_K)
    return d_0 * (1 - k2/kM)**2 + d_M * k2/kM * (2-k2/kM) + k2 * (kM-k2) * (8 * d_M_p + c * (kM - k2)/m_K**3 )

def delta3(s):
    k2 = k_abs(s,m_K)
    k3 = k_abs(s,m_eta)
    val = d_0 + B * k2**2/m_K**2 + C * k2**4/m_K**4
    if s > 4 * m_eta**2: 
        val += D * k3**2/m_eta**2
    return val

d_M3 = delta3(1.42**2)
d_M3_p = (delta3(1.42**2+diff) - delta3(1.42**2-diff)) / (2*diff)

def delta4(s):
    a = 360
    alpha = d_M3_p/(a-d_M3)
    p = d_M3_p/(a*alpha)*math.exp(alpha*1.42**2)
    return a*(1.-p*math.exp(-alpha*s))

def delta(s):
    if s <= 4 * m_pi**2:
        return 0.0
    if s <= s_M:
        return delta1(s)
    if s <= 4*m_K**2:
        return delta2(s)
    if s <= 1.42**2:
        return delta3(s)
    #return delta3(1.42**2)   
    return delta4(s)    


## S2 wave
def cotDelta1_S2(s):
    w = omega(s,s_l)
    k = k_abs(s,m_pi)
    return math.sqrt(s)/(2*k) * m_pi**2/(s-2*z2**2) * ( B_0_S2 + B_1_S2 * w )

def delta1_S2(s):
    return  180./np.pi * ( np.pi/2 - np.arctan(cotDelta1_S2(s)) ) - 180.
 
def cotDelta2_S2(s):
    w = omega(s,s_h)
    w_S_M = omega(s_M,s_h)

    B_h0 = B_0_S2 + B_1_S2 * omega(s_M,s_l)
    B_h1 = B_1_S2 * s_l/s_h * math.sqrt(s_h - s_M) / math.sqrt(s_l - s_M) * ( (math.sqrt(s_M) + math.sqrt(s_h - s_M)) / (math.sqrt(s_M) + math.sqrt(s_l - s_M)) )**2

    k = k_abs(s,m_pi)
    return math.sqrt(s)/(2*k) * m_pi**2/(s-2*z2**2) * ( B_h0 + B_h1 * (w-w_S_M) + B_h2 * (w - w_S_M)**2 )

def delta2_S2(s):
    return  180./np.pi * ( np.pi/2 - np.arctan(cotDelta2_S2(s)) ) - 180.

def delta_S2(s):
    if s <= 4 * m_pi**2:
        return delta1_S2(4 * m_pi**2)
    if s <= s_M:
        return delta1_S2(s)
    if s <= 1.42**2:
        return delta2_S2(s)
    return delta2_S2(1.42**2)   
    #return delta4(s) 

# Plot delta(s)
m_values = np.linspace(2.001*m_pi, 2, 100)
#m_values = np.linspace(2.5*m_pi, 2, 100)
delta0_values = [delta(m**2) for m in m_values]
delta0_S2_values = [delta_S2(m**2) for m in m_values]

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(m_values, delta0_values, label='delta S0', color='red')
plt.plot(m_values, delta0_S2_values, label='delta S2', color='blue')
plt.xlabel('m')
plt.ylabel('delta(m)')
plt.grid(True)
plt.legend()
#plt.show()

## Omnes function from https://www.jlab.org/conferences/fdsa2014/talks/thursday/Dai.pdf
from scipy.integrate import quad

#ie = 1e-7j
sth = 4 * m_pi**2

def Integrand(s, zeta):
    return delta(zeta)*np.pi/180. / (zeta * (zeta - s))

def Integrand_S2(s, zeta):
    return delta_S2(zeta)*np.pi/180. / (zeta * (zeta - s))

def A_pipiS(s):
    integral_value1, err1 = quad(lambda zeta: Integrand(s, zeta), sth, s-diff, limit=500)  
    integral_value2, err2 = quad(lambda zeta: Integrand(s, zeta), s+diff, np.inf, limit=500)  

    integral_value =  integral_value1 + integral_value2
    # return integral_value
    return math.exp(integral_value * s / np.pi)

def A_pipiS2(s):
    integral_value1 = 0
    integral_value1, err1 = quad(lambda zeta: Integrand_S2(s, zeta), sth, s-diff, limit=500)  
    integral_value2, err2 = quad(lambda zeta: Integrand_S2(s, zeta), s+diff, np.inf, limit=500)  

    integral_value =  integral_value1 + integral_value2
    #return integral_value
    return math.exp(integral_value * s / np.pi)

# Plot A_pipiS
A_pipiS_values = [A_pipiS(m**2) for m in m_values]
A_pipiS2_values = [A_pipiS2(m**2) for m in m_values]

norm_A_pipiS = np.trapz(A_pipiS_values, m_values)
norm_A_pipiS2 = np.trapz(A_pipiS2_values, m_values)
#A_pipiS_values = A_pipiS_values /  norm_A_pipiS
#A_pipiS2_values = A_pipiS2_values /norm_A_pipiS2
A_pipiS_values = np.array(A_pipiS_values)
A_pipiS2_values = np.array(A_pipiS2_values)
A_pipiS_values = A_pipiS_values /  max(A_pipiS_values)
A_pipiS2_values = A_pipiS2_values / max(A_pipiS2_values) * 0.5

plt.subplot(1, 2, 2)
plt.plot(m_values, A_pipiS_values, label='A_pipiS0', color='red')
plt.plot(m_values, A_pipiS2_values, label='A_pipiS2', color='blue')
plt.xlabel('m')
plt.ylabel('A_pipiS(m)')
plt.grid(True)
plt.legend()
#ymax = 1.2 * max(A_pipiS_values)
#plt.ylim(0, ymax)  
plt.savefig('omnes_pipiSwave.png', dpi=300, bbox_inches='tight')
plt.show()
