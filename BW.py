import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from math import sqrt

# Fixed Parameters 
mpi = 0.13957
mK = 0.4937
m_B = 5.27934
m_psi2S = 3.6861

r = 1.5
L = 0
m1 = m_psi2S
m2 = mK

@np.vectorize
def kFactor(m,g):
  gamma = m * np.sqrt( m * m + g * g )
  k = (2 * np.sqrt(2) / np.pi) * m * g * gamma / sqrt( m*m + gamma )
  return np.sqrt(k)

def Q2(s, M1sq, M2sq):
  num = s -  2*M1sq - 2*M2sq  + ( M1sq - M2sq ) * ( M1sq - M2sq ) / s 
  return num/4.

def BlattWeisskopf_Norm(z,z0,L):
    return np.where(L==1,(1+z0) / (1+z),1)

def width(s, s1, s2, m, g, r, L, useEffMass=True):
  
  m0 = np.where(useEffMass and m <= (m1 + m2), m0_eff(m,m_psi2S+mK, m_B-mpi), m)

  q2v = Q2(s,s1,s2) 
  q2  = np.where(q2v > 0,q2v,0)
  #q2  = np.abs(q2v)
  q20 = np.abs( Q2( m0 * m0, s1, s2 ) )
  BF  = BlattWeisskopf_Norm( q2 * r * r, q20 * r * r, L )
  q2r = q2 / q20

  rt = g * BF * m * np.sqrt( q2r / s ) * q2r**L 

  #rt += 0.1 * np.sqrt(s)/m
  #rt = g * np.sqrt(s)/m

  return rt

def m0_eff(m,min,max):
   return min + 0.5 * (max-min) * (1 + np.tanh( (m-0.5*(min+max)) / (max-min)  ))

def BW(s, m, g, useEffMass=True):

    m0 = np.where(useEffMass and m <= (m1 + m2), m0_eff(m,m_psi2S+mK, m_B-mpi), m)

    q2 = np.abs( Q2( s, m1**2, m2**2 ) )  
    #q2 = np.where( Q2( s, m1**2, m2**2 ) > 0, Q2( s, m1**2, m2**2 ), 0)  

    q20 = np.abs( Q2( m0 * m0, m1**2, m2**2 ) ) 
    FormFactor = np.sqrt( BlattWeisskopf_Norm( q2 * r * r, 0, L ) )
    runningWidth = width( s, m1**2, m2**2, m, g, r, L, useEffMass )

    #runningWidth = np.where(np.sqrt(s) < (m1 + m2), g / m * np.sqrt(s), runningWidth)
    #runningWidth = g / m * np.sqrt(s)
    #runningWidth = g * runningWidth / width( m0**2, m1**2, m2**2, m, g, r, L, useEffMass )
    # m = m + (1j)*(1.e-6) 
    BW = FormFactor / ( m * m - s  -1j * m * runningWidth )
    kf = kFactor(m,g)

    return kf*BW

# Plot lineshapes
m_Zs = 4.00
g_Zs = 0.150

smin = (m_psi2S+mK)**2
#smin = (3.5)**2
smax = (m_B-mpi)**2
s_values = np.linspace(smin, smax, 1000)

m_Zs_min = 3.9
m_Zs_max = 6.5
m_Zs_values = np.linspace(m_Zs_min, m_Zs_max, 1000)

## BW
BW_values = np.abs(BW(smin*0.99, m_Zs_values, g_Zs, True))**2
BW_values2 = np.abs(BW(smin*1.0, m_Zs_values, g_Zs, True))**2
BW_values3 = np.abs(BW(smin*1.01, m_Zs_values, g_Zs, True))**2
BW_values4 = np.abs(BW(smin+(smax-smin)/2., m_Zs_values, g_Zs, True))**2

plt.figure(figsize=(10,6))
plt.plot(m_Zs_values, BW_values, label="BW")
plt.plot(m_Zs_values, BW_values2, label="BW2")
plt.plot(m_Zs_values, BW_values3, label="BW3")
plt.plot(m_Zs_values, BW_values4, label="BW4", linestyle="--")

plt.xlabel("$m_{Zs}$ (GeV)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True)
plt.show()

##
Q2_values = Q2(s_values, m1**2, m2**2)
plt.figure(figsize=(10,6))
plt.plot(np.sqrt(s_values), Q2_values, label="BW")
plt.xlabel("$\sqrt{s}$ (GeV)")
plt.ylabel("Q2")
plt.legend()
plt.grid(True)
#plt.show()

width_values = width(smax, m1**2, m2**2, m_Zs_values, g_Zs, r, L)
width_values2 = width(smax, m1**2, m2**2, m_Zs_values, g_Zs, r, L, False)
plt.figure(figsize=(10,6))
plt.plot(m_Zs_values, width_values, label="BW")
plt.plot(m_Zs_values, width_values2, label="BW2", linestyle="--")
plt.xlabel("$m_Zs$ (GeV)")
plt.ylabel("Width")
plt.legend()
plt.grid(True)
#plt.show()

width_values = width(s_values, m1**2, m2**2, m_Zs, g_Zs, r, L)
width_values2 = width(s_values, m1**2, m2**2, m_Zs, g_Zs, r, L, False)
plt.figure(figsize=(10,6))
plt.plot(s_values, width_values, label="BW")
plt.plot(s_values, width_values2, label="BW2", linestyle="--")
plt.xlabel("$s_values$ (GeV)")
plt.ylabel("Width")
plt.legend()
plt.grid(True)
plt.show()

m0_eff_values = m0_eff(np.sqrt(s_values), m_psi2S+mK, m_B-mpi )
plt.figure(figsize=(10,6))
plt.plot(np.sqrt(s_values), m0_eff_values, label="BW")
plt.xlabel("$\sqrt{s}$ (GeV)")
plt.ylabel("m0_eff_values")
plt.legend()
plt.grid(True)
#plt.show()