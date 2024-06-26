import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from math import sqrt

####
# Calculate K1 BFs from fit fractions obtained from amplitude fits
# Correct for isopin factors
# Compute omega contribution to rho/omega mixing lineshape
###

# Fixed Parameters 
mpi = 0.13957
mK = 0.4937
m_B = 5.27934
m_psi2S = 3.6861

m_rho = 0.77526  
gamma_rho = 0.1474  
m_omega = 0.78266 
gamma_omega = 0.00868  

radius = 1.2

reso = "K1"
m_K1 = 1.279
g_K1 = 0.116
if(reso == "K1p"):
  m_K1 = 1.403
  g_K1 = 0.174
if(reso == "K"):
  m_K1 = 1.4824
  g_K1 = 0.3356
if(reso == "K2*"):
  m_K1 = 1.4273
  g_K1 = 0.100
if(reso == "K*"):
  m_K1 = 1.414
  g_K1 = 0.232
if(reso == "K*p"):
  m_K1 = 1.718
  g_K1 = 0.322

# Rho/Omega parameters from amp fit
delta = 0.0036
sigma_delta = 0.0009
phi = 0.0 
if(reso != "K1"):
  delta = 0.0038 
  sigma_delta = 0.0013


# Implements Rho/Omega mixing from https://arxiv.org/pdf/hep-ex/0112031.pdf
def safe_sqrt(x):
  return np.where(x > 0, np.sqrt(np.abs(x)), 0)

def q_2(s):
  return s - 4*mpi*mpi

def q(s):
  return safe_sqrt(q_2(s))

def logTerm(s):
  return np.log( ( np.sqrt(s) + safe_sqrt(q_2(s)) ) / ( 2*mpi ) )

def kFactor(m,g):
  gamma = m * np.sqrt( m * m + g * g )
  k = (2 * np.sqrt(2) / np.pi) * m * g * gamma / sqrt( m*m + gamma )
  return np.sqrt(k)

def Q2(s, M1sq, M2sq):
  num = s -  2*M1sq - 2*M2sq  + ( M1sq - M2sq ) * ( M1sq - M2sq ) / s 
  return num/4.

def BlattWeisskopf_Norm(z,z0,L):
    return np.where(L==1,(1+z0) / (1+z),1)

def width(s, s1, s2, m, g, r, L):
  q2v = Q2(s,s1,s2) 
  q2  = np.where(q2v > 0,q2v,0)
  q20 = np.abs( Q2( m * m, s1, s2 ) )
  BF  = BlattWeisskopf_Norm( q2 * r * r, q20 * r * r, L )
  q2r = q2 / q20

  rt = g * BF * m * np.sqrt( q2r / s ) * q2r**L 
  return rt

def GS(s, m, g):
    s0 = m * m
    t1 = 2 * s0 * q_2(s) * q(s) *logTerm(s)/(np.sqrt(s)* q(s0)**3)
    t2 = logTerm(s0)*( q_2(s)*(s - 3*s0) + s*(s0-s) )/(m*q_2(s0)) 
    t3 = (s0-s)/q(s0)
    M2 = s0 + g * (t1 + t2 + t3 )/np.pi; 

    s1 = mpi**2
    s2 = mpi**2
    D = M2 - 1j*m*width(s, s1, s2, m, g, 0, 1)
    q2  = np.abs( Q2( s, s1, s2 ) )
    FormFactor = np.sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, 1 ) )
    GS = kFactor(m,g) * FormFactor / (D-s)

    return GS

def BW(s, m, g):
    return m**2 / (m**2 - s - 1j*m*g)

def RhoOmega(s, m_rho, gamma_rho, m_omega, gamma_omega, delta, phi=0, a=1):
    d = delta * ( np.cos(phi*np.pi/180.) + 1j * np.sin(phi*np.pi/180.)  )
    q2  = np.abs( Q2( s, mpi**2, mpi**2) )
    FormFactor = np.sqrt( q2 )
    return FormFactor * GS(s,m_rho,gamma_rho) * (a + d * s / m_omega**2 * BW(s, m_omega, gamma_omega) )

def integrand_GS(s, m, g):
    val = GS(s, m, g)
    return np.abs(val)**2

def integrand_omega(s,m,g, delta, phi=0):
       d = delta * ( np.cos(phi*np.pi/180.) + 1j * np.sin(phi*np.pi/180.)  )
       val = d * s / m**2 * BW(s, m, g)
       return  np.abs(val)**2

def integrand_BW(s, m, g):
    val = BW(s, m, g)
    return np.abs(val)**2

def integrand_RhoOmega_s12(s12, s123, m_rho, gamma_rho, m_omega, gamma_omega, delta, phi=0, a=1):
    q2  = np.abs( Q2( s123, s12, mK**2 ) )
    FormFactor = np.sqrt( q2 )
    val = FormFactor*RhoOmega(s12, m_rho, gamma_rho, m_omega, gamma_omega, delta, phi, a)
    return np.abs(val)**2

def integrand_RhoOmega_s12_s123(s123, m_rho, gamma_rho, m_omega, gamma_omega, delta, phi=0, a=1):
    integral_rhoOmega, _ = quad(integrand_RhoOmega_s12, (2*mpi)**2, (np.sqrt(s123)-mK)**2, args=(s123, m_rho, gamma_rho, m_omega, gamma_omega, delta, phi, a), limit=500)
    bw = BW(s123, m_K1, g_K1/m_K1*np.sqrt(s123)) #TODO::Implement K1 running width

    return integral_rhoOmega * np.abs(bw)**2

# Plot lineshapes
smin = (2*mpi)**2
smax = (m_K1-mK)**2
s_values = np.linspace(smin, smax, 1000)

rho_values = np.abs(RhoOmega(s_values, m_rho, gamma_rho, m_omega, gamma_omega, 0, phi))**2
omega_values = np.abs(RhoOmega(s_values, m_rho, gamma_rho, m_omega, gamma_omega, 1, phi, 0))**2
rhoOmega_values = np.abs(RhoOmega(s_values, m_rho, gamma_rho, m_omega, gamma_omega, delta, phi))**2

rho_values_norm = rho_values / np.max(rho_values)
omega_values_values_norm = omega_values / np.max(omega_values)
rhoOmega_values_norm = rhoOmega_values / np.max(rhoOmega_values)

plt.figure(figsize=(10,6))
plt.plot(np.sqrt(s_values), rho_values_norm, label="GS for rho")
plt.plot(np.sqrt(s_values), omega_values_values_norm, label="BW for omega", linestyle="--")
plt.plot(np.sqrt(s_values), rhoOmega_values_norm, label="RhoOmega", linestyle="--")
plt.xlabel("$\sqrt{s}$ (GeV)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True)
#plt.show()

# Integrate squared amplitudes for s123=m_K1**2
integral_rho, _ = quad(integrand_RhoOmega_s12, smin, smax, args=(m_K1**2, m_rho, gamma_rho,m_omega, gamma_omega, 0, phi))
integral_omega, _ = quad(integrand_RhoOmega_s12, smin, smax, args=(m_K1**2, m_rho, gamma_rho,m_omega, gamma_omega, delta, phi, 0))
integral_rhoOmega, _ = quad(integrand_RhoOmega_s12, smin, smax, args=(m_K1**2, m_rho, gamma_rho,m_omega, gamma_omega, delta, phi))

print("Integrate squared amplitudes for s123=m_K1**2")
print(f"I_tot= {integral_rhoOmega}")
print(f"I_rho= {integral_rho}")
print(f"I_omega= {integral_omega}")

R_all_fixed = (1-integral_rho/integral_rhoOmega) * 100
R_0_fixed = integral_omega/integral_rhoOmega * 100 
print(f"R^all_omega= {R_all_fixed}")
print(f"R^0_omega= {R_0_fixed}")

# Integrate squared amplitudes weighted by K1 BW
s123_min = (m_K1-3*g_K1)**2  #(mK+2*mpi)**2
s123_max = (m_B-m_psi2S)**2 #(m_K1+10*g_K1)**2
integral_rho, _ = quad(integrand_RhoOmega_s12_s123, s123_min, s123_max, args=(m_rho, gamma_rho,m_omega, gamma_omega, 0, phi), limit=500)
integral_omega, _ = quad(integrand_RhoOmega_s12_s123, s123_min, s123_max, args=(m_rho, gamma_rho,m_omega, gamma_omega, delta, phi, 0), limit=500)
integral_rhoOmega, _ = quad(integrand_RhoOmega_s12_s123, s123_min, s123_max, args=(m_rho, gamma_rho,m_omega, gamma_omega, delta, phi), limit=500)

print("Integrate squared amplitudes weighted by K1 BW")
print(f"I_tot= {integral_rhoOmega}")
print(f"I_rho= {integral_rho}")
print(f"I_omega= {integral_omega}")

R_all = (1-integral_rho/integral_rhoOmega) * 100
R_0 = integral_omega/integral_rhoOmega * 100 
print(f"R^all_omega= {R_all}")
print(f"R^0_omega= {R_0}")

## Shift delta by sigma(delta)
integral_omega, _ = quad(integrand_RhoOmega_s12_s123, s123_min, s123_max, args=(m_rho, gamma_rho,m_omega, gamma_omega, delta-sigma_delta, phi, 0), limit=500)
integral_rhoOmega, _ = quad(integrand_RhoOmega_s12_s123, s123_min, s123_max, args=(m_rho, gamma_rho,m_omega, gamma_omega, delta-sigma_delta, phi), limit=500)
print("Integrate squared amplitudes weighted by K1 BW: delta-sigma_delta")
R_all_m = (1-integral_rho/integral_rhoOmega) * 100
R_0_m = integral_omega/integral_rhoOmega * 100 
print(f"R^all_omega= {R_all_m}")
print(f"R^0_omega= {R_0_m}")

integral_omega, _ = quad(integrand_RhoOmega_s12_s123, s123_min, s123_max, args=(m_rho, gamma_rho,m_omega, gamma_omega, delta+sigma_delta, phi, 0), limit=500)
integral_rhoOmega, _ = quad(integrand_RhoOmega_s12_s123, s123_min, s123_max, args=(m_rho, gamma_rho,m_omega, gamma_omega, delta+sigma_delta, phi), limit=500)

print("Integrate squared amplitudes weighted by K1 BW: delta+sigma_delta")
R_all_p = (1-integral_rho/integral_rhoOmega) * 100
R_0_p = integral_omega/integral_rhoOmega * 100 
print(f"R^all_omega= {R_all_p}")
print(f"R^0_omega= {R_0_p}")

print("\n")
print(f"R^0_omega= {R_0} pm {np.abs(R_0_p-R_0_m)/2}")
print(f"R^0_rho= {100-R_0} pm {np.abs(R_0_p-R_0_m)/2}")
print("\n")

## Plot integrand to check valid integration range
s_values = np.linspace(s123_min, s123_max, 1000)
integrand_RhoOmega_s12_s123_values = [integrand_RhoOmega_s12_s123(s, m_rho, gamma_rho, m_omega, gamma_omega, delta, phi)**2 for s in s_values]
plt.figure(figsize=(10,6))
plt.plot(np.sqrt(s_values), integrand_RhoOmega_s12_s123_values)
plt.xlabel("$\sqrt{s}$ (GeV)")
plt.grid(True)
#plt.show()


## Calculate BFs

## Check isospin factors with results from https://arxiv.org/pdf/1009.5256v2.pdf
ffs = [0.383, 0.0157 ,0.232, 0.0045 ] # Belle fit 1
#ffs = [0.430, 0.0184 ,0.168, 0.00758]  # Belle fit 2
#ffs = [, 0.0184 ,0.168, 0.00758]  # Belle fit 3

total = sum(ffs)
normalized_ffs = [x/total * 100 for x in ffs]
print("Belle Fi: K rho, K pipi_S, K* pi, K omega")
print(normalized_ffs)

BF_K0s2Kpi = 0.93
BF_omega2pipi = 0.0153

bfs = [ffs[0], ffs[1] / BF_K0s2Kpi * 3./4., ffs[2] * 3./4., ffs[3] / BF_omega2pipi * 1./3. ]
total = sum(bfs)
normalized_bfs = [x/total * 100 for x in bfs]
print("Belle BFs: K rho, K pipi_S, K* pi, K omega")
print(normalized_bfs)

#Amp fit v1
# f_rho = 57.07
# f_rho_err = 2.44
# f_Ks_S = 13.70
# f_Ks_S_err = 1.57
# f_Ks_D = 5.89
# f_Ks_D_err = 0.72
# f_Kpi = 8.89
# f_Kpi_err = 1.29
#Sum_Ks_S_D = 97.44/100

# f_rho = 0.5516803 *100
# f_rho_err = 0.0091285881 *100
# f_Ks_S = 0.17489425 *100
# f_Ks_S_err = 0.01225596 *100
# f_Ks_D = 0.066988013 *100
# f_Ks_D_err = 0.0064358485 *100
# f_Kpi = 0.08634484 *100
# f_Kpi_err = 0.001427722 *100
#Sum_Ks_S_D = 97.44/100

#Amp fit v2
# K1(1270)
f_rho = 50.71
f_rho_err = sqrt(2.18**2 + 3.19**2)
f_Ks_S = 19.86
f_Ks_S_err = sqrt(1.44**2 + 2.05**2)
f_Ks_D = 8.32
f_Ks_D_err = sqrt(0.85**2 + 1.54**2)
f_Kpi = 11.35
f_Kpi_err = sqrt(1.45**2 + 2.11**2)
Sum_Ks_S_D = 98.64/100

# K1(1400)
if(reso == "K1p"):
  f_rho = 22.08
  f_rho_err = sqrt(4.4**2 + 6.25**2)
  f_Ks_S = 86.8
  f_Ks_S_err = sqrt(3.78**2 + 5.34**2)
  f_Ks_D = 0
  f_Ks_D_err = 0
  f_Kpi = 0
  f_Kpi_err = 0
  Sum_Ks_S_D = 1

# K(1460)
if(reso == "K"):
  f_rho = 0
  f_rho_err = 0
  f_Ks_S = 35.41
  f_Ks_S_err = sqrt(4.08**2 + 9.72**2)
  f_Ks_D = 0
  f_Ks_D_err = 0
  f_Kpi = 45.13
  f_Kpi_err = sqrt(4.22**2 + 10.91**2)
  Sum_Ks_S_D = 1

# K2*(1430)
if(reso == "K2*"):
  f_rho = 12.71
  f_rho_err = sqrt(2.3**2 + 1.8**2)
  f_Ks_S = 76.7
  f_Ks_S_err = sqrt(3.04**2 + 2.43**2)
  f_Ks_D = 0
  f_Ks_D_err = 0
  f_Kpi = 0
  f_Kpi_err = 0
  Sum_Ks_S_D = 1

# K*(1410)
if(reso == "K*"):
  f_rho = 38.36
  f_rho_err = sqrt(10.46**2 + 19.07**2)
  f_Ks_S = 88.5
  f_Ks_S_err = sqrt(8.39**2 + 12.65**2)
  f_Ks_D = 0
  f_Ks_D_err = 0
  f_Kpi = 0
  f_Kpi_err = 0
  Sum_Ks_S_D = 1

# K*(1680)
if(reso == "K*p"):
  f_rho = 31.16
  f_rho_err = sqrt(6.11**2 + 11.27**2)
  f_Ks_S = 49.69
  f_Ks_S_err = sqrt(6.72**2 + 13.32**2)
  f_Ks_D = 0
  f_Ks_D_err = 0
  f_Kpi = 0
  f_Kpi_err = 0
  Sum_Ks_S_D = 1

ffs = [ (1-R_0/100) * f_rho, f_Kpi, (f_Ks_S+f_Ks_D)/Sum_Ks_S_D, R_0/100 * f_rho ]  
ffs_err = [ sqrt( ((1-R_0/100)* f_rho_err)**2 + (np.abs(R_0_p-R_0_m)/2/100 *ffs[0])**2 ), f_Kpi_err, sqrt(f_Ks_S_err**2 + f_Ks_D_err**2)/Sum_Ks_S_D, sqrt( (np.abs(R_0_p-R_0_m)/2/100 * f_rho)**2 + (R_0/100 *f_rho_err)**2 ) ] 

bfs = [ffs[0], ffs[1] * 3./4., ffs[2] * 3./4., ffs[3] / BF_omega2pipi * 1./3.]
total = sum(bfs)
normalized_bfs = [x/total * 100 for x in bfs]

bfs_err = [ffs_err[0], ffs_err[1] * 3./4., ffs_err[2] * 3./4. , ffs_err[3] / BF_omega2pipi * 1./3.]
normalized_bfs_err = [x/total * 100 for x in bfs_err]

print("Baseline: K rho, K pipi_S, K* pi, K omega")
print(normalized_bfs)
print(normalized_bfs_err)