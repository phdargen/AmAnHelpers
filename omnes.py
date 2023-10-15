import math
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

# masses from the paper
m_pi = 0.13957
m_rho = 0.7736
Gamma_rho = 0.146
m_K = 0.496

# fit parameters (constrained)
PWavePars = namedtuple("PWavePars", ["B0", "B1", "l1", "l2", "epsilon1", "epsilon2", "e0"])
FWavePars = namedtuple("FWavePars", ["B0", "B1", "l"])

p_wave_pars = PWavePars(B0=1.043, B1=0.19, l1=1.39, l2=-1.70, epsilon1=0.00, epsilon2=0.07, e0=1.05)
f_wave_pars = FWavePars(B0=1.09e5, B1=1.41e5, l=0.051e5)

def k(s):
    return math.sqrt(s / 4 - m_pi**2)

def conformal_w(s, s0=1.45**2):
    return (math.sqrt(s) - math.sqrt(s0 - s)) / (math.sqrt(s) + math.sqrt(s0 - s))

# P-wave
def cot_delta1_less_1050(s, pars):
    if s < 4 * m_pi**2:
        return 0.0
    e0, B0, B1 = pars.e0, pars.B0, pars.B1
    return math.sqrt(s) / (2 * k(s)**3) * (m_rho**2 - s) * (2 * m_pi**3 / (m_rho**2 * math.sqrt(s)) + B0 + B1 * conformal_w(s, s0=e0**2))

def delta1_more_1050(s, pars):
    if s < 4 * m_pi**2:
        return 0.0
    l1, l2, epsilon1, epsilon2 = pars.l1, pars.l2, pars.epsilon1, pars.epsilon2
    l0 = np.pi/2 - np.arctan(cot_delta1_less_1050(4 * m_K**2, pars))
    return l0 + l1 * (math.sqrt(s) / (2 * m_K) - 1) + l2 * (math.sqrt(s) / (2 * m_K) - 1)**2

def _delta1(s, pars=p_wave_pars):
    if s < 4 * m_pi**2:
        return 0.0
    e0 = pars.e0
    if s < e0**2:
        v = np.pi/2 - np.arctan(cot_delta1_less_1050(s, pars))
        return v if v >= 0 else v + math.pi
    return delta1_more_1050(s, pars) if s < 1.4**2 else delta1_more_1050(1.4**2, pars)

def delta1(s, pars=p_wave_pars):
    v = _delta1(s, pars=pars)
    return v if s < 0.8**2 else (v if v > 0 else v + math.pi)

def cot_delta1(s, pars=p_wave_pars):
    return 1 / math.tan(delta1(s, pars))

# D-wave
def cot_delta3(s, pars=f_wave_pars):
    B0, B1, l_ = pars.B0, pars.B1, pars.l
    return math.sqrt(s) / (2 * k(s)**7) * m_pi**6 * (2 * l_ * m_pi / math.sqrt(s) + B0 + B1 * conformal_w(s))

def delta3(s, pars=f_wave_pars):
    return np.pi/2 - np.arctan(cot_delta3(s, pars))

# Plotting delta1(s)
s_values = np.linspace(0, 2, 500)
delta1_values = [delta1(s) for s in s_values]

plt.plot(s_values, delta1_values, label='delta1(s)')
plt.title('delta1(s) vs s')
plt.xlabel('s')
plt.ylabel('delta1(s)')
plt.grid(True)
plt.legend()
plt.show()
