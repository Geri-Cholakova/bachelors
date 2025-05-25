# -*- coding: utf-8 -*-
"""
Created on Fri May 16 11:18:20 2025

@author: Geri"""
import numpy as np
from astropy import constants as const
from astropy import units as u
#
x = 1
y = 0
z = 0 #here in AU's 
s = 0.000001 #here in m
r = np.sqrt(x**2 + y**2)* const.au
z = z * const.au
T = 200*u.K * (const.au / r)**0.5
Sigma = 200*u.g/u.cm**2 * (const.au / r)
mu = 2.3
c_s = (const.k_B.to(u.kg*u.m**2/u.s**2/u.K) * T / (mu * const.m_p))**0.5
Omega = (const.G * const.M_sun / r**3)**0.5
h = c_s/Omega
rho = (Sigma/(np.sqrt(2 * np.pi) * h))* np.e**((-z**2/(2*h)).value)
rho2 = rho.to(u.kg/u.m**3)
t_fric = 1.5 * u.g/u.cm**3 / rho2 * s*u.m / (c_s * np.sqrt(8/np.pi))
t_fric2 = t_fric.to(u.s)
print(rho2, t_fric2)
