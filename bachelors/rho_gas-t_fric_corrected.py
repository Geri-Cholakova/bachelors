# -*- coding: utf-8 -*-
"""
Created on Fri May 16 11:18:20 2025

@author: Geri"""
import numpy as np
import math
from astropy import constants as const
from astropy import units as u
def velocities(s, x, y, z=0):
    x = (x + 1) *1000 * const.au
    y = y *1000 * const.au
    z = z *1000 * const.au #here in AU's 
    s *= u.m
    n = 1.5 #for P proport. to r**-n
    rho_s = 1.5* u.g/u.cm**3 
    r = np.sqrt(x**2 + y**2 + z**2)
    T = 200*u.K * (const.au / r)**0.5
    Sigma = 200*u.g/u.cm**2 * (const.au / r)
    mu = 2.3
    c_s = (const.k_B.to(u.kg*u.m**2/u.s**2/u.K) * T / (mu * const.m_p))**0.5
    Omega = (const.G * const.M_sun / r**3)**0.5
    h = c_s/Omega
    rho = (Sigma/(np.sqrt(2 * np.pi) * h))* np.e**((-z**2/(2*h)).value)
    rho2 = rho.to(u.kg/u.m**3)
    t_fric = rho_s / rho2 * s / (c_s * np.sqrt(8/np.pi))
    t_fric2 = t_fric.to(u.s)
    #print(rho2, t_fric2)
    St = t_fric2 * Omega
    #we assume homogenous gas disk with P prop. to r**-n
    nu = -0.5 * (h/r)**2 * n
    eps = (rho/ rho_s).si
    v_dust_r = float(-2/ (St + (1 + eps)**2/St) * nu * Omega*r *u.s/u.m)
    v_dust_phi = float(-(1+eps)/ ((1+eps)**2 +St**2) *nu * Omega*r *u.s/u.m)
    v_gas_r = float(-2*eps/(St + (1+eps)**2/St)*nu* Omega*r *u.s/u.m)
    v_gas_phi = float((eps/ ((1+eps)*(1+St**2/(1+eps)**2))-1) *nu * Omega*r *u.s/u.m)
    v_dust_x = v_dust_r*math.sin(v_dust_phi)
    v_dust_y = v_dust_r*math.cos(v_dust_phi)
    v_gas_x = v_gas_r*math.sin(v_gas_phi)
    v_gas_y = v_gas_r*math.cos(v_gas_phi)
    return v_dust_x, v_dust_y, v_gas_x, v_gas_y
vel = velocities(s=1, x=20, y=30)
print(vel)
