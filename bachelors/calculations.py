# -*- coding: utf-8 -*-
"""
Created on Fri May 16 11:18:20 2025

@author: Geri"""
import rebound
import numpy as np
import math
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt

norm = float(np.sqrt(const.G.value * const.M_sun.value / const.au.value))
r_H = 0.01

def cs_rho(x=1.0, y=0, z=0):
    x = x * const.au.value
    y = y * const.au.value
    z = z * const.au.value  # here in AU's
    r = np.sqrt(x ** 2 + y ** 2)
    T = 200 * np.sqrt(const.au.value / r)
    Sigma = 2000 * (const.au.value / r)
    mu = 2.3
    c_s = np.sqrt(const.k_B.value * T / (mu * const.m_p.value))
    Omega = np.sqrt(const.G.value * const.M_sun.value / r ** 3)
    h = c_s / Omega
    exponent = -(z ** 2) / (2 * h**2)
    rho = (Sigma / (np.sqrt(2 * np.pi) * h)) * np.e ** exponent
    return (c_s, rho)


def velocities_rh(s, x, y, r_h=r_H, n=-2.75, z=0, eps=0.01):
    # everything should be in si
    x = (x * r_h + 1) * const.au.value
    y = y * r_h * const.au.value
    z = z * r_h * const.au.value  # here in AU's
    rho_s = 1500
    R = np.sqrt(x ** 2 + y ** 2)
    T = 200 * np.sqrt(const.au.value / R)
    Sigma = 2000 * (const.au.value / R)
    mu = 2.3
    c_s = np.sqrt(const.k_B.value * T / (mu * const.m_p.value))
    Omega = np.sqrt(const.G.value * const.M_sun.value / R ** 3)
    h = c_s / Omega
    exponent = -(z ** 2) / (2 * h**2)
    rho = (Sigma / (np.sqrt(2 * np.pi) * h)) * np.e ** exponent
    t_fric = rho_s * s / (rho * c_s * np.sqrt(8 / np.pi))
    St = t_fric * Omega
    # we assume homogenous gas disk with P prop. to r**n
    # Omega prop r**-1.5, c_s prop r**-0.25, Sigma prop r**-1, P prop Omega*Sigma*c_s
    # So n = dlnP/dlnr = -2.75
    nu = -0.5 * (h / R) ** 2 * n
    # we use units for m=1, g=1, so we need to norm the result
    v_k = Omega * R
    v_dust_r = -2 / (St + (1 + eps) ** 2 / St) * nu * v_k
    v_dust_phi = v_k - (1 + eps) / ((1 + eps) ** 2 + St ** 2) * nu * v_k
    v_gas_r = -2 * eps / (St + (1 + eps) ** 2 / St) * nu * v_k
    v_gas_phi = (
        v_k + (eps / ((1 + eps) * (1 + St ** 2 / (1 + eps) ** 2)) - 1) * nu * v_k
    )

    sina = y / R
    cosa = x / R

    v_kep_x = -v_k * sina
    v_kep_y = v_k * cosa
    v_dust_x = v_dust_r * cosa - v_dust_phi * sina
    v_dust_y = v_dust_r * sina + v_dust_phi * cosa
    v_gas_x = v_gas_r * cosa - v_gas_phi * sina
    v_gas_y = v_gas_r * sina + v_gas_phi * cosa
    return (
        v_kep_x / norm,
        v_kep_y / norm,
        v_dust_x / norm,
        v_dust_y / norm,
        v_gas_x / norm,
        v_gas_y / norm,
    )

def velocities_cart(s, x, y, n=-2.75, z=0, eps=0.01, rho_s = 1500):
    # coord should be in rebound units
    x = x * const.au.value
    y = y * const.au.value
    z = z * const.au.value  # here in AU's
    R = np.sqrt(x ** 2 + y ** 2) #answer in si
    T = 200 * np.sqrt(const.au.value / R)
    Sigma = 2000 * (const.au.value / R)
    mu = 2.3
    c_s = np.sqrt(const.k_B.value * T / (mu * const.m_p.value))
    Omega = np.sqrt(const.G.value * const.M_sun.value / R ** 3)
    h = c_s / Omega
    exponent = -(z ** 2) / (2 * h**2)
    rho = (Sigma / (np.sqrt(2 * np.pi) * h)) * np.e**exponent
    t_fric = rho_s * s / (rho * c_s) * np.sqrt( np.pi/ 8)
    St = t_fric * Omega
    # we assume homogenous gas disk with P prop. to r**n
    # Omega prop r**-1.5, c_s prop r**-0.25, Sigma prop r**-1, P prop Omega*Sigma*c_s
    # So n = dlnP/dlnr = -2.75
    nu = -0.5 * (h / R) ** 2 * n
    # we use units for m=1, g=1, so we need to norm the result
    v_k = Omega * R
    
    A = (1 + eps)/((1+eps)**2 + St**2)
    B = -2/(St + (1+eps)**2/ St)
    
    v_dust_r = B * nu * v_k
    v_dust_phi = v_k - A * nu * v_k
    v_gas_r = -eps * B * nu * v_k
    v_gas_phi = v_k + (eps * A - 1) * nu * v_k

    sina = y / R
    cosa = x / R

    v_kep_x = -v_k * sina
    v_kep_y = v_k * cosa
    v_dust_x = v_dust_r * cosa - v_dust_phi * sina
    v_dust_y = v_dust_r * sina + v_dust_phi * cosa
    v_gas_x = v_gas_r * cosa - v_gas_phi * sina
    v_gas_y = v_gas_r * sina + v_gas_phi * cosa
    return (
        St,
        nu,
        v_dust_x / norm,
        v_dust_y / norm,
        v_gas_x / norm,
        v_gas_y / norm,
    )
"""
x2 = 1
y2 = 0

vel1 = velocities_rh(1e-6, x2, y2)          
vel2 = velocities_cart(s = 1e-6, x=x2, y = y2)  
cs, rho = cs_rho(x = 2, y = 0.,z=0.)
print(vel2[0], cs, rho)

sim = rebound.Simulation()
sim.add(m=1)
sim.add(m=2.09e-30, x=x1, y=y1, vx=vel[2], vy=vel[3])
sim.add(m=2.09e-30, x=x2, y=y2, vx=vel2[4], vy=vel2[5])
sim.move_to_com()
# to not have a drift because of center of mass movement
sim.integrate(2.0 * np.pi * 1)
# seeing the simulation after time t=100 for 2*pi*t=1 yr for a=1, M_central=1
rebound.OrbitPlot(sim)
# plot the orbital data after t iterations
print(f' after:{sim.particles[1].a}, {sim.particles[2].a}', vel[0], vel2[0])
print(sim.particles[1].e, sim.particles[2].e, vel2[6], vel2[7])

"""