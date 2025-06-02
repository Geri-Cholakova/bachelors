# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 11:45:20 2025

@author: Geri
"""
import rebound
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
import calculations as calc

M_star = 1
M_planet = 3e-6
r_pebble = 1
dist_pl = 1
r_H = dist_pl * np.cbrt(M_planet / (3 * (M_planet + M_star)))

Re = 0.5

sim = rebound.Simulation()
sim.add(m=M_star, r=4.67e-3)
sim.add(m=M_planet, r=4.26e-5, a=dist_pl)
sim.add(
    m=4 * np.pi * r_pebble ** 3 * 1500 / (3 * const.M_sun.value),
    r=r_pebble / const.au.value,
    a=dist_pl * 1.05,
)
sim.move_to_com()
ps = sim.particles
px = ps[2].x
py = ps[2].y
pz = ps[2].z

T_rho = calc.T_rho(x=(px - 1) / r_H, y=py / r_H, z=pz / r_H)

# defining a non-conservative force F = m*v/tau


def Stokesdrag(reb_sim):
    if Re < 1:
        pwr = -1
        A = 24
    elif Re < 800:
        pwr = -0.6
        A = 24
    else:
        pwr = 0
        A = 0.44
    C_D = A * Re ** pwr
    ps[2].ax -= (
        C_D
        * np.pi
        * (r_pebble / const.au.value) ** 2
        * T_rho[1]
        * ps[2].vx
        * np.sqrt(8 * const.k_B.value * T_rho[0] / (np.pi * const.m_p.value * 2.3))
    )
    ps[2].ay -= (
        C_D
        * np.pi
        * (r_pebble / const.au.value) ** 2
        * T_rho[1]
        * ps[2].vy
        * np.sqrt(8 * const.k_B.value * T_rho[0] / (np.pi * const.m_p.value * 2.3))
    )
    ps[2].az -= (
        C_D
        * np.pi
        * (r_pebble / const.au.value) ** 2
        * T_rho[1]
        * ps[2].vz
        * np.sqrt(8 * const.k_B.value * T_rho[0] / (np.pi * const.m_p.value * 2.3))
    )


sim.additional_forces = Stokesdrag
sim.force_is_velocity_dependent = 1

Nout = 1000
a_s = np.zeros(Nout)
times = np.linspace(0.0, 200.0 * 2.0 * np.pi, Nout)
for i, time in enumerate(times):
    sim.integrate(time, exact_finish_time=0)
    a_s[i] = sim.particles[2].a
rebound.OrbitPlot(sim)
fig = plt.figure(figsize=(15, 5))
ax = plt.subplot(111)
plt.plot(times, a_s)
