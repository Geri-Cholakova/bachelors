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
r_pebble = 100
dist_pl = 1
norm_kg = 1/const.M_sun.value
norm_m = 1/const.au.value


r_H = dist_pl * np.cbrt(M_planet / (3 * (M_planet + M_star)))
m_pebble_normed = 4 * np.pi * r_pebble ** 3 * 1500 / 3 * norm_kg

Re = 2

sim = rebound.Simulation()
sim.G = 1
sim.add(m=M_star, r=4.67e-3)
sim.add(m=M_planet, r=4.26e-5, a=dist_pl)
sim.add(
    m=m_pebble_normed, r=r_pebble *norm_m, a=dist_pl * 1, f=0.05,
)
sim.move_to_com()
ps = sim.particles
# defining a non-conservative force F = m*v/tau

# def stokes for x, y only
def Stokesdrag(reb_sim):
    px = ps[2].x
    py = ps[2].y
    pz = ps[2].z

    T, rho = calc.T_rho(x=(px - 1) / r_H, y=py / r_H, z=pz / r_H)
    v_kep_x, v_kep_y, v_dust_x, v_dust_y, v_gas_x, v_gas_y = calc.velocities(
        s=r_pebble, x=(px - 1) / r_H, y=py / r_H, z=pz / r_H
    )

    v_th = np.sqrt(8 * const.k_B.value * T / (np.pi * const.m_p.value * 2.3))

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
        * (r_pebble * norm_m) ** 2
        * rho * norm_kg/ norm_m**3
        * (ps[2].vx - v_gas_x)
        * v_th * norm_m
        / m_pebble_normed
    )
    ps[2].ay -= (
        C_D
        * np.pi
        * (r_pebble * norm_m) ** 2
        * rho * norm_kg/ norm_m**3
        * (ps[2].vy - v_gas_y)
        * v_th *norm_m
        / m_pebble_normed
    )
    """ps[2].az -= (
        C_D
        * np.pi
        * (r_pebble / const.au.value) ** 2
        * T_rho[1]
        * (vel[] - vel[]
        * v_th)/ m_pebble_normed
    )"""


# we've defined epstein for x, y only!!!!
def Epsteindrag(reb_sim):
    px = ps[2].x
    py = ps[2].y
    pz = ps[2].z

    T, rho = calc.T_rho(x=(px - 1) / r_H, y=py / r_H, z=pz / r_H)
    v_kep_x, v_kep_y, v_dust_x, v_dust_y, v_gas_x, v_gas_y = calc.velocities(
        s=r_pebble, x=(px - 1) / r_H, y=py / r_H, z=pz / r_H
    )

    v_th = np.sqrt(8 * const.k_B.value * T / (np.pi * const.m_p.value * 2.3))

    ps[2].ax -= -(
        4
        / 3
        * np.pi
        * rho *norm_kg/norm_m**3
        * (r_pebble / const.au.value) ** 2
        * v_th *norm_m
        * (ps[2].vx - v_gas_x) /m_pebble_normed )
    ps[2].ay -= -(
        4
        / 3
        * np.pi
        * rho *norm_kg/norm_m**3
        * (r_pebble / const.au.value) ** 2
        * v_th *norm_m
        * (ps[2].vy - v_gas_y) /m_pebble_normed )

sim.additional_forces = Epsteindrag
sim.force_is_velocity_dependent = 1

Nout = 5000
a_s = np.zeros(Nout)
data = np.zeros((Nout, 2))
times = np.linspace(0.0, 50 * 2.0 * np.pi, Nout)
for i, time in enumerate(times):
    sim.integrate(time, exact_finish_time=0)
    data[i] = sim.particles[2].vx, sim.particles[2].vy
    a_s[i] = sim.particles[2].a
rebound.OrbitPlot(sim)
fig = plt.figure(figsize=(15, 5))
ax = plt.subplot(111)
plt.plot(times, a_s)
lines = 135
