# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 11:45:20 2025

@author: Geri
"""
import rebound
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u
import calculations as calc

r_pebble = 0.7958
dist_pl = 1
M_star = 1
M_planet = 3e-6
inclin = 0

if inclin > 4e-6:
    print("Your z is too large, there's a risk rho = 0")
else:
    
    norm_kg = 1/const.M_sun.value
    norm_m = 1/const.au.value
    norm_s = 2.0 * np.pi/(1*u.year).si.value
    
    
    r_H = dist_pl * np.cbrt(M_planet / (3 * (M_planet + M_star)))
    m_pebble= 4 * np.pi * r_pebble ** 3 * 1500 / 3 
    
    Re = 2
    
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.G = 1
    sim.add(m=M_star, r=4.67e-3)
    #sim.add(m=M_planet, r=4.26e-5, a=dist_pl)
    sim.add(m=m_pebble*norm_kg, r=r_pebble*norm_m, a=dist_pl * 1.05, inc=inclin)
    sim.move_to_com()
    ps = sim.particles
    # defining a non-conservative force F = m*v/tau
    # coeff for Stokes
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
    
    # def stokes for x, y only
    
    def Stokesdrag(reb_sim):
        px = ps[1].x
        py = ps[1].y
        pz = ps[1].z
    
        c_s, rho = calc.T_rho(x=(px - 1) / r_H, y=py / r_H, z=pz / r_H)
        v_kep_x, v_kep_y, v_dust_x, v_dust_y, v_gas_x, v_gas_y = calc.velocities(
            s=r_pebble, x=(px - 1) / r_H, y=py / r_H, z=pz / r_H
        )
    
        
        ps[1].ax -= (
            C_D
            * np.pi
            * (r_pebble) ** 2
            * rho 
            * (ps[1].vx - v_gas_x)
            * c_s * np.sqrt(8/np.pi) /norm_s
            / m_pebble
        )
        ps[1].ay -= (
            C_D
            * np.pi
            * (r_pebble) ** 2
            * rho
            * (ps[1].vy - v_gas_y)
            * c_s * np.sqrt(8/np.pi) / norm_s
            / m_pebble
        )
        ps[1].az -= (
            C_D
            * np.pi
            * (r_pebble) ** 2
            * rho
            * ps[1].vz / norm_s
            * c_s * np.sqrt(8/np.pi)/ m_pebble
        )
    
    coef_E = np.sqrt(8 / np.pi)/ (1500 * r_pebble)/ norm_s
    # we're norming that for rho_dust, 
    # we've defined epstein for x, y only!!!!
    
    def Epsteindrag(reb_sim):
        px = ps[1].x
        py = ps[1].y
        pz = ps[1].z
    
    
        c_s, rho = calc.cs_rho(x=(px - 1) / r_H, y=py / r_H, z=pz / r_H)
        v_kep_x, v_kep_y, v_dust_x, v_dust_y, v_gas_x, v_gas_y = calc.velocities(
            s = r_pebble, x = (px - 1)/r_H, y= py/r_H, z= pz/r_H)
        # everything thats not in those patameters is incl in coef_E
    
        ps[1].ax -= coef_E * c_s *rho * (ps[1].vx - v_gas_x)
        ps[1].ay -= coef_E * c_s *rho * (ps[1].vy - v_gas_y) 
        ps[1].az -= coef_E * c_s *rho * ps[1].vz
        
        #i wanna break here
    
    sim.additional_forces = Epsteindrag
    sim.force_is_velocity_dependent = 0
    
    Nout = 5000
    r = np.zeros(Nout)
    times = np.linspace(0.0, 100 * 2.0 * np.pi, Nout)
    for i, time in enumerate(times):
        sim.integrate(time, exact_finish_time=0 )
        if not np.isfinite(sim.particles[1].x):
           print(f"NaN encountered at t = {time}", ps[1].vx, ps[1].vy, ps[1].vz)
           break
        r[i] = np.hypot(sim.particles[1].x, sim.particles[1].y)
    rebound.OrbitPlot(sim)
    fig = plt.figure(figsize=(15,5))
    ax = plt.subplot(111)
    plt.plot(times, r);