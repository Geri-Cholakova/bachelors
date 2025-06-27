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

r_pebble = 1e-2 #here in m
dist_pl = 1
#distance is measured in AUs
M_star = 1
#mass is measured in solar units
inclin = 1e-7
eps = 0.01
rho_s = 1500

if inclin > 1:
    print("Your z is too large, there's a risk rho = 0")
else:
    
    norm_kg = 1/const.M_sun.value
    norm_m = 1/const.au.value #in rebound unit * m ^-1
    norm_s = 2.0 * np.pi/(1*u.year).si.value # in rebound unit * s^-1
    
    
    m_pebble= 4 * np.pi * r_pebble ** 3 * rho_s / 3 #in kg
    
    x = dist_pl * 1.0
    y = 0.0
    z = inclin * dist_pl
    
    reb_calc = calc.velocities_cart(r_pebble, x, y, z = z)
    St = reb_calc[0] #unitless
    v_dust_x = reb_calc[2] #rebound u
    v_dust_y = reb_calc[3] #rebound u
    R = np.sqrt(x**2 + y**2) #rebound u
    Re = 2
    
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.G = 1
    
    Omega = np.sqrt(sim.G * M_star/R**3) # in rebound u
    sim.add(m=M_star, r=4.67e-3)
    
    sim.add(m=m_pebble*norm_kg, r=r_pebble*norm_m, x=x, y=y, z=z, vx= v_dust_x, vy = v_dust_y, 
            vz = -St * z * Omega )
    sim.move_to_com()
    ps_pebble = sim.particles[1]
    # defining a non-conservative force F = m*v/tau
    # coeff for Stokes
    """
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
    coef_S = C_D * np.pi * (r_pebble) ** 2 * np.sqrt(8/np.pi) /norm_s / m_pebble
    # def stokes for x, y only
    
    def Stokesdrag(reb_sim):
        px = ps[1].x
        py = ps[1].y
        pz = ps[1].z
    
        c_s, rho = calc.cs_rho(x=px, y=py, z=pz)
        St, nu, v_dust_x, v_dust_y, v_gas_x, v_gas_y = calc.velocities_cart(
            s=r_pebble, x=px, y=py, z=pz
        )
    
        
        ps[1].ax -= coef_S * rho * (ps[1].vx - v_gas_x) * c_s 
        ps[1].ay -= coef_S * rho* (ps[1].vy - v_gas_y)* c_s
        ps[1].az -= coef_S * rho * ps[1].vz  * c_s """
    
    coef_E = np.sqrt(8 / np.pi) / (rho_s * r_pebble) /norm_s
    # we're norming that for rho_dust, 
    
    def Epsteindrag(reb_sim):
        px = ps_pebble.x
        py = ps_pebble.y
        pz = ps_pebble.z
    
    
        c_s, rho = calc.cs_rho(x=px, y=py, z=pz)
        St, nu, v_dust_x, v_dust_y, v_gas_x, v_gas_y = calc.velocities_cart(
            s = r_pebble, x = px, y= py, z= pz)
        # everything thats not in those parameters is incl in coef_E
    
        ps_pebble.ax -= coef_E * c_s *rho * (ps_pebble.vx - v_gas_x)
        ps_pebble.ay -= coef_E * c_s *rho * (ps_pebble.vy - v_gas_y) 
        ps_pebble.az -= coef_E * c_s *rho * ps_pebble.vz
        
        #i wanna break here
    
    sim.additional_forces = Epsteindrag
    sim.force_is_velocity_dependent = 1
    """
    
    #sim.integrate(sim.t + 2* np.pi )  # advance slightly so that a=/= 0 for first timestep
    
    Nout = 50000
    comparison = np.zeros((6, Nout))
    r = np.zeros(Nout)
    z = np.zeros(Nout)
    v_th = np.zeros(Nout)
    vz_th = np.zeros(Nout)
    v_code = np.zeros(Nout)    
    vz_code1 = np.zeros(Nout)
    v_code1 = np.zeros(Nout)
    vz_code = np.zeros(Nout)
    St = np.zeros(Nout)
    nu = np.zeros(Nout)
    z_th = np.zeros(Nout)
    v_gas_r = np.zeros(Nout)
    times = np.linspace(0.0, 10 * 2.0 * np.pi, Nout)
    for i, time in enumerate(times):
        sim.integrate(time, exact_finish_time=1 )
        if not np.isfinite(sim.particles[1].x):
           print(f"NaN encountered at t = {time}", ps[1].vx, ps[1].vy, ps[1].vz)
           break
        comparison[:, i] = calc.velocities_cart(s = r_pebble, x = ps[1].x, y = ps[1].y,
                                          z = ps[1].z)
        St[i] = comparison[0, i]
        nu[i] = comparison[1, i]  
        r[i] = np.hypot(ps[1].x, ps[1].y)#in rebound units
        z[i] = ps[1].z #in rebound units
        v_th[i] = -2/ (St[i] + (1+ eps)**2/St[i]) * nu[i] * np.sqrt(const.G.value * const.M_sun.value/(r[i]* const.au.value)) #in m/s
        vz_th[i] = -St[i] * z[i] * const.au.value * np.sqrt(const.G.value * const.M_sun.value/(r[i]* const.au.value)**3) #in m/s
        # approx only for St<1
        v_gas_r[i] = (comparison[4, i]* ps[1].x/r[i] + comparison[5, i]* ps[1].y/r[i])/norm_m #in m
        v_code1[i] = ps[1].vx* ps[1].x/r[i] + ps[1].vy* ps[1].y/r[i] #in rebound u
        vz_code1[i] = ps[1].vz #in rebound u
        
    #a_code = a_code * norm_s**2 / norm_s
    v_code = v_code1 * const.au.value *norm_s #in m/s
    vz_code = vz_code1 * const.au.value *norm_s #in m/s
    
    plt.figure(figsize=(15, 5))
    plt.plot(times/(2*np.pi), v_code, 'r-', label='v_code')
    plt.plot(times/(2*np.pi), v_th, 'g--', label='v_th')
    plt.legend(); plt.show()
    
    plt.figure(figsize=(15, 5))
    plt.plot(times/(2*np.pi), vz_code, 'r-', label='vz_code')
    plt.plot(times/(2*np.pi), vz_th, 'g--', label='vz_th')
    #plt.plot(times[:Nout-1]/(2*np.pi), k, 'b-', label='vz_calc_steps')
    plt.legend(); plt.show()

    fig = plt.figure(figsize=(15,5))
    ax = plt.subplot(111)
    plt.plot(times/(2*np.pi), r);
    
    fig = plt.figure(figsize=(15,5))
    ax = plt.subplot(111)
    plt.plot(times[1:]/(2*np.pi), (v_code[1:]-v_th[1:])/v_th[1:]);
    
    fig = plt.figure(figsize=(15,5))
    ax = plt.subplot(111)
    plt.plot(times/(2*np.pi), z * const.au.value);
    
    t_stop = St[0] / Omega
    del_v = np.mean(v_code - v_th)
    print("Delta x: ", del_v * t_stop)

   """ 
    