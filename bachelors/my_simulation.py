# -*- coding: utf-8 -*-
"""
Created on Sun May 18 12:34:14 2025

@author: Geri
"""
import rebound
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
import calculations as calc
from force_example import make_epstein_drag

class DistanceException(Exception):
    """Particle with index 2 was too far away from
    particle with index 1.
    """      
    def __str__(self):
        return("Particle gets really far away")
    pass

def distance_check(sim):
    "heartbeat function to stop simulation"
    #sim = sim.contents
    p0 = sim.particles[1]
    p1 = sim.particles[2]
    d = np.sqrt((p0.x - p1.x)**2 + (p0.y - p1.y)**2 + (p0.z - p1.z)**2)
    if d > 0.6:
        raise DistanceException
        
def my_simulation( r_pebble, x_h, y_h, force, z_h=0, nt=50000, n_orb=200, dist_pl=1, M_star=1, M_planet=3e-6, E_coef = 5343.21, check_distance = True ):# in meters / Hill spheres, here 0.01au
    collision = False
    sim = rebound.Simulation()
    sim.collision = "line"
    sim.units = ("yr", "AU", "Msun")
    sim.G = 1

    r_H = dist_pl * np.cbrt(M_planet / (3 * (M_planet + M_star)))
    St, eta, v_dust_x, v_dust_y, v_gas_x, v_gas_y = calc.velocities_rh(s=r_pebble, x=x_h, y=y_h, r_h=r_H, z=z_h)
    
    R = np.sqrt((x_h * r_H + 1)**2 + (y_h * r_H)**2)
    Omega = np.sqrt(sim.G * M_star/R**3)

    sim.add(m=M_star, r=4.67e-3)
    sim.add(m=M_planet, r=4.26e-5, a=dist_pl)
    sim.add(
        m=4 * np.pi * r_pebble ** 3 * 1500 / (3 * const.M_sun.value),
        r=r_pebble / const.au.value,
        x=x_h * r_H + 1,
        y=y_h * r_H,
        z=z_h * r_H,
        vx=v_dust_x,
        vy=v_dust_y, vz= -St * Omega * z_h * r_H 
    )  # see rebound.orbit.Orbit?
    # vel[0/1] - kepl velocity, vel[2/3] - dust velocity
    sim.move_to_com()


    data = np.zeros([nt + 1, 2 * sim.N])
    time = np.linspace(0, n_orb, nt)
    
    if force == True:
        sim.additional_forces = make_epstein_drag(ps_pebble = sim.particles[2], E_coef = E_coef) 
        sim.force_is_velocity_dependent = 1

    data[0, :] = [i / r_H for p in sim.particles for i in [p.x, p.y]]
    try:
        for i, t in enumerate(time[1:]):
            sim.integrate(t)
            if check_distance:
               distance_check(sim)
            data[i + 1, :] = [coord / r_H for p in sim.particles for coord in [p.x, p.y]]
        # using coord instead of i for clarity
    except (rebound.Collision, DistanceException) as e:
        if isinstance(e, rebound.Collision):
            collision = True
        for j in range(6):    
            data[i + 1, j] = float(collision)
        data = data[:i + 2, :]
        print(f"{e} for x_start = {x_h}, y_start = {y_h}")
    return data


def worker(arg, r_pebble, E_coef, n_orb = 20, nt = 50000, force = True):
    return my_simulation(r_pebble= r_pebble, x_h=arg[0], y_h=arg[1], z_h=0, n_orb=n_orb, nt=nt, force = force, E_coef = E_coef)


"""
data = my_simulation(r_pebble = 1, x_h = 1., y_h = 0.5, force = False) #IN HILL RADII FOR X,Y

width = 300
f, ax = plt.subplots()
ax.plot(data[:-1, 0], data[:-1, 1], 'k+-')
ax.plot(data[:-1, 2], data[:-1, 3], '+-')
ax.plot(data[:-1, 4], data[:-1, 5], '+--')
ax.set(xlim=[-width, width], ylim=[-width, width])
ax.set_aspect(1)
#ax.plot(x, y, style)


phi = np.arctan2(data[:-1, 3], data[:-1, 2])

def rotate(x, y, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    return c * x - s * y, s * x + c * y
#rotation matrix [cosa -sina] [sina cosa]
#the fkt returns 2 elements 

x0, y0 = rotate(data[:-1, 0], data[:-1, 1], -phi)
x1, y1 = rotate(data[:-1, 2], data[:-1, 3], -phi)
x2, y2 = rotate(data[:-1, 4], data[:-1, 5], -phi)

width = 1e2
f, ax = plt.subplots()
ax.plot(x0, y0, 'k+-')
ax.plot(x1, y1, '+-')
ax.plot(x2, y2, 'r--', linewidth=0.5)

ax.set_aspect(1)
ax.set(xlim=[100-width, 100+width], ylim=[-width, width]) """

