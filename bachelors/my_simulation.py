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


def my_simulation(
    r_pebble, x_h, y_h, z_h=0, nt=50000, n_orb=200, dist_pl=1, M_star=1, M_planet=3e-6
):  # in meters / Hill spheres, here 0.01au
    sim = rebound.Simulation()
    sim.collision = "line"
    sim.units = ("yr", "AU", "Msun")
    sim.G = 1

    r_H = dist_pl * np.cbrt(M_planet / (3 * (M_planet + M_star)))
    vel = calc.velocities_rh(s=r_pebble, x=x_h, y=y_h, r_h=r_H, z=z_h)

    sim.add(m=M_star, r=4.67e-3)
    sim.add(m=M_planet, r=4.26e-5, a=dist_pl)
    sim.add(
        m=4 * np.pi * r_pebble ** 3 * 1500 / (3 * const.M_sun.value),
        r=r_pebble / const.au.value,
        x=x_h * r_H + 1,
        y=y_h * r_H,
        z=z_h * r_H,
        vx=vel[2],
        vy=vel[3],
    )  # see rebound.orbit.Orbit?
    # vel[0/1] - kepl velocity, vel[2/3] - dust velocity
    sim.move_to_com()

    p0 = sim.particles[0]
    p1 = sim.particles[1]
    p2 = sim.particles[2]

    data = np.zeros([nt, 2 * sim.N])
    time = np.linspace(0, n_orb, nt)

    data[0, :] = [i / r_H for p in sim.particles for i in [p.x, p.y]]
    try:
        for i, t in enumerate(time[1:]):
            sim.integrate(t)
            data[i + 1, :] = [
                coord / r_H for p in sim.particles for coord in [p.x, p.y]
            ]
        # using coord instead of i for clarity
    except rebound.Collision as error:
        print(f"{error} for x_start = {x_h}, y_start = {y_h}")
        data = data[: i + 1, :]
    return data


def worker(arg):
    return my_simulation(r_pebble=1, x_h=arg[0], y_h=arg[1], z_h=0, n_orb=20, nt=50000)


"""
data = my_simulation(r_pebble = 1, x_h = 10, y_h = 40) #IN HILL RADII FOR X,Y

width = 1100
f, ax = plt.subplots()
ax.plot(data[:, 0], data[:, 1], 'k+-')
ax.plot(data[:, 2], data[:, 3], '+-')
ax.plot(data[:, 4], data[:, 5], '+--')
ax.set(xlim=[-width, width], ylim=[-width, width])
ax.set_aspect(1)
#ax.plot(x, y, style)


phi = np.arctan2(data[:, 3], data[:, 2])

def rotate(x, y, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    return c * x - s * y, s * x + c * y
#rotation matrix [cosa -sina] [sina cosa]
#the fkt returns 2 elements 

x0, y0 = rotate(data[:, 0], data[:, 1], -phi)
x1, y1 = rotate(data[:, 2], data[:, 3], -phi)
x2, y2 = rotate(data[:, 4], data[:, 5], -phi)

width = 1e2
f, ax = plt.subplots()
ax.plot(x0, y0, 'k+-')
ax.plot(x1, y1, '+-')
ax.plot(x2, y2, 'r--', linewidth=0.5)

ax.set_aspect(1)
ax.set(xlim=[1000-width, 1000+width], ylim=[-width, width]) 
"""
