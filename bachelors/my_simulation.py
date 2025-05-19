# -*- coding: utf-8 -*-
"""
Created on Sun May 18 12:34:14 2025

@author: Geri
"""

import rebound
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
def my_simulation(r_pebble, a_pebble, f_pebble ): #in meters / AU, y=0.025 does the work
    sim = rebound.Simulation()
    sim.collision = "line"
    sim.units = ('yr', 'AU', 'Msun')
    sim.add(m=1, r=4.67e-3)
    sim.add(m=3e-6,r=4.26e-5, a=1)
    sim.add(m=4*np.pi*r_pebble**3*1500/(3*const.M_sun.value) ,r=r_pebble/const.au.value, a = a_pebble, f = f_pebble) # see rebound.orbit.Orbit?
    sim.move_to_com()
    nt = 4000
    data = np.zeros([nt, 2 * sim.N])
    time = np.linspace(0, 80, nt)
    
    data[0, :] = [i for p in sim.particles for i in [p.x, p.y]]
    try: 
        for i, t in enumerate(time[1:]):
            sim.integrate(t) 
            data[i + 1, :] = [coord for p in sim.particles for coord in [p.x, p.y]] 
        #using coord instead of i for clarity
    except rebound.Collision as error:
        print(error)
        data = data[:i+1, :]
    return data
data = my_simulation(r_pebble = 1, a_pebble = 0.993, f_pebble = 0.04)

width = 1.1
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

width = 1e-1
f, ax = plt.subplots()
ax.plot(x0, y0, 'k+-')
ax.plot(x1, y1, '+-')
ax.plot(x2, y2, '+--', linewidth=0.5)

ax.set_aspect(10)
ax.set(xlim=[1-width, 1+width], ylim=[-width/10, width/10])

