# -*- coding: utf-8 -*-
"""
Created on Thu May  8 15:27:38 2025

@author: Geri
"""

import rebound
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const
if 'sim' in globals():
    globals()['sim'].stop_server()
    
#def collision_print_only(sim_pointer, collision):
    #sim = sim_pointer.contents           # get simulation object from pointer
    #print(sim.t)                         # print time 
    #print(sim.particles[collision.p1].x) # x position of particle 1
    #print(sim.particles[collision.p2].x) # x position of particle 2
    #return 0                             # Don't remove either particle

sim = rebound.Simulation()
sim.collision = "line"
#sim.collision_resolve = collision_print_only
#sim.exit_min_distance = 0.01
sim.units = ('yr', 'AU', 'Msun')
sim.add(m=1, r=4.67e-3)
sim.add(m=3e-6,r=4.26e-5, a=1)
r_pebble = 100 #in m
sim.add(m=4*np.pi*r_pebble**3*1500/(3*const.M_sun.value) ,r=r_pebble/const.au.value, a = 1.005, f = 0.04) # see rebound.orbit.Orbit?
sim.move_to_com()

p0 = sim.particles[0]
p1 = sim.particles[1]

#sim.widget()

nt = 4000
time = np.linspace(0, 80, nt) # returns evenly spaced numbers for (start, stop, number of samples)

data = np.zeros([nt, 2 * sim.N]) # returns array w zeroes with size nt rows and 2* sim.N columns 
# 2* sim.N because we need the x and y values of all particles
data[0, :] = [i for p in sim.particles for i in [p.x, p.y]]

#numpy arrays get def as np.array([1,2,3],[2,3,4]) for 2D
#you can print things like .size and .ndim
#you can do math like data*2 directly for all elements
#data[0:5] means slicing from index 0 to 4 incl
#data[0, :] means 0 row, all columns
#data[0:2, 2:10] means 0-1 incl elements frm row and 2-9 incl el from column

#the idea of the for for loop is that p takes values 0,1,2 
#and i takes values x, y 3 times and we just print 0.x, 0.y, 1.x, 1.y ... 
try: 
    for i, t in enumerate(time[1:]):
        sim.integrate(t) 
        data[i + 1, :] = [coord for p in sim.particles for coord in [p.x, p.y]] 
    #using coord instead of i for clarity
except rebound.Collision as error:
    print(error)
    data = data[:i+1, :]
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

ax.set_aspect(1)
ax.set(xlim=[1-width, 1+width], ylim=[-width, width])