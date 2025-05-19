# -*- coding: utf-8 -*-
"""
Created on Sun May 18 15:56:02 2025

@author: Geri
"""

import rebound
import matplotlib.pyplot as plt
import numpy as np
import my_simulation as ms

all_data = [] # a list
initial_data = ms.my_simulation(r_pebble=100, a_pebble=0.99, f_pebble=0.03)
all_data.append(initial_data)
a_upper_limit = 1.01
a_lower_limit = 0.99
da = 0.005
    
phi = np.arctan2(initial_data[:, 3], initial_data[:, 2])

def rotate(x, y, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    return c * x - s * y, s * x + c * y
#rotation matrix [cosa -sina] [sina cosa]
#the fkt returns 2 elements 

width = 1e-1
x0, y0 = rotate(initial_data[:, 0], initial_data[:, 1], -phi)
x1, y1 = rotate(initial_data[:, 2], initial_data[:, 3], -phi) 
f, ax = plt.subplots()
ax.plot(x0, y0, 'k+-')
ax.plot(x1, y1, '+-')

for x in np.linspace(a_upper_limit, a_lower_limit, int((a_upper_limit - a_lower_limit)/da)):
    data = ms.my_simulation(r_pebble = 100, a_pebble = x, f_pebble = 0.03)
    xp, yp = rotate(data[:, 4], data[:, 5], -phi)
    ax.plot(xp, yp, '+--', linewidth= 0.5)

ax.set_aspect(10)
ax.set(xlim=[1-width, 1+width], ylim=[-width/10, width/10])