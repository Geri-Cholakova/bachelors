# -*- coding: utf-8 -*-
"""
Created on Sun May 18 15:56:02 2025

@author: Geri
"""

import matplotlib.pyplot as plt
import numpy as np
import my_simulation as ms
from multiprocessing.pool import Pool
import time

start = time.time()

all_data = [] # a list
initial_data = ms.my_simulation(r_pebble=100, x_h=0, y_h=100)
#here the assumption is that initial data doesn't have a collision; that in diff simulations the difference in planet placement is negligible
all_data.append(initial_data)

x_upper = 40
x_lower = 20
dx = 2 #we'll use it for dy as well

y_upper = 40
y_lower = 20
dy = 2


phi = np.arctan2(initial_data[:, 3], initial_data[:, 2])

def rotate_h(x, y, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    return c * x - s * y - 999.985, s * x + c * y
#rotation matrix [cosa -sina] [sina cosa]
#the fkt returns 2 elements 

width = 50
x0, y0 = rotate_h(initial_data[:, 0], initial_data[:, 1], -phi)
x1, y1 = rotate_h(initial_data[:, 2], initial_data[:, 3], -phi) 
f, ax = plt.subplots()
ax.plot(x0 , y0, 'k+-')
ax.plot(x1 , y1, 'k--')
r_hill = plt.Circle((0, 0), 1, color='c')
r_earth = plt.Circle((0, 0), 4.267e-2, color='b')
ax.add_patch(r_hill)
ax.add_patch(r_earth)

float_err = 1e-14
x_range = np.array([m for m in range(-x_upper, x_upper, dx) if abs(m) >= 20])
y_range = np.array([l for l in range(-y_upper, y_upper, dy) if abs(l) >= 20])

for x in x_range: 
    for y in y_range: 
        data = ms.my_simulation(r_pebble = 100, x_h = x, y_h = y) 
        phi_pebble = np.arctan2(data[:, 3], data[:, 2]) 
        xp, yp = rotate_h(data[:, 4], data[:, 5], -phi_pebble) 
        if data.shape == (4000, 6):
            ax.plot(xp , yp, '-', color='gray' , linewidth= 0.5, alpha =0.25)
        else: 
            ax.plot(xp , yp, '-', color='red' , linewidth= 0.5, zorder=10)

ax.set_aspect(1)
ax.set(xlim=[-width, width], ylim=[-width, width])
end = time.time()
print(end-start)