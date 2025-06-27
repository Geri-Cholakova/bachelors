# -*- coding: utf-8 -*-
"""
Created on Wed May  7 17:33:45 2025

@author: Geri
"""
import rebound
import numpy as np
import matplotlib.pyplot as plt
def setup_simulation():
    sim = rebound.Simulation()
    sim.integrator = "ias15"
    sim.collision = "line"
    sim.collision_resolve = "merge"   
    sim.add(m=1., r=4.67e-3)                   # Sun
    sim.add(m=3e-6, r=4.26e-5, a=1)            # Planet
    sim.add(m=3e-38, r=1e-13, a=1.02)          # Pebble
    return sim
n = float(input("Orbit no: "))
def center_on_planet(sim, planet_index=1):
    planet = sim.particles[planet_index]
    for p in sim.particles:
        p.x -= planet.x
        p.y -= planet.y
        p.z -= planet.z
        p.vx -= planet.vx
        p.vy -= planet.vy
        p.vz -= planet.vz
sim = setup_simulation()
sim.integrate(n * 2 * np.pi)
center_on_planet(sim)  # shift coordinate system
x_vals = [p.x for p in sim.particles]
y_vals = [p.y for p in sim.particles]

colors = ['orange', 'blue', 'green']
labels = ['Sun', 'Planet', 'Pebble']

plt.figure(figsize=(6, 6))
for i in range(len(sim.particles)):
    plt.scatter(x_vals[i], y_vals[i], label=labels[i], color=colors[i], s=50)

plt.xlabel("x [AU]")
plt.ylabel("y [AU]")
plt.title(f"Particle Positions after {n} Orbits (Centered on Planet)")
plt.legend()
plt.axis()
plt.grid(False)
plt.show()
