# -*- coding: utf-8 -*-
"""
Created on Sun May 18 15:56:02 2025

@author: Geri
"""

import matplotlib.pyplot as plt
import numpy as np
import my_simulation as ms
from multiprocessing.pool import Pool
from tqdm import tqdm
from astropy import units as u
from functools import partial


all_data = []  # a list
force = True

initial_data = ms.my_simulation(r_pebble=1, x_h=0, y_h=10, force = False, check_distance = False)
# here the assumption is that initial data doesn't have a collision; that in diff simulations the difference in planet placement is negligible
if initial_data.shape == (50001, 6):
    initial_data = initial_data[:-1, :]
    all_data.append(initial_data)
else:
    print("Pick another starting value")
    
phi = np.arctan2(initial_data[:, 3], initial_data[:, 2])


def rotate_h(x, y, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    return c * x - s * y - 99.9985, s * x + c * y


# rotation matrix [cosa -sina] [sina cosa]
# the fkt returns 2 elements

width = 5
x0, y0 = rotate_h(initial_data[:, 0], initial_data[:, 1], -phi)
x1, y1 = rotate_h(initial_data[:, 2], initial_data[:, 3], -phi)
f, ax = plt.subplots()
ax.plot(x0, y0, "k+-")
ax.plot(x1, y1, "k--")
r_hill = plt.Circle((0, 0), 1, color="c")
r_earth = plt.Circle((0, 0), 4.267e-2, color="b")
ax.add_patch(r_hill)
ax.add_patch(r_earth)

if __name__ == "__main__":
    number_x_values = 46
    number_y_values = 11
    x_values = np.linspace(1.0, 1.9, number_x_values)
    y_values = np.linspace(9, 10, number_y_values)
    # make long array with all parameter combinations, here x,z
    xy = np.array(np.meshgrid(x_values, y_values, indexing="ij")).reshape(2, -1).T

    norm_s = 2.0 * np.pi/(1*u.year).si.value # in rebound unit * s^-1
    rho_s = 1500 # in kg/m**3
    s = 1  # radius of the pebble, in m
    worker_with_s = partial(ms.worker, r_pebble = s, E_coef = np.sqrt(8 / np.pi) / (rho_s * s) /norm_s, force = force)
    #I am only incl change of parameter r_pebble, but the worker function could have other parameters changed as well
    histogr_data = np.zeros((number_x_values*number_y_values, 3))
    histogr_data[:, 0] = xy[:, 0]
    histogr_data[:, 1] = xy[:, 1]

    with Pool(4) as p:
        #this is with Epstein force
        simulations = list(tqdm(p.imap(worker_with_s, xy), total=len(xy)))
        for i, data in enumerate(simulations):
            phi_pebble = np.arctan2(data[:-1, 3], data[:-1, 2])
            xp, yp = rotate_h(data[:-1, 4], data[:-1, 5], -phi_pebble)
            if bool(data[-1, 1]) == False:
                ax.plot(xp,yp,"-",color="gray",linewidth=0.5,alpha=0.25,rasterized=True,)
            else:
                ax.plot(xp, yp, "-", color="red", linewidth=0.5, zorder=10, rasterized=False)
            histogr_data[i, 2] = data[-1,1]

    ax.set_aspect(1)
    ax.set(xlim=[-width, width], ylim=[-2*width, 2*width])
    filename = f"r_pebble_{s:.4f}_x_{x_values[0]:.2f}_{x_values[-1]:.2f}_y_{y_values[0]:.2f}_{y_values[-1]:.2f}_delx_{x_values[1]-x_values[0]}".replace('.', 'p')
    # Result: x1p00_y2p50

    if force: 
        plt.savefig(f"Epstein_{filename}_rH.pdf", dpi=300,)
    else: 
        plt.savefig(f"Noforce_{filename}_rH.pdf", dpi=300,)
        
    bins = [37, 1]
    range2d = [[x_values[0], x_values[-1]], [y_values[0], y_values[-1]]]
    
    
    H_total, xedges, yedges = np.histogram2d(histogr_data[:, 0], histogr_data[:,1], bins=bins, range=range2d)
    H_coll, _, _ = np.histogram2d(histogr_data[:,0], histogr_data[:,1], bins=[xedges, yedges], weights=histogr_data[:,2])
    
    with np.errstate(divide='ignore', invalid='ignore'):
        collision_fraction = np.nan_to_num(H_coll / H_total)

    plt.figure(figsize=(6, 5))
    plt.imshow(collision_fraction.T, origin='lower', extent=[
        xedges[0], xedges[-1], yedges[0], yedges[-1]],
        cmap='Greys', interpolation='nearest', aspect=0.1)
    plt.xlabel("x start (Hill units)")
    plt.ylabel("y start (Hill units)")
    plt.title("Collision Fraction per Starting Position")
    plt.colorbar(label="Collision fraction")
    plt.tight_layout()
    plt.savefig(f"heatmap_{filename}.pdf", dpi=300)
    plt.show()
