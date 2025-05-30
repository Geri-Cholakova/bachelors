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


all_data = []  # a list

initial_data = ms.my_simulation(r_pebble=1, x_h=0, y_h=100)
# here the assumption is that initial data doesn't have a collision; that in diff simulations the difference in planet placement is negligible
all_data.append(initial_data)

phi = np.arctan2(initial_data[:, 3], initial_data[:, 2])


def rotate_h(x, y, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    return c * x - s * y - 999.9985, s * x + c * y


# rotation matrix [cosa -sina] [sina cosa]
# the fkt returns 2 elements

width = 50
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
    x_values = np.linspace(1, 51, 50)
    y_values = np.linspace(1, 51, 50)
    s = 1  # radius of the pebble, in m
    # make long array with all parameter combinations, here x,z
    xy = np.array(np.meshgrid(x_values, y_values, indexing="ij")).reshape(2, -1).T

    with Pool(4) as p:
        simulations = list(tqdm(p.imap(ms.worker, xy), total=len(xy)))
        for data in simulations:
            phi_pebble = np.arctan2(data[:, 3], data[:, 2])
            xp, yp = rotate_h(data[:, 4], data[:, 5], -phi_pebble)
            if data.shape == (50000, 6):
                ax.plot(
                    xp,
                    yp,
                    "-",
                    color="gray",
                    linewidth=0.5,
                    alpha=0.25,
                    rasterized=True,
                )
            else:
                ax.plot(
                    xp, yp, "-", color="red", linewidth=0.5, zorder=10, rasterized=False
                )
    ax.set_aspect(1)
    ax.set(xlim=[-width, width], ylim=[-width, width])
    plt.savefig("Pebble_with_vdust_xy_1-51_rH_dx_1_rH.pdf", dpi=300)
