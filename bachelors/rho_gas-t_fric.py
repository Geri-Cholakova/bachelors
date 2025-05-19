# -*- coding: utf-8 -*-
"""
Created on Wed May 14 16:05:21 2025

@author: Geri
"""
import math
from astropy import astropy.constants
#all in SI
AU = 149597870700 
M =  1.988416*10**30
G = 6.67430/10**11
k_b = 1.380649/10**23
m = 1.66/10**27
x = float(input("x = "))
y = float(input("y = "))
z = float(input("z = "))
a = float(input("a = [m]"))
T = 200*AU**0.5
#T is prop to r**-0.5
Sigma = 2000*AU
#prop to r**-1
c_s = math.sqrt(k_b*T/(2.3*m))
#prop to r**-0.25
Omega = math.sqrt(G*M)
#prop to r**-1.5
h = c_s/Omega
#h is prop to r**1.25
rho_coeff = Sigma / (math.sqrt(2*math.pi)*h)
#coeff is prop to r**-2.25
rho_exp = -0.5/h
#prop to r**-1.25, z**2
t_coeff = 1500*a/(rho_coeff*math.sqrt(8*k_b*T/(math.pi*2.3*m)))
#prop to r**2.5
t_exp = -rho_exp
#prop to r**-1.25, z**2
t_fric = t_coeff*(x**2 + y**2)**1.25 * math.exp(t_exp*z**2/(x**2 + y**2)**0.625)
rho_gen = rho_coeff/(x**2+y**2)**(1.125) * math.exp((rho_exp*z**2/(x**2+y**2)**0.625))
print(f'rho is {rho_coeff} *(x**2+y**2)**(-0.8125) * exp( {rho_exp}*r**(-1.25)*z**2')
print(f't_fric is {t_coeff}*(x**2 + y**2)**0.9375*{t_exp}*r**-1.25*z**2')
print(t_fric, rho_gen)
          