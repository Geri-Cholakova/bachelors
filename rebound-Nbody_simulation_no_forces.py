import rebound
import numpy as np #to use pi
sim = rebound.Simulation()
#creating a simulation
a = 2
b = 1
c = 2
r = np.sqrt(a**2+b**2+c**2)
sim.add(m=1.)                # Central object
sim.add(m=1e-3, x=a, y=b, z=c, vx=1/np.sqrt(3)*a*r**0.5, vy=1/np.sqrt(3)*b*r**0.5, vz=1/np.sqrt(3)*c*r**0.5) #Jupyter mass planet 
sim.add(m=1e-5, x=1, vy=1)# 3* Earth mass planet
#sim.add(m=1e-4, a=1.6, e=0.1)# 30*Earth mass planet
sim.move_to_com()
#to not have a drift because of center of mass movement
t=float(input("iterations of the simulation? (local years): ")) #enter how many years 
sim.integrate(2.*np.pi*t)
#seeing the simulation after time t=100 for 2*pi*t=1 yr for a=1, M_central=1
sim.status()
#seeing the data on which basis the plot is made
rebound.OrbitPlot(sim)
#plot the orbital data after t iterations
