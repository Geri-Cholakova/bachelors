import rebound
import numpy as np #we need it so that we can use pi later
def setupSimulation():
    sim = rebound.Simulation()
    sim.integrator = "ias15" # IAS15 is the default integrator, so we don't need this line
    sim.collision = "line" #this scales as o(n^2), use a diff one for a bigger simulation
    #sim.collision_resolve = "merge" #just for example
    sim.add(m=1.,r=4.67e-3 )  #parameters of the Sun
    sim.add(m=3e-6,r=4.26e-5, a=1) #parameters of Earth
    sim.add(m=3e-32, r=1e-11, a=1.02) #example dust stone w radius several km
    return sim
sim = setupSimulation()
t=float(input("iterations of the simulation? (local years): ")) #enter how many years
sim.integrate(t*2.*np.pi) #we have collision between 925 and 950 iterations
sim.move_to_com()
for o in sim.orbits():
    print(o)
    rebound.OrbitPlot(sim)
