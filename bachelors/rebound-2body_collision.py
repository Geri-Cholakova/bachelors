import rebound
import numpy as np #we need it so that we can use pi later
def setupSimulation():
    sim = rebound.Simulation()
    sim.integrator = "ias15" # IAS15 is the default integrator, so we don't need this line
    sim.collision = "line" #this scales as o(n^2), use a diff one for a bigger simulation
    #sim.collision_resolve = "merge" #just for example
    sim.add(m=1.,r=0.05 ) #for m in solar masses, a,r in AU 
    sim.add(m=1e-3,r=0.005, a=1.)
    sim.add(m=3e-6,r=0.0005, a=1.1)
    sim.move_to_com()
    return sim
sim = setupSimulation()
t=float(input("iterations of the simulation? (local years): ")) #enter how many years
sim.integrate(t*2.*np.pi) #we have collision between 200 and 250 iterations 
for o in sim.orbits():
    print(o)
    rebound.OrbitPlot(sim)