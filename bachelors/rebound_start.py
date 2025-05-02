import rebound
import matplotlib.pyplot as plt
sim = rebound.Simulation()
sim.add(m=1)
sim.add(m=1e-3, x=1, vy=1)
sim.add(m=2e-3, x=2, vy=0.66)
sim.integrate(800000)
fig = rebound.OrbitPlot(sim)
plt.show()
