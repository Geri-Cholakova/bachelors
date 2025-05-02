import rebound
sim = rebound.Simulation()
sim.add(m=1)
sim.add(m=1e3, x=1, vy=1)
sim.integrate(10000)
sim.status()

