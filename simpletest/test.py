import rebound

sim = rebound.Simulation()
sim.G = 6.67430e-11  # SI units
sim.units = ('s', 'm', 'kg')
tf = 350000
# Add central body (Sun-like mass)
sun_mass = 1.989e30
sim.add(m=sun_mass, x=0, y=0, z=0, vx=0, vy=0, vz=0)

# Add Earth-mass particle at 1 AU with circular velocity
earth_mass = 5.972e24
AU = 1.496e11
v_circ = (sim.G * sun_mass / AU)**0.5
sim.add(m=earth_mass, x=AU, y=0, z=0, vy=v_circ)

sim.integrate(tf*3.154e7)  # 1 year

sim2 = rebound.Simulation()
sim2.G = 6.67430e-11
sim2.units = ('s', 'm', 'kg')

sim2.add(m=earth_mass, x=AU, vy=v_circ)

def custom_sun_force(_):
    for p in sim2.particles:
        dx = -p.x
        dy = -p.y
        dz = -p.z
        r3 = (dx**2 + dy**2 + dz**2)**1.5
        if r3 == 0: continue
        a = sim2.G * sun_mass / r3
        p.ax += a * dx
        p.ay += a * dy
        p.az += a * dz

sim2.additional_forces = custom_sun_force

sim2.integrate(tf*3.154e7)

earth1 = sim.particles[1]
earth2 = sim2.particles[0]

print("Final positions:")
print(f"  Real:   ({earth1.x:.3e}, {earth1.y:.3e})")
print(f"  Custom: ({earth2.x:.3e}, {earth2.y:.3e})")

