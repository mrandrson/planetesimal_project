import rebound
import matplotlib.pyplot as plt
import imageio.v2 as imageio  # Use imageio version 2
import os
import shutil

# Initialize simulation
sim = rebound.Simulation()
sim.add(m=1)
sim.add(m=1e-2, e=0.041, a=0.4, inc=0.2, f=0.43, Omega=0.82, omega=2.98)
sim.add(m=2e-3, e=0.24, a=1.0, pomega=2.14)
sim.add(m=2e-3, e=0.24, a=1.5, omega=1.14, l=2.1)
op = rebound.OrbitPlot(sim, unitlabel="[AU]", color=True, periastron=True)

# Create a directory to store frames
frames_dir = "frames"
os.makedirs(frames_dir, exist_ok=True)

# Generate and save frames
filenames = []
for i in range(4000):
    op.sim.integrate(sim.t + 0.2)
    op.update(updateLimits=True)  # update data

    # Save frame
    frame_filename = f"{frames_dir}/frame_{i:03d}.png"
    op.fig.savefig(frame_filename)
    filenames.append(frame_filename)

