import rebound
import matplotlib.pyplot as plt
import imageio.v2 as imageio  # Use imageio version 2
import os
import shutil
import numpy as np
# Initialize simulation
sim = rebound.Simulation()
sim.add(m=1)
sim.add(a=4, e=0, m=3e-6)
sim.add(a=10, e=.6, m=3e-6)
rebound.OrbitPlotSet(sim, color=True);
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

