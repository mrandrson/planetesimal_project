import rebound
import matplotlib.pyplot as plt
import imageio.v2 as imageio  # Use imageio version 2
import os
import shutil
import numpy as np
sim = rebound.Simulation()
sim.add(["Sun","Mercury","Venus","Earth","Mars"],date="2018-02-10 00:00")
earlier_M = np.radians(51.158) - np.radians(100)  # Adjust this value as needed
sim.add(a=-1.2723, e=1.20113, inc=np.radians(122.74),
        Omega=np.radians(24.597), omega=np.radians(241.811),
        M=earlier_M)

rebound.OrbitPlotSet(sim, color=True, xlim=[-3,3], ylim=[-3,3]);
op = rebound.OrbitPlot(sim, unitlabel="[AU]", color=True, periastron=True)

frames_dir = "frames"
os.makedirs(frames_dir, exist_ok=True)

filenames = []
for i in range(10000):
    op.sim.integrate(sim.t + 0.1)
    op.update(updateLimits=True)  # update data

    frame_filename = f"{frames_dir}/frame_{i:03d}.png"
    op.fig.savefig(frame_filename)
    filenames.append(frame_filename)


