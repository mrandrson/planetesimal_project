import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.animation as animation
from scipy.interpolate import interp1d
from matplotlib.patches import Circle

cs = np.sqrt(3.36e4)  # Sound speed in m/s (as per your previous definition)

# Load all position files
data_list = []

for f in os.listdir('.'):
    if f.startswith('positions_') and f.endswith('.txt'):
        dat = pd.read_csv(f, delim_whitespace=True, header=None, names=["time", "x", "y", "particle_id"])
        dat["v0"] = float(f.split('_')[1].replace('.txt',''))  # Extract v0 from filename
        data_list.append(dat)

if not data_list:
    raise ValueError("No position files found.")

df = pd.concat(data_list, ignore_index=True)
df = df.sort_values(by=["time"])

# Create a denser time grid for smooth animation
original_times = df["time"].unique()
interpolated_times = np.linspace(df["time"].min(), df["time"].max(), len(original_times) * 25)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-3e15, 3e15)
ax.set_ylim(-3e15, 3e15)
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title("Particle Motion Over Time")

scat = ax.scatter([], [], s=10)
line_traces = {}
interpolated_paths = {}

# Interpolate positions for each unique (v0, particle_id) combination
for (v0_val, pid), pdat in df.groupby(["v0", "particle_id"]):
    pdat = pdat.drop_duplicates(subset=["time"]).sort_values(by="time").reset_index(drop=True)
    if len(pdat) < 2:
        continue  # Need at least two points to interpolate

    interp_x = interp1d(pdat["time"], pdat["x"], kind='linear', fill_value="extrapolate", assume_sorted=True)
    interp_y = interp1d(pdat["time"], pdat["y"], kind='linear', fill_value="extrapolate", assume_sorted=True)
    interpolated_paths[(v0_val, pid)] = {
        "x": interp_x(interpolated_times),
        "y": interp_y(interpolated_times)
    }
    line_traces[(v0_val, pid)], = ax.plot([], [], linewidth=1, alpha=0.6)

def update(frame):
    t = interpolated_times[frame]
    xs, ys = [], []

    for key, path in interpolated_paths.items():
        x, y = path["x"][frame], path["y"][frame]
        xs.append(x)
        ys.append(y)
        line_traces[key].set_data(path["x"][:frame], path["y"][:frame])

    scat.set_offsets(np.c_[xs, ys])
    ax.set_title(f"Time = {3.171e-8 * t:.2e} yr")

ani = animation.FuncAnimation(fig, update, frames=len(interpolated_times), interval=50)

plt.show()




