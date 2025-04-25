import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.animation as animation
from scipy.interpolate import interp1d
from matplotlib.patches import Circle

cs = np.sqrt(3.36e4)

data_list = []
file_counter = 0

for f in os.listdir('.'):
    if f.endswith('.csv'):
        dat = pd.read_csv(f)
        dat = dat.dropna(axis=1, how='all')
        if not dat.empty:
            dat["particle_id"] = dat["particle_id"] + file_counter * 1000000
            file_counter += 1
            data_list.append(dat)

if not data_list:
    raise ValueError("No CSV files found in the directory.")

df = pd.concat(data_list, ignore_index=True)
df = df.sort_values(by=["time"])
original_time_steps = df["time"].unique()
interpolated_time_steps = np.linspace(df["time"].min(), df["time"].max(), num=len(original_time_steps) * 10)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-1e17, 1e17)
ax.set_ylim(-1e17, 1e17)
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title("Particle Motion Over Time")
scat = ax.scatter([], [], s=10)
line_traces = {}
interpolated_positions = {}

for p in df["particle_id"].unique():
    pdat = df[df["particle_id"] == p].copy()
    pdat = pdat.drop_duplicates(subset=["particle_id", "time"]).reset_index(drop=True)
    if len(pdat) < 3:
        continue
    if any(pdat['energy'] < 0) and any(np.diff(pdat['energy']) > 10): 
        interp_x = interp1d(pdat["time"], pdat["x"], kind="quadratic", fill_value="extrapolate", assume_sorted=True)
        interp_y = interp1d(pdat["time"], pdat["y"], kind="quadratic", fill_value="extrapolate", assume_sorted=True)
        interpolated_positions[p] = {"x": interp_x(interpolated_time_steps), "y": interp_y(interpolated_time_steps)}
        line_traces[p], = ax.plot([], [], linewidth=1.0, alpha=0.7)

circle = Circle((0, 0), 0, fill=False, color='maroon', linewidth=1.5)
ax.add_patch(circle)

def update(frame):
    current_time = interpolated_time_steps[frame]
    x_positions = []
    y_positions = []
    for p in interpolated_positions:
        x_positions.append(interpolated_positions[p]["x"][frame])
        y_positions.append(interpolated_positions[p]["y"][frame])
        line_traces[p].set_data(
            interpolated_positions[p]["x"][:frame],
            interpolated_positions[p]["y"][:frame]
        )
    scat.set_offsets(np.c_[x_positions, y_positions])
    circle.set_radius(cs * current_time)
    ax.set_title(f"Planetesimal Motion at Time t={3.171e-8*current_time:.2e} (yr)")

ani = animation.FuncAnimation(fig, update, frames=len(interpolated_time_steps), interval=50)
#ani.save('/Users/richardanderson/Downloads/shuproblemplanetesimals.mp4')
plt.show()
