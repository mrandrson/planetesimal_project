import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
import energyplot as ep

Rmax = ep.au_to_meters(ep.rmax)

def plot_and_save(x_data, y_data, index, output_dir):
    x_data = np.atleast_1d(x_data)  # Ensure x_data is at least 1D
    y_data = np.atleast_1d(y_data)  # Ensure y_data is at least 1D

    plt.figure()
    plt.scatter(x_data, y_data)
    plt.title(f'Scatter Plot {index}')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.xlim(-Rmax, 3e15)
    plt.ylim(-Rmax, 3e15)

    circle = plt.Circle((0, 0), Rmax, fill=False)
    plt.gca().add_artist(circle)

    filename = f'plot_{index:03}.png'
    plt.savefig(os.path.join(output_dir, filename))
    plt.close()

# Example usage
data_dir = './testposition'
output_dir = './testframe'
gif_output_dir = os.path.expanduser('~/Downloads/planetesimal_project/')
os.makedirs(output_dir, exist_ok=True)
os.makedirs(gif_output_dir, exist_ok=True)

for i in range(0, 99):  # Adjust the range according to the number of file pairs
    x_file = f'x{i}.txt'
    y_file = f'y{i}.txt'

    # Construct full file paths
    x_path = os.path.join(data_dir, x_file)
    y_path = os.path.join(data_dir, y_file)

    # Read data and plot
    x_data = np.loadtxt(x_path)
    y_data = np.loadtxt(y_path)

    plot_and_save(x_data, y_data, i, output_dir)

# Create GIF using imageio
gif_path = os.path.join(gif_output_dir, 'animation.gif')
with imageio.get_writer(gif_path, mode='I', duration=0.1) as writer:
    for i in range(0, 99):  # Adjust the range according to the number of file pairs
        filename = os.path.join(output_dir, f'plot_{i:03}.png')
        image = imageio.imread(filename)
        writer.append_data(image)

