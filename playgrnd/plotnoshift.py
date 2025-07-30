import pandas as pd
import matplotlib.pyplot as plt
from ast import literal_eval
import re

df = pd.read_csv('dtnoshift.csv')


df.columns = df.columns.str.strip()


def fix_energy_string(s):
    return literal_eval(re.sub(r'(?<=\d)\s+(?=-?\d)', ', ', s.strip()))

df['Energy'] = df['Energy'].apply(fix_energy_string)


time_values = []
energy_values = []
for t, energies in zip(df['Time'], df['Energy']):
    for e in energies:
        time_values.append(t)
        energy_values.append(e)

dt = df['dt']


fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.scatter(time_values, energy_values, color='tab:red', label='Energy')
ax1.set_yscale('symlog', linthresh=1e0)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Energy', color='tab:red')
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.grid(True, linestyle='--', linewidth=0.5)

ax2 = ax1.twinx()
ax2.set_yscale('log')
ax2.plot(df['Time'], dt, color='tab:blue', linewidth=2, label='dt')
ax2.set_ylabel('dt (s)', color='tab:blue')
ax2.tick_params(axis='y', labelcolor='tab:blue')

fig.tight_layout()
plt.show()


