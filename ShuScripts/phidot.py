import energyplot as ep
import numpy as np
import matplotlib.pyplot as plt

h = 1e7
phi_dot = np.zeros(len(ep.t)-2)  
energy_dot = np.zeros(len(ep.t)-2)
ke_dot = np.zeros(len(ep.t)-2)
for i in range(1, len(ep.t)-2):
    phi_dot[i] = (ep.getPhi(ep.meters_to_au(ep.r[i]), ep.seconds_to_years(ep.t[i]+h))-ep.getPhi(ep.meters_to_au(ep.r[i]), ep.seconds_to_years(ep.t[i]-h)))/(2*h)
for i in range(1, len(ep.t)-1):
    #phi_dot[i-1] = (ep.PE[i+1] - ep.PE[i-1]) / (ep.t[i+1] - ep.t[i-1])
    energy_dot[i-1] = (ep.Ecluster[i+1] - ep.Ecluster[i-1]) / (ep.t[i+1] - ep.t[i-1])
    ke_dot[i-1] = (ep.KE[i+1] - ep.KE[i-1]) / (ep.t[i+1] - ep.t[i-1])

plt.figure(figsize=(21, 6))

plt.subplot(1, 3, 1)
plt.plot(ep.t[1:-1], phi_dot, label="$\dot{\phi}$")
plt.ylim(min(phi_dot), max(phi_dot))
plt.xlabel('Time (s)')
plt.ylabel('d$\Phi$/dt ($m^2/s^3$)')
plt.title('Derivative of Potential Energy')
plt.legend()

plt.subplot(1, 3, 2)
plt.ylim(min(phi_dot), max(phi_dot))
plt.plot(ep.t[1:-1], energy_dot, label="$\dot{E}_{tot}$")
plt.xlabel('Time (s)')
plt.ylabel('d$E_{tot}$/dt ($m^2/s^3$)')
plt.title('Derivative of Total Energy')
plt.legend()

plt.subplot(1, 3, 3)
plt.plot(ep.t[1:-1], ke_dot, label="$\dot{K}$")
#plt.ylim(-7e-8, 7e-8)
plt.xlabel('Time (s)')
plt.ylabel('d$K$/dt ($m^2/s^3$)')
plt.title('Derivative of Kinetic Energy')
plt.legend()

plt.tight_layout()
plt.show()
