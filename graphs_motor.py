import subfunctions as sf
import numpy as np
import matplotlib.pyplot as plt

planet, power_subsys, science_payload, chassis, motor, speed_reducer, wheel, wheel_assembly, rover = sf.create_dictionary()

omega_NL = .38
omega = np.linspace(0, omega_NL, 100)
tau = sf.tau_dcmotor(omega, motor)
power = omega * tau

fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
fig.subplots_adjust(hspace = .5)
fig.tight_layout()
fig.suptitle('Motor Shaft Speed / Torque / Power Curves', fontsize = 10)

ax1.plot(tau, omega)
ax1.set_xlabel('Shaft Torque [Nm]')
ax1.set_ylabel('Shaft Speed [rad/s]')


ax2.plot(tau, power)
ax2.set_xlabel('Shaft Torque [Nm]')
ax2.set_ylabel('Shaft Power [W]')

ax3.plot(omega, power)
ax3.set_xlabel('Shaft Speed (rad/s)')
ax3.set_ylabel('Shaft Power [W]')

print(power)