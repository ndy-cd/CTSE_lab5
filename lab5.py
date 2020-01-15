# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

m = 9.11e-28
e = 4.8e-10
c = 3e10
erstd = 4*np.pi/1e3
sgse = np.sqrt(4*np.pi*8.854e-12)
B = 10 * erstd
E = 10 * sgse

C = (m, e, c, B, E)

# function
def system(u, t, C):
    dvxdt = (C[1]*C[4] + C[1]/C[2] * u[1]*C[3])/C[0]
    dvydt = C[1]/C[2]/C[0] * u[0]*C[3]
    drxdt = u[0]
    drydt = u[1]
    return [dvxdt, dvydt, drxdt, drydt]

# initial condition
u0 = [10, 0, 1, 1]

# time points
t = np.linspace(0, 1e-6, 101)

# solve
sol = odeint(system, u0, t, args=(C,))
print(sol)

#analytic
g = e*B/m/c
c3 = -c*E/B
c2 = 5 + c3/2
c1 = 5 - c3/2

def analytic_vx(t):
    return c1 * np.exp(g*t) + c2 * np.exp(-g*t)

def analytic_vy(t):
    return c1 * np.exp(g*t) - c2 * np.exp(-g*t) + c3

# plot trajectory
plt.plot(sol[:, 2], sol[:, 3])
plt.xlabel('$x$ [см]')
plt.ylabel('$y$ [см]')
plt.savefig('plotTrajectory.pdf')
plt.clf()

# plot speed
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.plot(t, sol[:, 0], label='$v_x$')
plt.plot(t, sol[:, 1], label='$v_y$')
plt.plot(t, analytic_vx(t), ':', label='$v_x$ analytic', linewidth=3)
plt.plot(t, analytic_vy(t), ':', label='$v_y$ analytic', linewidth=3)
plt.xlabel('$t$ [сек]')
plt.ylabel('$v$ [см/с]')
plt.legend()
plt.savefig('plotSpeed.pdf')
plt.clf()

#plot.mistake
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.plot(t, abs(analytic_vx(t)-sol[:, 0])/sol[:, 0], 'r-+', label='$\delta v_x$')
plt.plot(t, abs(analytic_vy(t)-sol[:, 1])/sol[:, 1], 'm-+', label='$\delta v_y$')
plt.legend()
plt.savefig('plotMistake.pdf')