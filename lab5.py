# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

m = 9.11e-28
e = 4.8e-10
c = 3e10
erstd = 4*np.pi/1e3
sgse = np.sqrt(4*np.pi*8.854e-12)
B = -100 * erstd
E = 10e-1 * sgse

C = (m, e, c, B, E)

# function
def system(u, t, C):
    dvxdt = (C[1]*C[4] + C[1]/C[2] * u[1]*C[3])/C[0]
    dvydt = -C[1]/C[2]/C[0] * u[0]*C[3]
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
xi = c*E/B
c3 = -xi
c2 = 5 + c3/2
c1 = 5 - c3/2

c1 = 10
c2 = xi

def analytic_vx(t):
    return c1 * np.cos(g*t) + c2 * np.sin(g*t)

def analytic_vy(t):
    return c2 * np.cos(g*t) - c1 * np.sin(g*t) + c3

def analytic_rx(t):
    return 10/g * np.sin(g*t) - xi/g * np.cos(g*t) + xi/g + 1

def analytic_ry(t):
    return xi/g * np.sin(g*t) + 10/g * np.cos(g*t) - xi*t + 1 - 10/g

# plot trajectory
plt.plot(sol[:, 2], sol[:, 3], label='$\mathbf{r}\,(t)$')
plt.plot(analytic_rx(t), analytic_ry(t), ':', 
         label='$\mathbf{r}\,(t)\:analytic$',
         linewidth=3)
plt.xlabel('$x$ [см]')
plt.ylabel('$y$ [см]')
plt.legend()
# plt.savefig('plotTrajectory.pdf')
plt.show()
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
# plt.savefig('plotSpeed.pdf')
plt.show()
plt.clf()

#plot.mistake
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.yscale("log")
plt.plot(t, abs(analytic_vx(t)-sol[:, 0])/sol[:, 0], 'r-+', label='$\delta v_x$')
plt.plot(t, abs(analytic_vy(t)-sol[:, 1])/sol[:, 1], 'm-+', label='$\delta v_y$')
plt.plot(t, abs(analytic_rx(t)-sol[:, 2])/sol[:, 2], '-+', label='$\delta r_x$')
plt.plot(t, abs(analytic_ry(t)-sol[:, 3])/sol[:, 3], '-+', label='$\delta r_y$')
plt.legend()
plt.savefig('plotMistake.pdf')