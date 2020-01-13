# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

m = 9.11e-28
e = 4.8e-10
c = 3e10
erstd = 0.01256637
sgse = 3.3e-5
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

# plot trajectory
plt.plot(sol[:, 2], sol[:, 3])
plt.xlabel('$x$ [см]')
plt.ylabel('$y$ [см]')
plt.title('Тракетория за интервал $t$ = [0, 1e-3] сек')
plt.savefig('plotTrajectory.pdf')
plt.clf()

# plot speed
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.plot(sol[:, 0], t, label='$v_x$')
plt.plot(sol[:, 1], t, label='$v_y$')
plt.xlabel('$t$ [сек]')
plt.ylabel('$v$ [см/с]')
plt.title('Зависимость скорости от времени')
plt.legend()
plt.savefig('plotSpeed.pdf')