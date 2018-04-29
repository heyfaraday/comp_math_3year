#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
# t, x, u = np.genfromtxt("analytical_solution.dat").T
t, x, u = np.genfromtxt("numerical_solution.dat").T
# t, x, u = np.genfromtxt("error.dat").T

Nx = 400
Nt = 1800

u = u.reshape(Nt + 1, Nx + 1)
x = x.reshape(Nt + 1, Nx+ 1)
t = t.reshape(Nt + 1, Nx + 1)

#%%

plt.figure(figsize=(10, 10))
plt.pcolor(x, t, u)
plt.title("breather, numerical, omega = 0.4, boundary conditions")
plt.show()

#%%
plt.figure(figsize=(10, 10))
plt.title("breather, numerical, omega = 0.4, boundary conditions, n = 400/800/1400")
a = plt.plot(x[400], u[400], 'b--', label = 'analytical')
a = plt.plot(x[800], u[800], 'r--', label = 'analytical')
a = plt.plot(x[1400], u[1400], 'g--', label = 'analytical')

plt.show()
