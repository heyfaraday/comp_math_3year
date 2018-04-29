#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
# t, x, u = np.genfromtxt("computational_map.dat").T
t, x, u = np.genfromtxt("computational_map_explicit.dat").T

Nx = 20
Nt = 2000

u = u.reshape(Nt + 1, Nx + 1)
x = x.reshape(Nt + 1, Nx+ 1)
t = t.reshape(Nt + 1, Nx + 1)

plt.figure(figsize=(10, 10))
plt.title("explicit, 1/20, 1/2000")
plt.pcolor(x, t, u)
plt.show()

#%%
t, u, u1, u2= np.genfromtxt("solutions_at_the_middle.dat").T

plt.figure(figsize=(10, 10))
a, = plt.plot(t, u, 'b--', label = 'crank-nicolson')
e, = plt.plot(t, u1, 'r--', label = 'analytical')
ex, = plt.plot(t, u2, 'g--', label = 'explicit')
first_legend = plt.legend(handles=[a, e, ex], loc=3)
plt.show()

