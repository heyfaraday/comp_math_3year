#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
t, x, u = np.genfromtxt("4/error.dat").T

Nx = 200
Nt = 200

u = u.reshape(Nt + 1, Nx + 1)
t = t.reshape(Nt + 1, Nx + 1)
x = x.reshape(Nt + 1, Nx + 1)

plt.figure(figsize=(7, 7))
plt.pcolor(x, t, u)
plt.show()

#%%
x, y, u = np.genfromtxt("4/bin/data.dat").T

N = 100

u = u.reshape(N + 1, N + 1)
x = x.reshape(N + 1, N + 1)
y = y.reshape(N + 1, N + 1)

plt.figure(figsize=(7, 7))
plt.pcolor(x, y, u)
plt.show()

#%%
t, e1, e2 = np.genfromtxt("4/bin/error.dat").T

plt.figure(figsize=(10, 10))
e1, = plt.plot(t[1:], e1[1:], 'g--')
plt.show()

plt.figure(figsize=(10, 10))
e2, = plt.plot(t[3:], e2[3:], 'g--')
plt.show()