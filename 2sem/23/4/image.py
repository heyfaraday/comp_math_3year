#%%
import numpy as np
import matplotlib.pyplot as plt


#%%
# t, x, u = np.genfromtxt("computational_map.dat").T
t, x, u = np.genfromtxt("data.dat").T

Nx = 20
Nt = 20

u = u.reshape(Nt + 1, Nx + 1)
x = x.reshape(Nt + 1, Nx+ 1)
t = t.reshape(Nt + 1, Nx + 1)

plt.figure(figsize=(10, 10))
plt.title("")
plt.pcolor(x, t, u, cmap=plt.cm.jet)
plt.show()

#%% 
t, e1, e2 = np.genfromtxt("error.dat").T

plt.figure(figsize=(10, 10))
e1, = plt.plot(t[1:], e1[1:], 'g--')
plt.show()

plt.figure(figsize=(10, 10))
e2, = plt.plot(t[3:], e2[3:], 'g--')
plt.show()