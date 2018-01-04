import numpy as np

dim = 1000
alpha = 0.01
beta = 10.0

A = np.empty([dim, dim])

A[0][0] = 2.0 + beta
for i in range(1, dim):
    A[i][i] = 2.0 + alpha
    A[i][i - 1] = - 1.0
for i in range(0, dim - 1):
    A[i][i + 1] = - 1.0

v = np.linalg.eigvals(A)

print '\lambda_{max}:', v.max()
print '\lambda_{min}:', v.min()