import random
import numpy as np

n_max = 100
n_point=100000
maxq = 2 * np.pi / 3.0
minq = 2 * np.pi / 200.0

q = np.zeros(1000000)
N = np.zeros(1000000)
for i in range(1, n_point + 1):
    rand1 = random.random()
    rand2 = random.random()
    rand3 = random.random()
    nx =(n_max + 1) * rand1
    ny =((2 * n_max) + 1) * rand2 - n_max
    nz =((2 * n_max) + 1) * rand3 - n_max
    N2 = (nx * nx + ny * ny + nz * nz)
    if N2 <= (n_max*n_max):
        qx = nx * minq
        qy = ny * minq
        qz = nz * minq
        q2 = qx*qx+qy*qy+qz*qz
        q[i] = np.sqrt(q2)
        a = int(q[i]/0.1)
        N[a] = N[a] + 1

for i in range(1, 401):
    if N[i] != 0:
        q_value = float(i * 0.1)
        print(q_value)