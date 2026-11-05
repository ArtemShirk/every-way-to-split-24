# Verifies Theorem A: eta(az) eta((24-a)z) = (theta[1,a,6a] - theta[4,a-4,a+1]) / 2
# for a = 1..12, coefficient by coefficient up to N.  Needs only numpy.
import numpy as np
from math import isqrt

def pentagonal(s, N):                 # prod (1 - q^{s n}) = sum_k (-1)^k q^{s k(3k-1)/2}
    e, sg, k = [], [], 0
    while True:
        new = [kk for kk in ((k, -k) if k else (0,)) if s * kk * (3 * kk - 1) // 2 <= N]
        if not new: return np.array(e), np.array(sg)
        for kk in new:
            e.append(s * kk * (3 * kk - 1) // 2); sg.append(-1 if kk % 2 else 1)
        k += 1

def eta_pair(a, b, N):                # coefficients of q * prod (1 - q^{an})(1 - q^{bn})
    c = np.zeros(N + 1, dtype=np.int64)
    E1, S1 = pentagonal(a, N); E2, S2 = pentagonal(b, N)
    for e2, s2 in zip(E2, S2):
        i = 1 + E1 + e2; m = i <= N
        np.add.at(c, i[m], S1[m] * s2)
    return c

def theta(A, B, C, N):                # number of (x, y) with A x^2 + B x y + C y^2 = n
    t = np.zeros(N + 1, dtype=np.int64); D = 4 * A * C - B * B
    R = isqrt(4 * max(A, C) * N // D) + 2 + isqrt(N) * (abs(B) + 1)
    for y in range(-isqrt(4 * A * N // D) - 2, isqrt(4 * A * N // D) + 3):
        x = np.arange(-R, R + 1, dtype=np.int64)
        v = A * x * x + B * x * y + C * y * y
        np.add.at(t, v[(v >= 0) & (v <= N)], 1)
    return t

N = 100000
for a in range(1, 13):
    lhs = 2 * eta_pair(a, 24 - a, N)
    rhs = theta(1, a, 6 * a, N) - theta(4, a - 4, a + 1, N)
    print(a, 'identity holds for n <=', N, ':', bool((lhs == rhs).all()))
