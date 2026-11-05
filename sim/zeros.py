"""Prediction of positions: rebuild the collapse record from the zeros of zeta.

Riemann-von Mangoldt explicit formula (x > 1, not a prime power):
    psi(x) = x - sum_rho x^rho / rho - ln(2 pi) - (1/2) ln(1 - x^-2)
psi(x) = sum of ln p over prime powers p^k <= x. It jumps by ln p at every prime power.
Using only the first K zeros rho = 1/2 + i gamma, the jump at n is read off as
    J_K(n) = psi_K(n + 1/2) - psi_K(n - 1/2).
A prime power shows a jump near ln p; anything else shows a jump near 0.
"""
import json, sys
import numpy as np
from math import log


def load_zeros():
    z = {}
    for f in ('z_odd.json', 'z_even.json', 'z_odd2.json', 'z_even2.json'):
        try:
            z.update({int(k): v for k, v in json.load(open(f)).items()})
        except FileNotFoundError:
            pass
    ks = sorted(z)
    assert ks == list(range(1, len(ks) + 1)), 'zero list has gaps'
    return np.array([z[k] for k in ks])


def psi_K(x, g, taper=None):
    x = np.asarray(x, dtype=np.float64)
    lx = np.log(x)[:, None]
    rho = 0.5 + 1j * g[None, :]
    w = np.ones_like(g) if taper is None else taper
    terms = np.exp(rho * lx) / rho * w[None, :]
    s = 2 * terms.real.sum(axis=1)
    return x - s - log(2 * np.pi) - 0.5 * np.log(1 - x ** -2.0)


def psi_true(x):
    x = np.asarray(x, dtype=np.float64)
    N = int(x.max()) + 1
    lam = np.zeros(N + 1)
    sv = np.ones(N + 1, dtype=bool); sv[:2] = False
    for i in range(2, N + 1):
        if sv[i]:
            sv[i * i::i] = False
            pk = i
            while pk <= N:
                lam[pk] = log(i); pk *= i
    cum = np.cumsum(lam)
    return cum[np.floor(x).astype(int)], lam


def detect(g, X, taper=None, chunk=400):
    n = np.arange(2, X + 1)
    J = np.empty(len(n))
    for i in range(0, len(n), chunk):
        nn = n[i:i + chunk].astype(np.float64)
        J[i:i + chunk] = psi_K(nn + 0.5, g, taper) - psi_K(nn - 0.5, g, taper)
    return n, J


def main():
    g_all = load_zeros()
    print('zeros', len(g_all), 'largest', g_all[-1])
    Xmax = 1200
    _, lam = psi_true(np.array([Xmax + 1.0]))
    out = dict(zeros=len(g_all), gamma_max=float(g_all[-1]))
    rows = []
    for K in (10, 30, 100, 300, 1000, 2000, 3000, 4000):
        if K > len(g_all):
            break
        g = g_all[:K]
        for name, taper in (('plain', None), ('tapered', np.cos(np.pi * g / (2 * g[-1] * 1.02)) ** 2)):
            n, J = detect(g, Xmax, taper)
            truth = lam[n] > 0
            guess = J > log(2) / 2
            wrong = n[truth != guess]
            clean = int(wrong[0] - 1) if len(wrong) else Xmax
            # among the clean range, tell primes from higher prime powers by the size of the jump
            m = (n <= clean) & guess
            is_prime_guess = J[m] > 0.75 * np.log(n[m])     # a prime jumps by ln n, a power p^k by ln n / k
            is_prime_true = np.array([all(k % d for d in range(2, int(k ** 0.5) + 1)) for k in n[m]])
            primes_ok = bool((is_prime_guess == is_prime_true).all())
            rows.append(dict(K=K, method=name, clean_up_to=clean, primes_in_range=int(is_prime_true.sum()),
                             prime_vs_power_ok=primes_ok, max_err_in_clean=float(np.abs(J[n <= clean] - lam[n[n <= clean]]).max())))
            print(rows[-1], flush=True)
    out['detection'] = rows
    # curves for the figure
    xs = np.linspace(1.6, 60, 2400)
    tru, _ = psi_true(xs)
    curves = {'x': xs.tolist(), 'true': tru.tolist()}
    for K in (3, 30, min(1000, len(g_all))):
        curves[str(K)] = psi_K(xs, g_all[:K]).tolist()
    out['curves'] = curves
    Kmax = len(g_all)
    gK = g_all[:Kmax]
    tap = np.cos(np.pi * gK / (2 * gK[-1] * 1.02)) ** 2
    n, J = detect(gK, 200, tap)
    out['spikes'] = dict(K=Kmax, n=n.tolist(), J=J.tolist(), lam=lam[n].tolist())
    json.dump(out, open('zeros.json', 'w'))


if __name__ == '__main__':
    main()
