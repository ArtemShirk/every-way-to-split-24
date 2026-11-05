"""Prediction stage: can the Multiplicative Oasis predict its own collapses?

1. Count law. From the record of its first million days the Oasis fits a density
   law for collapses, 1 / (A ln x + B), the way Gauss read his tables. It then
   predicts the number of collapses in windows it has never seen, up to 10^16,
   and the windows are counted with a segmented sieve.
2. Error envelope. pi(x) - li(x) and pi(x) - R(x) up to 10^8, against
   Schoenfeld's bound |pi(x) - li(x)| < sqrt(x) ln x / (8 pi), x >= 2657, which
   holds for all x if and only if the Riemann hypothesis is true.

Usage: python3 predict.py [SIEVE_LIMIT]
"""
import sys, json, time
from math import isqrt, log, sqrt, pi
import numpy as np
from scipy.special import expi
import mpmath

mpmath.mp.dps = 30


def sieve(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    s[4::2] = False
    for i in range(3, isqrt(n) + 1, 2):
        if s[i]:
            s[i * i::2 * i] = False
    return s


def window_count(lo, W, base):
    """Number of primes in [lo, lo + W), lo > sqrt(lo + W)."""
    arr = np.ones(W, dtype=bool)
    hi = lo + W
    lim = isqrt(hi) + 1
    b = base[base <= lim]
    small = b[b < W]
    for p in small.tolist():
        st = (-lo) % p
        arr[st::p] = False
    big = b[b >= W]
    if len(big):
        st = (-lo) % big
        st = st[st < W]
        arr[st] = False
    return int(arr.sum())


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 8
    t0 = time.time()
    s = sieve(N + 1000)
    primes = np.nonzero(s)[0].astype(np.int64)
    print('sieve', len(primes), round(time.time() - t0, 1), flush=True)
    out = dict(N=N, primes=int(len(primes)))

    # ------------------------------------------------ 1. the Oasis writes a density law from days <= 10^6
    blk = 10 ** 4
    edges = np.arange(10 ** 4, 10 ** 6 + 1, blk)
    cnt = np.diff(np.searchsorted(primes, edges))
    mid = edges[:-1] + blk / 2
    inv_density = blk / cnt                     # average gap in the block
    A, B = np.polyfit(np.log(mid), inv_density, 1)
    out['fit'] = dict(A=float(A), B=float(B), blocks=int(len(cnt)))
    print('fit 1/density = %.4f ln x + %.4f' % (A, B), flush=True)

    def fitted(lo, hi):
        f = lambda x: 1 / (A * mpmath.log(x) + B)
        return float(mpmath.quad(f, [lo, hi]))

    li = lambda x: mpmath.li(x)
    R = lambda x: mpmath.riemannr(x)
    W = 10 ** 6
    rows = []
    for k in range(7, 17):
        lo = 10 ** k
        t = time.time()
        actual = window_count(lo, W, primes)
        p_fit = fitted(lo, lo + W)
        p_li = float(li(lo + W) - li(lo))
        p_R = float(R(lo + W) - R(lo))
        rows.append(dict(k=k, lo=lo, actual=actual, fit=p_fit, li=p_li, R=p_R,
                         z_li=(actual - p_li) / sqrt(p_li), err_fit_pct=100 * (actual - p_fit) / actual,
                         err_li_pct=100 * (actual - p_li) / actual, seconds=round(time.time() - t, 1)))
        print(rows[-1], flush=True)
    out['windows'] = rows

    # ------------------------------------------------ 2. the error envelope up to N
    P = primes[primes <= N]
    idx = np.arange(1, len(P) + 1, dtype=np.float64)          # pi(p) at p
    lix = expi(np.log(P.astype(np.float64)))                 # li(x) = Ei(ln x)
    env = np.sqrt(P) * np.log(P) / (8 * pi)
    m = P >= 2657
    r_hi = np.abs(idx[m] - lix[m]) / env[m]                   # at p
    r_lo = np.abs(idx[m] - 1 - lix[m]) / env[m]               # the moment before p arrives
    r_lo[0] = 0.0                                             # x = 2657 is where the bound starts
    rr = np.maximum(r_hi, r_lo)
    j = int(np.argmax(rr))
    out['schoenfeld'] = dict(checked=int(m.sum()), max_ratio=float(rr[j]), at=int(P[m][j]),
                             holds=bool((rr < 1).all()))
    print('schoenfeld', out['schoenfeld'], flush=True)
    # curve samples for the figure
    xs = np.unique(np.round(np.logspace(np.log10(3000), np.log10(N), 400)).astype(np.int64))
    pis = np.searchsorted(P, xs, side='right')
    curve = []
    for x, pv in zip(xs.tolist(), pis.tolist()):
        curve.append((x, pv, float(li(x)), float(R(x)), sqrt(x) * log(x) / (8 * pi)))
    out['curve'] = curve
    # pi(10^k) table
    out['pi_table'] = [(10 ** k, int(np.searchsorted(P, 10 ** k, side='right')),
                        float(li(10 ** k)), float(R(10 ** k))) for k in range(1, 9)]
    out['seconds'] = round(time.time() - t0, 1)
    json.dump(out, open('predict.json', 'w'), indent=1)
    print('done', out['seconds'])


if __name__ == '__main__':
    main()
