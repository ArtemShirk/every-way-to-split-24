#!/usr/bin/env python3
"""The Prime Oasis: run every stage and print a readable report.

IDLE: open this file and press F5 to run everything (about a minute).
To run one stage: Run > Run... Customized (Shift+F5) and type one of
    1        multiplication: every collapse lands on a prime
    2        Z[sqrt(-5)]: the first illogical collapse (6 = 2*3 = (1+sqrt-5)(1-sqrt-5))
    3        which primes are x^2 + Dy^2: congruences, then polynomials
    24       the 24-split rule, checked for all twelve splits
    predict  the Oasis fits its own counting law and forecasts blind windows
    zeros    rebuild the prime powers from zeta zeros
    all      everything (same as pressing F5)
Add a second word "full" to predict for the run out to 10^16 (needs about 1 GB of memory).
Terminal: python3 run_all.py 3
"""
import os, sys, time
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
import numpy as np
from math import isqrt, log, sqrt, pi

from oasis import run


def header(title):
    print()
    print('=' * 72)
    print(title)
    print('=' * 72)


# ---------------------------------------------------------------- stage 1
def stage1(days=20021):
    from env_mult import MultOasis, sieve, smooth_count_curve
    header(f'STAGE 1  Multiplication  ({days:,} days)')
    t = time.time()
    res = run(MultOasis(), 0, days)
    col = [e for e in res['events'] if e.kind not in ('constitution',)]
    last = 12 + days - 21
    primes = [int(p) for p in np.nonzero(sieve(last))[0] if p > 11]
    items = [int(e.items[0]) for e in col]
    print('Founders (day 1): 2..11.  Constitution (day 10): generators 2, 3, 5, 7, 11.')
    print('First collapses (day: newcomer that could not be built):')
    print('   ' + ',  '.join(f'day {e.day}: {e.items[0]}' for e in col[:10]))
    print(f'Collapses: {len(col):,}, all logical: {all(e.kind == "logical" for e in col)}, '
          f'C1-C3 passed every time: {all(all(e.checks) for e in col)}')
    print(f'The collapsing numbers are exactly the primes 13..{primes[-1]:,}: {items == primes}')
    print()
    print('Capacity test (2,021 days each): bigger C waits longer, adds the same primes')
    base = None
    for C in (0, 3, 10, 30):
        r = run(MultOasis(), C, 2021)
        g = r['rule'].gens
        base = base or g
        same = g == [x for x in base if x <= g[-1]]
        n = sum(1 for e in r['events'] if e.kind == 'logical')
        print(f'   C = {C:>2}: {n:>3} collapses, same primes added: {same}')
    counts, total = smooth_count_curve([2, 3, 5, 7, 11], 10 ** 6, [10 ** 6])
    print()
    print(f'Frozen constitution (never adds a generator): at 1,000,000 it can build only '
          f'{counts[0]:,} numbers, {100 * (1 - counts[0] / (10 ** 6 - 1)):.2f}% absurd.')
    print("Euclid's newcomer 2*3*5*7*11 + 1 = 2311 has no generator dividing it.")
    print(f'({time.time() - t:.1f} s)')


# ---------------------------------------------------------------- stage 2
def stage2(norm=10000):
    from env_ring import RingOasis, fmt, NORM
    header(f'STAGE 2  Z[sqrt(-5)]  (all numbers a + b*sqrt(-5) up to norm {norm:,})')
    t = time.time()
    env = RingOasis(norm)
    res = run(env, 0, 10 ** 9)
    ev = res['events']
    x, fs = env.first_contradiction
    day = [e.day for e in ev if e.kind == 'contradiction'][0]
    print('Founders: ' + ', '.join(fmt(e) for e in env.founders()))
    print(f'Logical collapses before the break: {sum(1 for e in ev if e.kind == "logical" and e.day < day)}')
    print(f'Day {day}: newcomer {fmt(x)} is built two ways from atoms that cannot be split:')
    for f in fs:
        print('      ' + ' x '.join(f'({fmt(a)})' for a in f))
    print('   No list of atoms fixes it. The engine ascends to ideals (Kummer 1847, Dedekind 1871).')
    after = sum(1 for e in ev if e.kind == 'contradiction' and e.day > day)
    print(f'After the ascent: {len(env.els) - 10 - (day - 20):,} more days, {after} contradictions, '
          f'{env.ideal_checks:,} factorizations checked against ideals.')
    print('In ideals:  2 = P2*P2   3 = P3*P3\'   1+sqrt-5 = P2*P3   1-sqrt-5 = P2*P3\'')
    print("            so 6 = P2*P2*P3*P3' either way.")
    print(f'({time.time() - t:.1f} s)')


# ---------------------------------------------------------------- stage 3
def stage3(limit=10 ** 6, M=1000):
    from env_forms import sieve, represented, FormOasis
    header(f'STAGE 3  Which primes are x^2 + Dy^2?  (all primes to {limit:,})')
    t = time.time()
    primes = np.nonzero(sieve(limit))[0]
    for D in (1, 2, 3, 5, 27, 64):
        env = FormOasis(D, primes, represented(D, limit), M)
        res = run(env, 0, 10 ** 9)
        ev = res['events']
        path = [e.note for e in ev if e.kind in ('constitution', 'logical', 'ascent')]
        ex = [e.day for e in ev if e.kind == 'exhausted']
        last = max(e.day for e in ev if e.kind in ('constitution', 'logical', 'ascent'))
        n = len(env.items) - 10
        print(f'D = {D:>2}: founders wrote {path[0]:<10} final law: {env.describe(res["rule"]):<16} '
              f'collapses: {sum(1 for e in ev if e.kind == "logical"):>2}')
        if ex:
            print(f'        every modulus up to {M} refuted by day {ex[0]}; ascended to polynomials '
                  f'after rejecting {env.rejected} smaller ones')
        print(f'        no errors in the {21 + n - 1 - max(last, 20):,} primes after the last change')
    print('Rules found: D=1 mod 4, D=2 mod 8, D=3 mod 3, D=5 mod 20 (Fermat, Euler, Lagrange);')
    print('D=27: x^3 - 2 has 3 roots mod p (Euler, Gauss); D=64: x^4 - 2 has 4 roots mod p (Euler, Gauss).')
    print(f'({time.time() - t:.1f} s)')


# ---------------------------------------------------------------- the 24 rule
def theta(A, B, C, N):
    t = np.zeros(N + 1, dtype=np.int64); D = 4 * A * C - B * B
    R = isqrt(4 * max(A, C) * N // D) + 2 + isqrt(N) * (abs(B) + 1)
    for y in range(-isqrt(4 * A * N // D) - 2, isqrt(4 * A * N // D) + 3):
        x = np.arange(-R, R + 1, dtype=np.int64)
        v = A * x * x + B * x * y + C * y * y
        np.add.at(t, v[(v >= 0) & (v <= N)], 1)
    return t


def stage24(N=100000):
    from eta import eta_pair
    from env_forms import sieve
    header(f'THE 24 RULE  eta(az) eta((24-a)z) = (theta[1,a,6a] - theta[4,a-4,a+1]) / 2   (n <= {N:,})')
    t = time.time()
    s = sieve(N)
    print(' a +  b   identity   primes p (not dividing ab) with coefficient 2 = primes of x^2+axy+6ay^2')
    for a in range(1, 13):
        b = 24 - a
        c = eta_pair(a, b, N)[0]
        t0, t4 = theta(1, a, 6 * a, N), theta(4, a - 4, a + 1, N)
        ident = bool((2 * c == t0 - t4).all())
        ps = [p for p in np.nonzero(s)[0].tolist() if (a * b) % p]
        rule = all((c[p] == 2) == (t0[p] > 0) for p in ps)
        first = [p for p in ps if c[p] == 2][:6]
        print(f'{a:>2} + {b:>2}   {str(ident):<9}  {str(rule):<5}  first: {first}')
    print(f'({time.time() - t:.1f} s)')


# ---------------------------------------------------------------- prediction
def stage_predict(full=False):
    import mpmath
    from predict import sieve, window_count
    top = 16 if full else 14
    header(f'PREDICTION  The Oasis writes its own counting law, then forecasts to 10^{top}')
    t = time.time()
    base = 10 ** 8 if full else 10 ** 7
    primes = np.nonzero(sieve(base + 1000))[0].astype(np.int64)
    edges = np.arange(10 ** 4, 10 ** 6 + 1, 10 ** 4)
    cnt = np.diff(np.searchsorted(primes, edges))
    A, B = np.polyfit(np.log(edges[:-1] + 5000), 10 ** 4 / cnt, 1)
    print(f'Law fitted on days up to 1,000,000:  average gap = {A:.4f} ln x {B:+.4f}')
    print('Window of 1,000,000 numbers starting at   forecast    counted   off by')
    for k in range(7, top + 1):
        lo = 10 ** k
        f = float(mpmath.quad(lambda x: 1 / (A * mpmath.log(x) + B), [lo, lo + 10 ** 6]))
        n = window_count(lo, 10 ** 6, primes)
        print(f'   10^{k:<2}                                {f:>9,.0f}  {n:>9,}   {100 * abs(n - f) / n:.2f}%')
    print(f'({time.time() - t:.1f} s)')


# ---------------------------------------------------------------- zeros
def stage_zeros():
    from zeros import load_zeros, psi_true, detect
    header('ZEROS  Rebuild the prime powers from zeta zeros alone')
    t = time.time()
    g_all = load_zeros()
    X = 700
    _, lam = psi_true(np.array([X + 1.0]))
    for K in (100, 1000, len(g_all)):
        g = g_all[:K]
        taper = np.cos(np.pi * g / (2 * g[-1] * 1.02)) ** 2
        n, J = detect(g, X, taper)
        wrong = n[(lam[n] > 0) != (J > log(2) / 2)]
        clean = int(wrong[0] - 1) if len(wrong) else X
        print(f'{K:>5,} zeros: every n read correctly up to {clean}')
    print(f'({time.time() - t:.1f} s)')


STAGES = {'1': stage1, '2': stage2, '3': stage3, '24': stage24, 'predict': stage_predict, 'zeros': stage_zeros}

if __name__ == '__main__':
    args = sys.argv[1:] or ['all']
    if args[0] == 'all':
        for f in STAGES.values():
            f()
    elif args[0] == 'predict':
        stage_predict(full=len(args) > 1 and args[1] == 'full')
    elif args[0] in STAGES:
        STAGES[args[0]]()
    else:
        print(__doc__)
