"""Language 3: pentagonal products.

Euler's pentagonal theorem: prod_{n>=1} (1 - q^n) = sum_k (-1)^k q^{k(3k-1)/2}, k over all integers.
So eta(a z) eta(b z) with a + b = 24 is the sparse double sum

    q * prod (1 - q^{a n}) (1 - q^{b n}) = sum_{j,k} (-1)^{j+k} q^{1 + a j(3j-1)/2 + b k(3k-1)/2}

and its coefficients can be computed exactly by adding +-1 at a few hundred
thousand positions. No zeta function, no class field theory, no floating point.

This script:
  1. computes the coefficients for every split 24 = a + b,
  2. checks the D = 27 rule three ways for every prime up to 10^6:
     p = x^2 + 27y^2   <=>   2 is a cube mod p (p = 1 mod 3)   <=>   coefficient of q^p equals 2,
  3. searches all reduced binary quadratic forms with |discriminant| <= 800 for the one whose
     primes are exactly the primes with coefficient 2, for each split a + b = 24.
"""
import sys, json, time
from math import isqrt, gcd
import numpy as np
from env_forms import sieve, represented


def pentagonal(scale, L):
    e, s = [], []
    k = 0
    while True:
        added = False
        for kk in ((k, -k) if k else (0,)):
            v = scale * kk * (3 * kk - 1) // 2
            if v <= L:
                e.append(v); s.append(-1 if kk % 2 else 1); added = True
        if not added:
            break
        k += 1
    return np.array(e, dtype=np.int64), np.array(s, dtype=np.int64)


def eta_pair(a, b, L):
    c = np.zeros(L + 1, dtype=np.int64)
    E1, S1 = pentagonal(a, L)
    E2, S2 = pentagonal(b, L)
    for e2, s2 in zip(E2.tolist(), S2.tolist()):
        idx = 1 + E1 + e2
        m = idx <= L
        np.add.at(c, idx[m], S1[m] * s2)
    return c, len(E1) * len(E2)


def reduced_forms(maxdisc):
    out = []
    for N in range(3, maxdisc + 1):
        if N % 4 not in (0, 3):
            continue
        Dl = -N
        A = 1
        while 3 * A * A <= N:
            for B in range(-A, A + 1):
                if (B * B - Dl) % (4 * A):
                    continue
                C = (B * B - Dl) // (4 * A)
                if C < A:
                    continue
                if B < 0 and (abs(B) == A or A == C):
                    continue
                if gcd(gcd(A, abs(B)), C) != 1:
                    continue
                out.append((A, B, C, Dl))
            A += 1
    return out


def form_primes(A, B, C, L, isprime):
    Dl = B * B - 4 * A * C
    N = -Dl
    xb = isqrt(4 * C * L // N) + 2
    yb = isqrt(4 * A * L // N) + 2
    x = np.arange(-xb, xb + 1, dtype=np.int64)[:, None]
    y = np.arange(0, yb + 1, dtype=np.int64)[None, :]
    v = (A * x * x + B * x * y + C * y * y).ravel()
    v = v[(v >= 2) & (v <= L)]
    v = np.unique(v)
    return v[isprime[v]]


def fmt_form(A, B, C):
    def term(c, mono, first):
        if c == 0:
            return ''
        sign = '' if first else (' + ' if c > 0 else ' − ')
        if first and c < 0:
            sign = '−'
        mag = abs(c)
        return sign + ('' if mag == 1 else str(mag)) + mono
    s = term(A, 'x²', True) + term(B, 'xy', False) + term(C, 'y²', False)
    return s


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 6
    L2 = int(sys.argv[2]) if len(sys.argv) > 2 else 2 * 10 ** 5
    t0 = time.time()
    s = sieve(L)
    primes = np.nonzero(s)[0]
    out = {}

    # ---------------- 1. the D = 27 rule three ways
    c, pairs = eta_pair(6, 18, L)
    rep = represented(27, L)
    ps = primes[primes >= 5]
    a_p = c[ps]
    cubic = np.array([(p % 3 == 1) and pow(2, (p - 1) // 3, int(p)) == 1 for p in ps.tolist()])
    form = rep[ps]
    eta2 = a_p == 2
    out['d27'] = dict(primes_checked=int(len(ps)), L=L, pairs=pairs,
                      form_count=int(form.sum()), cubic_count=int(cubic.sum()), eta_count=int(eta2.sum()),
                      form_eq_cubic=bool((form == cubic).all()), form_eq_eta=bool((form == eta2).all()),
                      values={str(v): int((a_p == v).sum()) for v in np.unique(a_p).tolist()},
                      minus1_rule=bool(((a_p == -1) == ((ps % 3 == 1) & ~form)).all()),
                      zero_rule=bool(((a_p == 0) == (ps % 3 == 2)).all()),
                      first_coeffs=c[:200].tolist(),
                      first_forms=[int(p) for p in ps[form][:24]])
    print('d27', {k: v for k, v in out['d27'].items() if k != 'first_coeffs'}, round(time.time() - t0, 1), flush=True)

    # ---------------- 2. the whole family a + b = 24
    s2 = s[:L2 + 1]
    pr2 = primes[primes <= L2]
    forms = reduced_forms(800)
    small = 4000
    fp_small = {}
    fam = []
    for a in range(1, 13):
        b = 24 - a
        cc, _ = eta_pair(a, b, L2)
        coeff = cc[pr2]
        bad = {p for p in range(2, 30) if (6 * a * b) % p == 0}
        keep = np.array([int(p) not in bad for p in pr2])
        S = set(pr2[keep & (coeff == 2)].tolist())
        vals = sorted(set(coeff[keep].tolist()))
        # multiplicative check: a(m n) = a(m) a(n) for small coprime m, n
        mult = all(cc[m * n] == cc[m] * cc[n] for m in range(2, 60) for n in range(m + 1, 60)
                   if gcd(m, n) == 1 and m * n <= L2)
        match = None
        nmatch = 0
        Ssmall = {p for p in S if p <= small}
        for (A, B, C, Dl) in forms:
            key = (A, B, C)
            if key not in fp_small:
                fp_small[key] = set(form_primes(A, B, C, small, s2).tolist())
            badD = bad | {p for p in range(2, 1000) if (-Dl) % p == 0}
            got = {p for p in fp_small[key] if p not in badD}
            want = {p for p in Ssmall if p not in badD}
            if got != want:
                continue
            full = set(form_primes(A, B, C, L2, s2).tolist())
            full = {p for p in full if p not in badD}
            wantf = {p for p in S if p not in badD}
            if full == wantf:
                nmatch += 1
                if match is None:
                    match = dict(form=fmt_form(A, B, C), A=A, B=B, C=C, disc=Dl, primes=len(wantf))
        fam.append(dict(a=a, b=b, multiplicative=mult, values=vals[:12], count2=len(S), match=match, matches=nmatch))
        print(a, b, 'mult', mult, 'values', vals[:8], 'a_p=2 count', len(S), 'match', match, flush=True)
    out['family'] = fam
    out['L2'] = L2
    out['seconds'] = round(time.time() - t0, 1)
    json.dump(out, open('eta.json', 'w'), indent=1, ensure_ascii=False)


if __name__ == '__main__':
    main()
