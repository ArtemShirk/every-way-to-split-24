"""Check the 24-split rule for every split up to LIMIT:
coefficient of q^p in eta(az)eta((24-a)z) is 2  <=>  p is represented by the principal form of discriminant -a(24-a)
(for primes p not dividing 6a(24-a)). Also records what the other coefficient values mean."""
import sys, json
import numpy as np
from env_forms import sieve
from eta import eta_pair, form_primes
L = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 6
s = sieve(L); primes = np.nonzero(s)[0]
def kron(d, p):  # Legendre symbol (d/p) for odd prime p
    d %= p
    return 0 if d == 0 else (1 if pow(d, (p - 1) // 2, p) == 1 else -1)
rows = []
for a in range(1, 13):
    b = 24 - a; N = a * b
    A, B, C = (1, 1, (N + 1) // 4) if N % 4 == 3 else (1, 0, N // 4)
    c, _ = eta_pair(a, b, L)
    bad = set(p for p in range(2, 50) if (6 * N) % p == 0)
    ps = np.array([p for p in primes.tolist() if p not in bad], dtype=np.int64)
    ap = c[ps]
    rep = np.zeros(L + 1, dtype=bool); rep[form_primes(A, B, C, L, s)] = True
    two = ap == 2
    ok2 = bool((two == rep[ps]).all())
    ch = np.array([kron(-N, int(p)) for p in ps.tolist()])
    ok0 = bool(((ap == 0) == (ch == -1)).all())
    vals = {str(v): int((ap == v).sum()) for v in np.unique(ap).tolist()}
    rows.append(dict(a=a, b=b, disc=-N, form=(A, B, C), primes=int(len(ps)), rule_two=ok2, zero_iff_inert=ok0,
                     values=vals, first=[int(p) for p in ps[two][:8]]))
    print(rows[-1], flush=True)
json.dump(rows, open('family.json', 'w'), indent=1)
