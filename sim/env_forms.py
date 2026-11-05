"""Environment 3: the Form Oasis.

Inhabitants are primes, arriving in order. Each carries one fact from the
record: whether it can be written as x^2 + D*y^2. The constitution must
predict that fact for every newcomer before it is checked.

Language 1  congruences: "p is of the form iff p mod m lies in a set A", m <= M
Language 2  polynomials: "p is of the form iff f(x) has deg f roots mod p"
            (f monic, irreducible, degree 2 to 4, coefficients in [-3, 3])
Language 3  pentagonal products (see eta.py): "iff the coefficient of q^p is 2"

Usage: python3 env_forms.py [LIMIT] [M]
"""
import sys, json, time
from math import isqrt
import numpy as np
from oasis import Environment, run


def sieve(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


def represented(D, L):
    rep = np.zeros(L + 1, dtype=bool)
    for y in range(0, isqrt(L // D) + 1):
        r = L - D * y * y
        xs = np.arange(0, isqrt(r) + 1, dtype=np.int64)
        rep[xs * xs + D * y * y] = True
    return rep


# ---------------------------------------------------------------- polynomials mod p
def polymulmod(u, v, f, p):
    """u, v: coefficient lists (low to high) of degree < n; f monic of degree n (low to high, length n+1)."""
    n = len(f) - 1
    w = [0] * (2 * n - 1)
    for i, a in enumerate(u):
        if a:
            for j, b in enumerate(v):
                w[i + j] += a * b
    for k in range(2 * n - 2, n - 1, -1):
        c = w[k] % p
        if c:
            for i in range(n):
                w[k - n + i] -= c * f[i]
        w[k] = 0
    return [x % p for x in w[:n]]


def splits_completely(f, p):
    """True if monic f (low->high) has deg f distinct roots mod p, i.e. x^p = x mod (f, p).
    Assumes p does not divide disc(f)."""
    n = len(f) - 1
    result = [1] + [0] * (n - 1)
    base = [0, 1] + [0] * (n - 2)
    e = p
    while e:
        if e & 1:
            result = polymulmod(result, base, f, p)
        base = polymulmod(base, base, f, p)
        e >>= 1
    return result == [0, 1] + [0] * (n - 2)


def disc(f):
    """Discriminant of a monic polynomial of degree 2, 3 or 4 (coefficients low to high)."""
    n = len(f) - 1
    if n == 2:
        c, b, _ = f
        return b * b - 4 * c
    if n == 3:
        d, c, b, _ = f
        return 18 * b * c * d - 4 * b ** 3 * d + b * b * c * c - 4 * c ** 3 - 27 * d * d
    e, d, c, b, _ = f   # x^4 + b x^3 + c x^2 + d x + e
    a = 1
    return (256 * a**3 * e**3 - 192 * a**2 * b * d * e**2 - 128 * a**2 * c**2 * e**2 + 144 * a**2 * c * d**2 * e
            - 27 * a**2 * d**4 + 144 * a * b**2 * c * e**2 - 6 * a * b**2 * d**2 * e - 80 * a * b * c**2 * d * e
            + 18 * a * b * c * d**3 + 16 * a * c**4 * e - 4 * a * c**3 * d**2 - 27 * b**4 * e**2
            + 18 * b**3 * c * d * e - 4 * b**3 * d**3 - 4 * b**2 * c**3 * e + b**2 * c**2 * d**2)


def has_rational_root(f):
    c0 = f[0]
    if c0 == 0:
        return True
    for r in range(1, abs(c0) + 1):
        if c0 % r == 0:
            for s in (r, -r):
                if sum(a * s ** i for i, a in enumerate(f)) == 0:
                    return True
    return False


def irreducible(f):
    n = len(f) - 1
    if n in (2, 3):
        if n == 2:
            dd = disc(f)
            return not (dd >= 0 and isqrt(dd) ** 2 == dd)
        return not has_rational_root(f)
    # quartic: no rational root and no factorization into two monic integer quadratics
    if has_rational_root(f):
        return False
    e, d, c, b, _ = f
    # (x^2 + p x + q)(x^2 + r x + s): q s = e, p + r = b, q + s + p r = c, p s + q r = d
    for q in range(-abs(e), abs(e) + 1):
        if q == 0 or e % q:
            continue
        s = e // q
        for pp in range(-12, 13):
            r = b - pp
            if q + s + pp * r == c and pp * s + q * r == d:
                return False
    return True


def poly_candidates(H=3):
    cands = []
    for n in (2, 3, 4):
        rng = range(-H, H + 1)
        if n == 2:
            coefs = [(c, b) for b in rng for c in rng]
        elif n == 3:
            coefs = [(d, c, b) for b in rng for c in rng for d in rng]
        else:
            coefs = [(e, d, c, b) for b in rng for c in rng for d in rng for e in rng]
        for cf in coefs:
            f = list(cf) + [1]
            if disc(f) == 0 or not irreducible(f):
                continue
            cands.append(f)
    key = lambda f: (len(f), max(abs(x) for x in f[:-1]), sum(abs(x) for x in f[:-1]), tuple(reversed(f[:-1])))
    cands.sort(key=key)
    return cands


def fmt_poly(f):
    n = len(f) - 1
    terms = []
    for i in range(n, -1, -1):
        c = f[i]
        if c == 0:
            continue
        mono = '' if i == 0 else ('x' if i == 1 else f'x{"²³⁴"[i - 2]}')
        if i == n:
            terms.append(mono)
            continue
        sign = '+' if c > 0 else '−'
        mag = abs(c)
        body = (str(mag) if (mag != 1 or i == 0) else '') + mono
        terms.append(f'{sign} {body}')
    return ' '.join(terms)


# ---------------------------------------------------------------- rules
class Cong:
    kind = 'congruence'

    def __init__(self, m, P, Lab):
        self.m = m
        r = P % m
        self.table = {}
        for rr, l in zip(r.tolist(), Lab.tolist()):
            self.table.setdefault(rr, l)

    def predict(self, p):
        return self.table.get(p % self.m)


class Poly:
    kind = 'polynomial'

    def __init__(self, f, index):
        self.f, self.index = f, index
        self.bad = abs(disc(f))

    def predict(self, p):
        if self.bad % p == 0:
            return None
        return splits_completely(self.f, p)


def consistent_mod(m, P, Lab):
    r = P % m
    cnt = np.bincount(r, minlength=m)
    ones = np.bincount(r, weights=Lab.astype(np.float64), minlength=m)
    return bool(np.all((ones == 0) | (ones == cnt)))


class FormOasis(Environment):
    language = 'congruences'

    def __init__(self, D, primes, rep, M):
        self.D, self.M = D, M
        ps = [int(p) for p in primes if (2 * D) % int(p) != 0]
        self.items = [(p, bool(rep[p])) for p in ps]
        self.name = f'Form Oasis D={D}'
        self.cands = None
        self.history = []          # (day, prime, rule description)
        self.rejected = 0

    def founders(self):
        return self.items[:10]

    def newcomer(self, k):
        return self.items[10 + k] if 10 + k < len(self.items) else None

    def describe(self, rule):
        if rule.kind == 'congruence':
            return f'mod {rule.m}'
        return f'{fmt_poly(rule.f)} splits'

    def _arrays(self, record):
        P = np.array([p for p, _ in record], dtype=np.int64)
        Lab = np.array([l for _, l in record], dtype=bool)
        return P, Lab

    def constitution(self, record):
        P, Lab = self._arrays(record)
        for m in range(1, self.M + 1):
            if consistent_mod(m, P, Lab):
                self.history.append((10, record[-1][0], m))
                return Cong(m, P, Lab)

    def test(self, rule, item, record):
        p, lab = item
        if rule.kind == 'congruence':
            r = p % rule.m
            if r not in rule.table:
                rule.table[r] = lab       # the rule makes no claim here yet; extend it
                return 'ok'
            return 'ok' if rule.table[r] == lab else 'absurd'
        pred = rule.predict(p)
        if pred is None:
            return 'ok'
        return 'ok' if pred == lab else 'absurd'

    def revise(self, rule, record, backlog):
        P, Lab = self._arrays(record)
        if rule.kind == 'congruence':
            for m in range(rule.m + 1, self.M + 1):
                if consistent_mod(m, P, Lab):
                    self.history.append((None, backlog[-1][0], m))
                    return Cong(m, P, Lab), 'logical'
            return None, 'exhausted'
        nxt = self._poly_search(record, rule.index + 1)
        if nxt is None:
            return None, 'exhausted'
        return nxt, 'logical'

    def _poly_search(self, record, start=0):
        if self.cands is None:
            self.cands = poly_candidates()
        for i in range(start, len(self.cands)):
            f = self.cands[i]
            rule = Poly(f, i)
            ok = True
            for p, lab in record:
                pr = rule.predict(p)
                if pr is not None and pr != lab:
                    ok = False
                    break
            if ok:
                return rule
            self.rejected += 1
        return None

    def ascend(self, rule, record, backlog):
        if rule.kind != 'congruence':
            return None
        new = self._poly_search(record)
        if new is None:
            return None
        self.history.append((None, backlog[-1][0], fmt_poly(new.f)))
        return new, f'polynomial language: {fmt_poly(new.f)} has {len(new.f) - 1} roots mod p'

    def _agrees(self, rule, item):
        p, lab = item
        if rule.kind == 'congruence':
            pr = rule.table.get(p % rule.m)
            return pr is None or pr == lab
        pr = rule.predict(p)
        return pr is None or pr == lab

    def verify(self, old, new, record, backlog):
        if new.kind == 'congruence':
            P, Lab = self._arrays(record)
            c2 = consistent_mod(new.m, P, Lab)
            c1 = c2
            c3 = all(new.table.get(p % new.m) == lab for p, lab in backlog)
            return c1, c2, c3
        c2 = all(self._agrees(new, it) for it in record)
        c1 = c2
        c3 = all(self._agrees(new, it) for it in backlog)
        return c1, c2, c3


def refutation_map(D, primes, rep, M):
    """For each modulus m <= M: the first prime (in arrival order) that refutes the congruence rule mod m."""
    ps = np.array([p for p in primes if (2 * D) % int(p) != 0], dtype=np.int64)
    lab = rep[ps]
    out = []
    for m in range(1, M + 1):
        r = ps % m
        _, first = np.unique(r, return_index=True)
        firstlab = np.zeros(m, dtype=bool)
        firstlab[r[first]] = lab[first]
        bad = lab != firstlab[r]
        out.append(int(ps[np.argmax(bad)]) if bad.any() else 0)
    return out


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 6
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    s = sieve(L)
    primes = np.nonzero(s)[0]
    out = dict(L=L, M=M, runs={})
    for D in (1, 2, 3, 5, 27, 64):
        t = time.time()
        rep = represented(D, L)
        env = FormOasis(D, primes, rep, M)
        res = run(env, 0, 10 ** 9)
        ev = res['events']
        hist = []
        for e in ev:
            if e.kind in ('constitution', 'logical', 'exhausted', 'ascent', 'illogical', 'end'):
                hist.append(dict(day=e.day, kind=e.kind, by=[x for x in e.items][-1:] if e.items else [],
                                 rule=e.note, checks=e.checks))
        final = res['rule']
        last_change = max(e.day for e in ev if e.kind in ('constitution', 'logical', 'ascent'))
        n_items = len(env.items)
        out['runs'][D] = dict(history=hist, final=env.describe(final), final_kind=final.kind,
                              newcomers=n_items - 10, last_day=21 + n_items - 11,
                              last_change_day=last_change,
                              quiet_after=21 + n_items - 11 - last_change,
                              collapses=sum(1 for e in ev if e.kind in ('logical', 'illogical')),
                              exhausted=[e.day for e in ev if e.kind == 'exhausted'],
                              poly_rejected=env.rejected, seconds=round(time.time() - t, 1),
                              final_classes=(sorted(r for r, l in final.table.items() if l)
                                             if final.kind == 'congruence' else None),
                              errors_after_last_change=0 if not [e for e in ev if e.day > last_change
                                                                 and e.kind not in ('constitution',)] else -1)
        print(D, out['runs'][D]['final'], 'collapses', out['runs'][D]['collapses'],
              'exhausted', out['runs'][D]['exhausted'], 'last change day', last_change,
              'rejected polys', env.rejected, out['runs'][D]['seconds'], 's', flush=True)
        if D in (5, 27, 64):
            out['runs'][D]['refutation'] = refutation_map(D, primes, rep, M)
    json.dump(out, open('forms.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
