"""Environment 1: the Multiplicative Oasis.

Inhabitants are whole numbers. Newcomer k arrives on day 21+k and is the number 12+k.
The constitution is a list of generators plus one clause inherited from counting:
every inhabitant is a product of generators in exactly one way.

A newcomer with no factorization over the generators is absurd. A newcomer with
two different factorizations contradicts the constitution.
"""
import sys, json, time, bisect
import numpy as np
from oasis import Environment, run


def nfact(n, gens, gidx, i=0, cap=2):
    """Number of multisets of generators (index >= i) whose product is n, capped."""
    if n == 1:
        return 1
    t = 1 if gidx.get(n, -1) >= i else 0
    for j in range(i, len(gens)):
        g = gens[j]
        if g * g > n:
            break
        if n % g == 0:
            t += nfact(n // g, gens, gidx, j, cap)
            if t >= cap:
                return t
    return t


class Rule:
    def __init__(self, gens):
        self.gens = sorted(gens)
        self.gidx = {g: i for i, g in enumerate(self.gens)}

    def count(self, n):
        return nfact(n, self.gens, self.gidx)

    def add(self, g):
        # generators arrive in increasing order in this desert, so append keeps order
        if self.gens and g < self.gens[-1]:
            self.gens = sorted(self.gens + [g])
            self.gidx = {x: i for i, x in enumerate(self.gens)}
        else:
            self.gidx[g] = len(self.gens)
            self.gens.append(g)

    def copy(self):
        return Rule(list(self.gens))


class MultOasis(Environment):
    name = 'Multiplicative Oasis'
    language = 'products of generators, each in one way'

    def __init__(self, frozen=False):
        self.frozen = frozen

    def founders(self):
        return list(range(2, 12))

    def newcomer(self, k):
        return 12 + k

    def describe(self, rule):
        g = rule.gens
        return f'{len(g)} generators: ' + ', '.join(map(str, g[:6])) + (' ...' if len(g) > 6 else '')

    def constitution(self, record):
        r = Rule([])
        for x in sorted(record):
            if r.count(x) == 0:
                r.add(x)
        return r

    def test(self, rule, item, record):
        c = rule.count(item)
        return 'ok' if c == 1 else ('absurd' if c == 0 else 'contradiction')

    def revise(self, rule, record, backlog):
        if self.frozen:
            return None, 'exhausted'
        new = rule.copy()
        for x in sorted(backlog):
            if new.count(x) == 0:
                new.add(x)
        if any(new.count(x) != 1 for x in backlog):
            return None, 'contradiction'
        return new, 'logical'

    def verify(self, old, new, record, backlog):
        # C1: every factorization over old generators is still a factorization over new ones
        c1 = set(old.gens) <= set(new.gens)
        # C2: only record entries divisible by a new generator can change their count.
        top = record[-1]
        added = sorted(set(new.gens) - set(old.gens))
        c2 = all(new.count(m) == 1 for g in added for m in range(g, top + 1, g))
        # C3: every newcomer that caused the collapse now has exactly one factorization
        c3 = all(new.count(x) == 1 for x in backlog)
        return c1, c2, c3


def sieve(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


def smooth_count_curve(y_primes, N, points):
    """Count y-smooth numbers in [2, N] and return counts at the given points."""
    vals = [1]
    for p in y_primes:
        nxt = []
        for v in vals:
            w = v
            while w <= N:
                nxt.append(w)
                w *= p
        vals = nxt
    vals = np.sort(np.array([v for v in vals if v >= 2], dtype=np.int64))
    return [int(np.searchsorted(vals, x, side='right')) for x in points], len(vals)


def main():
    out = {}
    MAXD = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    # ---------------- main run, C = 0
    t = time.time()
    env = MultOasis()
    res = run(env, 0, MAXD)
    ev = res['events']
    collapses = [e for e in ev if e.kind in ('logical', 'illogical', 'contradiction', 'exhausted')]
    kinds = {}
    for e in collapses:
        kinds[e.kind] = kinds.get(e.kind, 0) + 1
    last_n = 12 + (res['last_day'] - 21)
    s = sieve(last_n)
    primes = np.nonzero(s)[0]
    collapse_items = [int(e.items[0]) for e in collapses]
    founders_primes = [g for g in res['rule'].gens if g <= 11]
    same = collapse_items == [int(p) for p in primes if p > 11]
    allchecks = all(all(e.checks) for e in collapses if e.checks)
    out['run'] = dict(days=res['last_day'], last_number=last_n, collapses=len(collapses), kinds=kinds,
                      founder_generators=founders_primes, collapse_items_equal_primes=same,
                      primes_up_to_last=int(len(primes)), all_C1_C2_C3=allchecks,
                      absurd_total=res['absurd_total'], seconds=round(time.time() - t, 1),
                      first_collapses=[(e.day, int(e.items[0])) for e in collapses[:40]])
    print(json.dumps(out['run'], indent=1)[:1500])

    # ---------------- capacity invariance
    cap = []
    for C in (0, 1, 3, 10, 30, 100):
        env = MultOasis()
        r = run(env, C, 20021)
        col = [e for e in r['events'] if e.kind == 'logical']
        delays = []
        for e in col:
            for x in e.items:
                delays.append(e.day - (int(x) - 12 + 21))
        cap.append(dict(C=C, collapses=len(col), generators=len(r['rule'].gens),
                        gens=list(r['rule'].gens), absurd=r['absurd_total'],
                        max_wait=max(delays) if delays else 0,
                        mean_wait=round(sum(delays) / len(delays), 1) if delays else 0,
                        illogical=sum(1 for e in r['events'] if e.kind != 'logical' and e.kind != 'constitution')))
    base = cap[0]['gens']
    for c in cap:
        g = c.pop('gens')
        c['decided_up_to'] = g[-1]
        c['same_generators_as_C0'] = g == [x for x in base if x <= g[-1]]
    out['capacity'] = cap
    print(cap)

    # ---------------- frozen constitution: generators fixed on day 10
    pts = [int(round(10 ** (k / 8))) for k in range(8, 57)]
    counts, total = smooth_count_curve([2, 3, 5, 7, 11], 10 ** 7, pts)
    out['frozen'] = dict(points=pts, smooth=counts,
                         absurd_share=[1 - c / (x - 1) for c, x in zip(counts, pts)], smooth_total_1e7=total)
    # Euclid's newcomer: product of the founders' generators plus one
    out['euclid_newcomer'] = 2 * 3 * 5 * 7 * 11 + 1
    out['euclid_newcomer_is_prime'] = bool(sieve(2311)[2311])

    # ---------------- day grid for the figure (days 1-130)
    grid = []
    for n in range(2, 12):
        grid.append((1, n, 'founder', bool(s[n])))
    for n in range(12, 122):
        grid.append((n - 12 + 21, n, 'collapse' if s[n] else 'explained', bool(s[n])))
    out['grid'] = grid
    json.dump(out, open('mult.json', 'w'), indent=1)
    print('frozen share at 1e7', out['frozen']['absurd_share'][-1], 'smooth', total)


if __name__ == '__main__':
    main()
