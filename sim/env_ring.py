"""Environment 2: the Oasis of Z[sqrt(-5)].

Inhabitants are numbers a + b*sqrt(-5) with whole a, b, counted up to sign.
They arrive in order of size (norm a^2 + 5b^2); among equal norms, whole
numbers first. The founders' constitution inherits the clause from the
Multiplicative Oasis: every inhabitant is a product of atoms in exactly one way.

Element language: the constitution is a list of atoms (irreducible inhabitants).
Ideal language (parent): the constitution is the list of prime ideals; an
inhabitant's record is how many times each prime ideal divides it.
"""
import sys, json, time
from functools import lru_cache
from oasis import Environment, run

NORM = lambda a, b: a * a + 5 * b * b


def canon(a, b):
    if a < 0 or (a == 0 and b < 0):
        return (-a, -b)
    return (a, b)


def mul(x, y):
    a, b = x; c, d = y
    return (a * c - 5 * b * d, a * d + b * c)


def div(x, y):
    """x / y if it lies in Z[sqrt(-5)], else None."""
    a, b = x; c, d = y
    n = c * c + 5 * d * d
    re, im = a * c + 5 * b * d, b * c - a * d
    if re % n or im % n:
        return None
    return (re // n, im // n)


def fmt(x):
    a, b = x
    if b == 0:
        return str(a)
    s = '√−5' if abs(b) == 1 else f'{abs(b)}√−5'
    if a == 0:
        return ('' if b > 0 else '−') + s
    return f'{a} {"+" if b > 0 else "−"} {s}'


def elements(X):
    out = []
    bmax = int((X / 5) ** 0.5) + 1
    for b in range(-bmax, bmax + 1):
        for a in range(0, int(X ** 0.5) + 2):
            if NORM(a, b) <= X and canon(a, b) == (a, b) and NORM(a, b) > 1:
                out.append((a, b))
    out.sort(key=lambda e: (NORM(*e), abs(e[1]), e[1], e[0]))
    return out


def divisors(n):
    ds = [1]
    m, p = n, 2
    while p * p <= m:
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p; e += 1
            ds = [d * p ** k for d in ds for k in range(e + 1)]
        p += 1
    if m > 1:
        ds = ds + [d * m for d in ds]
    return sorted(ds)


class Atoms:
    """Atoms sorted by arrival; factorization enumeration over them."""

    def __init__(self, atoms):
        self.list = []
        self.idx = {}
        self.by_norm = {}
        self._memo = {}
        for a in atoms:
            self.add(a)

    def add(self, a):
        self.idx[a] = len(self.list)
        self.list.append(a)
        self.by_norm.setdefault(NORM(*a), []).append(a)
        self._memo = {}

    def copy(self):
        return Atoms(list(self.list))

    def facts(self, x, start=0, cap=10 ** 6):
        """All factorizations of x (canonical) into atoms with index >= start,
        as tuples of atom indices in non-decreasing order."""
        key = (x, start)
        if key in self._memo:
            return self._memo[key]
        n = NORM(*x)
        res = []
        if n == 1:
            res = [()]
        else:
            i = self.idx.get(x, -1)
            if i >= start:
                res.append((i,))
            for d in divisors(n):
                if d == 1 or d * d > n:
                    continue
                for at in self.by_norm.get(d, ()):
                    j = self.idx[at]
                    if j < start:
                        continue
                    q = div(x, at)
                    if q is None:
                        continue
                    for f in self.facts(canon(*q), j, cap):
                        res.append((j,) + f)
                        if len(res) >= cap:
                            break
        self._memo[key] = res
        return res


# ------------------------------------------------------------ ideal language
def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def sqrt_mod(a, p):
    a %= p
    for r in range(p):
        if r * r % p == a:
            return r
    return None


def vp(n, p):
    e = 0
    while n and n % p == 0:
        n //= p; e += 1
    return e


def rational_primes(n):
    out, m, p = [], n, 2
    while p * p <= m:
        if m % p == 0:
            out.append(p)
            while m % p == 0:
                m //= p
        p += 1
    if m > 1:
        out.append(m)
    return out


_root = {}


def ideal_name(p, r=None):
    if p == 2:
        return '𝔭₂'
    if p == 5:
        return '𝔭₅'
    if r is None:
        return f'({p})'
    return f'𝔭{p},{r}'


def valuation(x):
    """Prime-ideal factorization of the principal ideal (x), as {name: exponent}."""
    a, b = x
    n = NORM(a, b)
    out = {}
    for p in rational_primes(n):
        if p == 2 or p == 5:
            out[ideal_name(p)] = vp(n, p)
            continue
        if legendre(-5, p) == -1:
            out[ideal_name(p)] = min(vp(a, p) if a else 99, vp(b, p) if b else 99)
            continue
        if p not in _root:
            r = sqrt_mod(-5, p)
            _root[p] = min(r, p - r)
        r = _root[p]
        k = min(vp(a, p) if a else 99, vp(b, p) if b else 99)
        a1, b1 = a // p ** k, b // p ** k
        e = vp(NORM(a1, b1), p)
        P, Q = ideal_name(p, r), ideal_name(p, p - r)
        vP = vQ = k
        if e:
            # sqrt(-5) -> -r in the map Z[sqrt(-5)] -> Z/p killing (p, sqrt(-5) + r) ... choose one side
            if (a1 + b1 * (p - r)) % p == 0:
                vP += e
            else:
                vQ += e
        if vP:
            out[P] = vP
        if vQ:
            out[Q] = vQ
    return out


def addv(u, v):
    w = dict(u)
    for k, e in v.items():
        w[k] = w.get(k, 0) + e
    return w


class ElementRule:
    kind = 'element'

    def __init__(self, atoms):
        self.atoms = atoms


class IdealRule:
    kind = 'ideal'

    def __init__(self, atoms):
        self.atoms = atoms          # inherited element atoms (the old record)
        self.vals = {}


class RingOasis(Environment):
    name = 'Oasis of Z[√−5]'

    def __init__(self, X):
        self.X = X
        self.els = elements(X)
        self.first_contradiction = None
        self.ideal_checks = 0

    def founders(self):
        return self.els[:10]

    def newcomer(self, k):
        return self.els[10 + k] if 10 + k < len(self.els) else None

    def describe(self, rule):
        if rule.kind == 'element':
            return f'{len(rule.atoms.list)} atoms; every number a product of atoms in one way'
        return 'prime ideals; every number a product of prime ideals in one way'

    def constitution(self, record):
        at = Atoms([])
        for x in record:
            if not at.facts(x):
                at.add(x)
        return ElementRule(at)

    def test(self, rule, item, record):
        fs = rule.atoms.facts(item)
        if rule.kind == 'element':
            if not fs:
                return 'absurd'
            return 'ok' if len(fs) == 1 else 'contradiction'
        # ideal language: every element factorization must give the same ideal record
        if not fs:
            rule.atoms.add(item)              # a new irreducible joins the inherited atoms
            fs = rule.atoms.facts(item)
        v = valuation(item)
        for f in fs:
            self.ideal_checks += 1
            s = {}
            for j in f:
                s = addv(s, self._val(rule, rule.atoms.list[j]))
            if s != v:
                return 'contradiction'
        return 'ok'

    def _val(self, rule, atom):
        if atom not in rule.vals:
            rule.vals[atom] = valuation(atom)
        return rule.vals[atom]

    def revise(self, rule, record, backlog):
        if rule.kind == 'ideal':
            return None, 'contradiction'
        new = ElementRule(rule.atoms.copy())
        for x in backlog:
            if not new.atoms.facts(x):
                new.atoms.add(x)
        if all(len(new.atoms.facts(x)) == 1 for x in backlog):
            return new, 'logical'
        # Every irreducible must be an atom (nothing else can explain it), so the atom
        # list is forced. Uniqueness fails for the backlog: no rule in this language passes C2.
        if self.first_contradiction is None:
            x = [y for y in backlog if len(new.atoms.facts(y)) > 1][0]
            self.first_contradiction = (x, [[new.atoms.list[j] for j in f] for f in new.atoms.facts(x)])
        return None, 'contradiction'

    def ascend(self, rule, record, backlog):
        return IdealRule(rule.atoms.copy()), 'ideal language (Kummer 1847, Dedekind 1871)'

    def verify(self, old, new, record, backlog):
        if new.kind == 'element':
            c1 = set(old.atoms.list) <= set(new.atoms.list)
            added = set(new.atoms.list) - set(old.atoms.list)
            c2 = all(len(new.atoms.facts(x)) == 1 for x in record
                     if any(NORM(*x) % NORM(*a) == 0 for a in added))
            c3 = all(len(new.atoms.facts(x)) == 1 for x in backlog)
            return c1, c2, c3
        # ideal language: C1 every old explanation maps to an ideal explanation;
        # C2 every record entry gets one ideal record; C3 the contradiction disappears
        def one_record(x):
            v = valuation(x)
            fs = new.atoms.facts(x)
            return all(self._sumv(new, f) == v for f in fs)
        c1 = all(len(new.atoms.facts(x)) >= 1 for x in record)
        c2 = all(one_record(x) for x in record)
        c3 = all(one_record(x) for x in backlog)
        return c1, c2, c3

    def _sumv(self, rule, f):
        s = {}
        for j in f:
            s = addv(s, self._val(rule, rule.atoms.list[j]))
        return s


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    t = time.time()
    env = RingOasis(X)
    res = run(env, 0, 10 ** 9)
    ev = res['events']
    out = dict(X=X, elements=len(env.els), days=res['last_day'])
    out['events_summary'] = {}
    for e in ev:
        out['events_summary'][e.kind] = out['events_summary'].get(e.kind, 0) + 1
    special = [e.as_dict() for e in ev if e.kind not in ('logical', 'constitution')]
    out['special_events'] = special
    fc = env.first_contradiction
    out['first_contradiction'] = dict(element=fmt(fc[0]), norm=NORM(*fc[0]),
                                      factorizations=[[fmt(a) for a in f] for f in fc[1]],
                                      day=[e.day for e in ev if e.kind == 'contradiction'][0])
    out['founders'] = [fmt(x) for x in env.founders()]
    out['ideal_checks'] = env.ideal_checks
    out['contradictions_after_ascent'] = sum(1 for e in ev if e.kind == 'contradiction') - 1
    # ---------------- statistics over all elements in the record
    rule = res['rule']
    at = rule.atoms
    buckets = [(2, 100), (101, 1000), (1001, 3000), (3001, 10000), (10001, 30000), (30001, 100000)]
    stats = []
    lengths_ok = True
    maxf = (0, None)
    nonunique_examples = []
    total_f = 0
    ideal_disagree = 0
    for lo, hi in buckets:
        if lo > X:
            break
        red = nonu = 0
        for x in env.els:
            n = NORM(*x)
            if n < lo or n > min(hi, X):
                continue
            fs = at.facts(x)
            if len(fs) == 1 and len(fs[0]) == 1:
                continue
            red += 1
            total_f += len(fs)
            if len(fs) > 1:
                nonu += 1
                if len(nonunique_examples) < 6 and n <= 200:
                    nonunique_examples.append((fmt(x), [[fmt(at.list[j]) for j in f] for f in fs]))
            if len({len(f) for f in fs}) != 1:
                lengths_ok = False
            if len(fs) > maxf[0]:
                maxf = (len(fs), fmt(x), n)
            v = valuation(x)
            for f in fs:
                s = {}
                for j in f:
                    s = addv(s, valuation(at.list[j]))
                if s != v:
                    ideal_disagree += 1
        stats.append(dict(lo=lo, hi=min(hi, X), reducible=red, nonunique=nonu,
                          share=round(nonu / red, 4) if red else None))
    # atoms that are prime ideals vs atoms that split into two non-principal prime ideals
    prime_atoms = sum(1 for a in at.list if sum(valuation(a).values()) == 1)
    split_atoms = sum(1 for a in at.list if sum(valuation(a).values()) == 2)
    out.update(stats=stats, all_lengths_equal=lengths_ok, max_factorizations=maxf,
               examples=nonunique_examples, factorization_records=total_f,
               ideal_disagreements=ideal_disagree, atoms=len(at.list),
               atoms_prime_ideal=prime_atoms, atoms_two_ideals=split_atoms,
               six_ideals={fmt(a): valuation(a) for a in [(2, 0), (3, 0), (1, 1), (1, -1), (6, 0)]},
               seconds=round(time.time() - t, 1))
    json.dump(out, open('ring.json', 'w'), indent=1, ensure_ascii=False)
    print(json.dumps({k: v for k, v in out.items() if k not in ('special_events',)}, indent=1, ensure_ascii=False)[:4000])


if __name__ == '__main__':
    main()
