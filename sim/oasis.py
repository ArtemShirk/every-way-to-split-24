"""Oasis-class engine.

One protocol for every desert. An environment supplies the inhabitants and the
language its constitutions are written in; the engine runs the days, counts
absurd potential, triggers collapse when P > C, and checks every revision
against the three criteria:

  C1 inheritance   the new constitution still explains everything the old one did
  C2 record        the new constitution agrees with every entry in the record
  C3 resolution    the new constitution explains the newcomers that caused the collapse

Protocol (fixed):
  day 1        10 founders arrive
  day 10       constitution written from the founders' record
  days 11-20   quiet
  day 21+      one newcomer per day

Outcomes of a collapse:
  logical        a revision in the current language passes C1, C2 and C3
  contradiction  the record holds two incompatible entries; no revision in the
                 language can pass C2
  exhausted      the language has no rule left that passes C2 and C3
In the last two cases the engine asks the environment to ascend to a parent
language. If there is none, the run ends there.
"""

FOUNDERS_DAY, CONSTITUTION_DAY, QUIET_END, FIRST_NEWCOMER = 1, 10, 20, 21


class Environment:
    """Interface. Subclasses fill these in."""
    name = 'environment'
    language = 'language'

    def founders(self):
        """Ten inhabitants present on day 1."""
        raise NotImplementedError

    def newcomer(self, k):
        """The k-th newcomer (k = 0 on day 21)."""
        raise NotImplementedError

    def constitution(self, record):
        """Rule written on day 10 from the founders' record."""
        raise NotImplementedError

    def test(self, rule, item, record):
        """Return 'ok', 'absurd' (rule cannot explain it) or 'contradiction'
        (rule explains it two incompatible ways). May extend the rule in place
        when the item falls where the rule makes no claim yet."""
        raise NotImplementedError

    def revise(self, rule, record, backlog):
        """Return (new_rule, outcome). outcome is 'logical', 'contradiction' or
        'exhausted'; new_rule is None unless outcome is 'logical'."""
        raise NotImplementedError

    def verify(self, old, new, record, backlog):
        """Return (c1, c2, c3) as booleans."""
        raise NotImplementedError

    def describe(self, rule):
        return str(rule)

    def ascend(self, rule, record, backlog):
        """Parent language. Return (new_env_state_rule, note) or None."""
        return None


class Event:
    __slots__ = ('day', 'kind', 'items', 'note', 'checks')

    def __init__(self, day, kind, items, note='', checks=None):
        self.day, self.kind, self.items, self.note, self.checks = day, kind, items, note, checks

    def as_dict(self):
        return dict(day=self.day, kind=self.kind, items=[str(i) for i in self.items],
                    note=self.note, checks=self.checks)


def run(env, capacity, max_day, on_day=None):
    """Run the protocol up to max_day. Returns dict with record, rule, events, stats."""
    record = list(env.founders())
    rule = env.constitution(record)
    events = [Event(CONSTITUTION_DAY, 'constitution', [], note=env.describe(rule))]
    P, backlog = 0, []
    absurd_total = 0
    for day in range(FIRST_NEWCOMER, max_day + 1):
        item = env.newcomer(day - FIRST_NEWCOMER)
        if item is None:
            break
        verdict = env.test(rule, item, record)
        record.append(item)
        if verdict != 'ok':
            P += 1
            absurd_total += 1
            backlog.append(item)
        if on_day:
            on_day(day, item, verdict, P)
        if P > capacity:
            new, outcome = env.revise(rule, record, backlog)
            if outcome == 'logical':
                checks = env.verify(rule, new, record, backlog)
                kind = 'logical' if all(checks) else 'illogical'
                events.append(Event(day, kind, list(backlog), env.describe(new), checks))
                rule = new
            else:
                events.append(Event(day, outcome, list(backlog), env.describe(rule)))
                up = env.ascend(rule, record, backlog)
                if up is None:
                    events.append(Event(day, 'end', [], 'no parent language'))
                    break
                new, note = up
                checks = env.verify(rule, new, record, backlog)
                events.append(Event(day, 'ascent', list(backlog), note, checks))
                rule = new
            P, backlog = 0, []
    return dict(record=record, rule=rule, events=events, absurd_total=absurd_total,
                last_day=day)
