"""Exact replay of the witness states in proof.md; a finite check, not the proof.

Each state is rebuilt from its gate list, the stated step is performed with
literal and constant propagation, and the change of every potential term is
recorded. It also checks the syntactic conditions showing that no earlier row
applies, and then, for every choice the variant may make at S1 and S2, the
certificate bound in exact rational arithmetic. It proves nothing about
other states.
"""
import json
import sys
from fractions import Fraction

AND = (0, 0, 0, 1)
XOR = (0, 1, 1, 0)


class State:
    def __init__(self, params, gates, output, pairs=(), equations=()):
        self.params = list(params)
        self.gates = {name: [table, a, b] for name, table, a, b in gates}
        self.order = [name for name, _, _, _ in gates]
        self.output = output
        self.pairs = {frozenset(p): dict(zip(p, corner)) for p, corner in pairs}
        self.equations = {name: set(terms) for name, terms in equations}
        self.aux = {t for _, terms in equations for t in terms} - set(self.params)
        self.deleted = 0

    def partner(self, u):
        for pair in self.pairs:
            if u in pair:
                return next(iter(pair - {u}))
        return None

    def degree(self, u):
        readers = {g for g, (_, a, b) in self.gates.items() if u in (a, b)}
        return len(readers) + sum(u in terms for terms in self.equations.values())

    def terms(self):
        """Potential terms: s, unpaired degree classes, pair degree profiles."""
        out = {"s": len(self.gates) + len(self.equations)}
        for u in self.params:
            d = min(self.degree(u), 3)
            if self.partner(u) is None and d:
                out[f"k{d}"] = out.get(f"k{d}", 0) + 1
        for pair in self.pairs:
            p, r = sorted(pair)
            prof = "".join(sorted(str(min(self.degree(v), 3)) for v in (p, r)))
            out[f"pi{prof}"] = out.get(f"pi{prof}", 0) + 1
        return out

    def a(self):
        return len(self.params) - len(self.pairs)

    def replace(self, old, new, negate=False):
        """Reroute every use of old to the constant or literal new."""
        for gate in self.gates.values():
            for pin in (1, 2):
                if gate[pin] == old:
                    gate[pin] = new
                    if negate:
                        t = list(gate[0])
                        gate[0] = tuple(t[i ^ (2 if pin == 1 else 1)] for i in range(4))
        if self.output == old:
            self.output = new

    def propagate(self):
        """Delete gates with a constant or repeated pin, rerouting to a literal."""
        changed = True
        while changed:
            changed = False
            for g in list(self.order):
                if g not in self.gates:
                    continue
                table, a, b = self.gates[g]
                if a in (0, 1) or b in (0, 1) or a == b:
                    if a in (0, 1) and b in (0, 1):
                        f = (table[2 * a + b],) * 2
                        other = None
                    elif a in (0, 1):
                        f = (table[2 * a], table[2 * a + 1])
                        other = b
                    elif b in (0, 1):
                        f = (table[b], table[2 + b])
                        other = a
                    else:
                        f = (table[0], table[3])
                        other = a
                    del self.gates[g]
                    self.order.remove(g)
                    self.deleted += 1
                    if f[0] == f[1]:
                        self.replace(g, f[0])
                    else:
                        self.replace(g, other, negate=(f == (1, 0)))
                    changed = True
            for g in list(self.order):
                if g != self.output and not any(g in (a, b) for _, a, b in self.gates.values()):
                    del self.gates[g]
                    self.order.remove(g)
                    self.deleted += 1
                    changed = True

    def fix(self, u, value):
        assert not any(u in terms for terms in self.equations.values())
        self.params.remove(u)
        for pair in [p for p in self.pairs if u in p]:
            del self.pairs[pair]
        self.replace(u, value)
        self.propagate()

    def pair_gate(self, g):
        """GP: pair the parents of AND-type g on its three majority points."""
        table, x, y = self.gates[g]
        corner = next(i for i in range(4) if table.count(table[i]) == 1)
        self.pairs[frozenset((x, y))] = {x: corner >> 1, y: corner & 1}
        majority = table[corner ^ 3]
        del self.gates[g]
        self.order.remove(g)
        self.deleted += 1
        self.replace(g, majority)
        self.propagate()

    def safe_value(self, u):
        pair = next(p for p in self.pairs if u in p)
        return 1 - self.pairs[pair][u]

    def move_affine(self, g, aux):
        """M: an XOR-type gate on two parameters or auxiliaries becomes an equation."""
        table, a, b = self.gates.pop(g)
        self.order.remove(g)
        assert table in (XOR, (1, 0, 0, 1)) and {a, b} <= set(self.params) | self.aux
        self.equations["E" + g] = {aux, a, b}
        self.aux.add(aux)
        self.replace(g, aux)

    def d4_exchange(self, p, eq, q):
        """D4: unpaired p is read only by eq; the auxiliary q in eq becomes a parameter.

        The new parameter keeps the name q (proof.md calls it w in S3).
        """
        assert self.partner(p) is None and self.degree(p) == 1
        assert p in self.equations[eq] and q in self.equations[eq] and q in self.aux
        assert not any(p in (x, y) for _, x, y in self.gates.values())
        self.params.remove(p)
        self.params.append(q)
        self.aux.discard(q)
        del self.equations[eq]
        self.deleted += 1


def delta(before, after):
    keys = set(before) | set(after)
    return {k: after.get(k, 0) - before.get(k, 0) for k in keys if after.get(k, 0) != before.get(k, 0)}


def step(build, action):
    st = build()
    t0, a0 = st.terms(), st.a()
    action(st)
    return st, delta(t0, st.terms()), st.a() - a0


def opaque(names):
    gates, params = [], []
    for i, h in enumerate(names):
        params += [f"t{2 * i + 1}", f"t{2 * i + 2}"]
        gates.append((h, AND, f"t{2 * i + 1}", f"t{2 * i + 2}"))
    return params, gates


def s1():
    tp, tg = opaque(["h1", "h2", "h3", "h4", "h5"])
    gates = [("G", AND, "x", "y")] + tg + [
        ("A", XOR, "x", "h1"), ("B", XOR, "x", "h2"), ("C", XOR, "y", "h3"),
        ("D", XOR, "y", "h4"), ("J", XOR, "G", "h5"), ("O1", AND, "A", "B"),
        ("O2", AND, "C", "D"), ("O3", AND, "O1", "O2"), ("OUT", AND, "O3", "J")]
    return State(["x", "y"] + tp, gates, "OUT")


def s2():
    tp, tg = opaque(["h1", "h3", "h6", "h7", "h8"])
    gates = [("G", AND, "x", "y")] + tg + [
        ("A", XOR, "x", "h1"), ("C", XOR, "y", "h3"), ("J", XOR, "G", "z"),
        ("J1", AND, "J", "h6"), ("J2", AND, "J", "h7"), ("J3", AND, "J", "h8"),
        ("O1", AND, "A", "C"), ("O2", AND, "J1", "J2"), ("O3", AND, "O2", "J3"),
        ("OUT", AND, "O1", "O3")]
    return State(["x", "y", "z"] + tp, gates, "OUT")


def s3():
    tp, tg = opaque(["ha", "hb", "hc", "hd", "he", "hf", "hg"])
    gates = tg + [
        ("H", XOR, "r", "ha"), ("V1", AND, "v", "hb"), ("V2", AND, "v", "hc"),
        ("V3", AND, "v", "hd"), ("Z1", AND, "z", "he"), ("Z2", AND, "z", "hf"),
        ("Z3", AND, "z", "hg"), ("O1", AND, "H", "V1"), ("O2", AND, "V2", "V3"),
        ("O3", AND, "Z1", "Z2"), ("O4", AND, "O1", "O2"), ("O5", AND, "O3", "Z3"),
        ("OUT", AND, "O4", "O5")]
    return State(["p", "r", "z"] + tp, gates, "OUT",
                 pairs=[(("p", "r"), (1, 1))], equations=[("E", ["v", "p", "z"])])


def m_applies(st):
    affine = set(st.params) | st.aux
    return any(t in (XOR, (1, 0, 0, 1)) and {x, y} <= affine for t, x, y in st.gates.values())


def no_deletion_or_m(st):
    """D1-D4 and M do not apply (syntactic checks of proof.md Section 2)."""
    affine = set(st.params) | st.aux
    readers = lambda w: [g for g, (_, x, y) in st.gates.items() if w in (x, y)]
    for g, (table, x, y) in st.gates.items():
        assert x not in (0, 1) and y not in (0, 1) and x != y, g
        assert table in (AND, XOR), g
        if table == XOR:
            assert not {x, y} <= affine, g
        if {x, y} <= affine:
            assert {x, y} <= set(st.params), g
        assert g == st.output or readers(g), g
    for e, terms in st.equations.items():
        assert len(terms) == 3, e
    for u in st.params:
        if st.partner(u) is None and not readers(u):
            assert sum(u in t for t in st.equations.values()) != 1, u


def run():
    a, b, c = s1(), s2(), s3()
    assert [a.degree(u) for u in "xy"] == [3, 3] and a.order[0] == "G"
    assert [b.degree(u) for u in "xyz"] == [2, 2, 1] and b.order[0] == "G"
    assert [c.degree(u) for u in ("p", "r", "z")] == [1, 1, 4] and c.degree("v") == 4
    for st in (a, b, c):
        no_deletion_or_m(st)
    assert not a.pairs and not b.pairs
    # S3: P0, P2 fail (degrees 1), PC fails (different readers), PE applies.
    assert c.degree("p") == c.degree("r") == 1
    assert [e for e, t in c.equations.items() if "p" in t] == ["E"]
    assert not any("p" in (x, y) for _, x, y in c.gates.values())
    assert [g for g, (_, x, y) in c.gates.items() if "r" in (x, y)] == ["H"]
    # S1': no deletion, M or P0; x and y are paired with degree 2, so P2 applies.
    s1p = s1()
    s1p.pair_gate("G")
    no_deletion_or_m(s1p)
    assert s1p.degree("x") == s1p.degree("y") == 2

    def gp(st):
        st.pair_gate("G")

    def after_gp(build):
        def make():
            st = build()
            st.pair_gate("G")
            return st
        return make

    measured = {
        "S1 GP": step(s1, gp),
        "S1' P2 on x": step(after_gp(s1), lambda st: st.fix("x", st.safe_value("x"))),
        "S1' P2 on y": step(after_gp(s1), lambda st: st.fix("y", st.safe_value("y"))),
        "S2 GP": step(s2, gp),
        "S3 PE": step(s3, lambda st: (st.fix("r", st.safe_value("r")), st.d4_exchange("p", "E", "v"))),
    }
    for name, build in (("S1", s1), ("S2", s2)):
        for u in "xy":
            for v in (0, 1):
                measured[f"{name} fix {u}={v}"] = step(build, lambda st, u=u, v=v: st.fix(u, v))

    for u in "xy":
        after = s1()
        after.fix(u, 1)
        assert not m_applies(after)

    def fix_then_m_d4(u):
        def act(st):
            st.fix(u, 1)
            assert m_applies(st)
            st.move_affine("J", "vJ")
            other = "y" if u == "x" else "x"
            assert st.degree("z") == 1 and st.degree(other) == 2
            st.d4_exchange("z", "EJ", "vJ")
        return act
    for u in "xy":
        measured[f"S2 fix {u}=1 with next M and D4"] = step(s2, fix_then_m_d4(u))
    expected = {
        "S1 GP": ({"s": -2, "k3": -2, "pi22": 1}, -1),
        "S1' P2 on x": ({"s": -2, "k2": 1, "pi22": -1}, 0),
        "S1' P2 on y": ({"s": -2, "k2": 1, "pi22": -1}, 0),
        "S2 GP": ({"s": -2, "k1": -1, "k2": -2, "k3": 1, "pi11": 1}, -1),
        "S3 PE": ({"s": -2, "pi11": -1, "k3": 1}, 0),
    }
    for u in "xy":
        expected[f"S1 fix {u}=0"] = ({"s": -4, "k3": -2, "k2": 1}, -1)
        expected[f"S1 fix {u}=1"] = ({"s": -3, "k3": -1}, -1)
        expected[f"S2 fix {u}=0"] = ({"s": -3, "k2": -2, "k3": 1}, -1)
        expected[f"S2 fix {u}=1"] = ({"s": -2, "k2": -1}, -1)
        expected[f"S2 fix {u}=1 with next M and D4"] = ({"s": -3, "k2": -2, "k3": 1}, -1)
    rows = {}
    for name, (_, d, da) in measured.items():
        assert (d, da) == expected[name], (name, d, da)
        rows[name] = {"dW_terms": d, "da": da}

    # Row r requires dW <= delta*da. Write delta = c + K with K = max(k1, k2, k3);
    # each row becomes a linear inequality  lhs . v <= rhs  over v = (c, K, k1, k2, k3, pi11, pi22).
    def ineq(*names):
        lhs, da_total = {}, 0
        for n in names:
            d, da = expected[n]
            da_total += da
            for t, coef in d.items():
                if t != "s":
                    lhs[t] = lhs.get(t, 0) + coef
        rhs = -sum(expected[n][0]["s"] for n in names)
        if da_total:
            lhs["c"] = lhs.get("c", 0) - da_total
            lhs["K"] = lhs.get("K", 0) - da_total
        return lhs, Fraction(rhs)

    def kmax(t):
        return {t: 1, "K": -1}, Fraction(0)

    def combine(parts):
        lhs, rhs = {}, Fraction(0)
        for lam, (l, r) in parts:
            for t, coef in l.items():
                lhs[t] = lhs.get(t, 0) + Fraction(lam) * coef
            rhs += Fraction(lam) * r
        lhs = {t: v for t, v in lhs.items() if v}
        assert set(lhs) == {"c"} and lhs["c"] > 0, lhs
        return rhs / lhs["c"]

    s1_rows = {"GP": ("S1 GP", "S1' P2 on x"), "fix0": ("S1 fix x=0",), "fix1": ("S1 fix x=1",)}
    bounds = {}
    for c1, r1 in s1_rows.items():
        for c2 in ("GP", "fix0", "fix1", "fix1+M+D4"):
            if c1 == "fix1":
                b = combine([(1, ineq(*r1)), (1, kmax("k3"))])
            elif c2 == "GP":
                b = combine([(2, ineq(*r1)), (1, ineq("S2 GP")), (1, ineq("S3 PE")),
                             (2, kmax("k3")), (1, kmax("k1"))])
            elif c2 == "fix0":
                b = combine([(2, ineq(*r1)), (1, ineq("S2 fix x=0")), (3, kmax("k3"))])
            elif c2 == "fix1+M+D4":
                b = combine([(2, ineq(*r1)), (1, ineq("S2 fix x=1 with next M and D4")), (3, kmax("k3"))])
            else:
                b = combine([(1, ineq("S2 fix x=1")), (1, kmax("k2"))])
            bounds[f"S1 {c1} / S2 {c2}"] = str(b)
    assert bounds["S1 GP / S2 GP"] == "4"
    assert max(Fraction(v) for v in bounds.values()) == 4
    # Flat weights (5, 5, 5; pairs 3) with delta = 9 meet every witness row of the rule.
    flat = {"k1": 5, "k2": 5, "k3": 5, "pi11": 3, "pi22": 3}
    for n in ("S1 GP", "S1' P2 on x", "S2 GP", "S3 PE"):
        d, da = expected[n]
        assert d["s"] + sum(coef * flat[t] for t, coef in d.items() if t != "s") == 9 * da, n
    return {"all_passed": True, "rows": rows, "coefficient_bounds_by_choice": bounds,
            "scope": "Replays only the explicit witness states of proof.md; the inequalities and their consequences are proved there."}


if __name__ == "__main__":
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError("the declared audit takes an empty input object")
    print(json.dumps(run(), sort_keys=True))
