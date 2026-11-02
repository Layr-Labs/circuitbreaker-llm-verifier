"""Executable audit of the ternary-pair gate elimination (elimination.md).

Runs the elimination procedure on random full-B2 circuits and checks, after
every row:
  * semantics: the state computes f(Phi(u)) on every point of Dom,
  * Phi is injective and the affine system is uniquely solvable,
  * every equation has at most three terms and no gate feeds an equation,
  * Delta W <= 9 Delta a with Delta a in {0,-1}, W = s + 5I + 3q,
  * each row's measured Delta W is at most its tabulated bound.
The loop runs until an exclusion of Section 1 would be needed (affine or
constant output, or the minimal AND is the output); for an affine disperser
those states cannot occur while a >= d+2. Random circuits are not dispersers,
so this tests the local accounting and semantics, not the disperser premise.
Usage: python audit.py [trials] [seed]
"""
import random
import sys

XOR_TABLES = (0b0110, 0b1001)


def tbl_val(t, a, b):
    return (t >> (2 * a + b)) & 1


def is_and_type(t):
    return bin(t).count("1") in (1, 3)


class Fail(Exception):
    pass


class State:
    def __init__(self, n, gates, out):
        self.n = n
        self.orig = (list(gates), out)
        self.nid = n + len(gates)
        self.kind = {}
        self.params = set(range(n))
        for i in range(n):
            self.kind[i] = "param"
        self.phi = [(1 << i, 0) for i in range(n)]
        self.partner, self.excl = {}, {}
        self.auxes = set()
        self.eqs, self.neq = {}, 0
        self.gates, self.gorder = {}, []
        for j, (t, a, b) in enumerate(gates):
            g = n + j
            self.kind[g] = "gate"
            self.gates[g] = [t, ("v", a), ("v", b)]
            self.gorder.append(g)
        self.out = ("v", out, 0)

    # ---------- semantics ----------
    def orig_eval(self, x):
        gates, out = self.orig
        val = list(x)
        for t, a, b in gates:
            val.append(tbl_val(t, val[a], val[b]))
        return val[out]

    def forms(self):
        """Solve the affine system: aux -> (param mask, const)."""
        auxl = sorted(self.auxes)
        if len(self.eqs) != len(auxl):
            raise Fail("equation/aux count mismatch")
        idx = {v: i for i, v in enumerate(auxl)}
        rows = []
        for terms, c in self.eqs.values():
            am, pm = 0, 0
            for v in terms:
                if self.kind[v] == "aux":
                    am ^= 1 << idx[v]
                elif self.kind[v] == "param":
                    pm ^= 1 << v
                else:
                    raise Fail("gate feeds an equation")
            rows.append([am, pm, c])
        piv = []
        r = 0
        for col in range(len(auxl)):
            sel = next((i for i in range(r, len(rows)) if rows[i][0] >> col & 1), None)
            if sel is None:
                raise Fail("affine system singular")
            rows[r], rows[sel] = rows[sel], rows[r]
            for i in range(len(rows)):
                if i != r and rows[i][0] >> col & 1:
                    rows[i] = [rows[i][0] ^ rows[r][0], rows[i][1] ^ rows[r][1], rows[i][2] ^ rows[r][2]]
            piv.append(col)
            r += 1
        return {auxl[c]: (rows[i][1], rows[i][2]) for i, c in enumerate(piv)}

    def check(self):
        F = self.forms()
        P = sorted(self.params)
        # Phi injective: rank of its columns equals m
        basis = []
        for p in P:
            col = 0
            for i, (mk, _) in enumerate(self.phi):
                if mk >> p & 1:
                    col |= 1 << i
            for bv in basis:
                col = min(col, col ^ bv)
            if col == 0:
                raise Fail("Phi not injective")
            basis.append(col)
        for terms, _ in self.eqs.values():
            if len(terms) > 3:
                raise Fail("equation with more than three terms")
        for bits in range(1 << len(P)):
            u = {p: bits >> i & 1 for i, p in enumerate(P)}
            if any(u[p] == self.excl[p] and u[self.partner[p]] == self.excl[self.partner[p]]
                   for p in self.partner):
                continue
            x = [(bin(mk & sum(1 << p for p in P if u[p])).count("1") + c) & 1 for mk, c in self.phi]
            val = dict(u)
            for v, (mk, c) in F.items():
                val[v] = (bin(mk & sum(1 << p for p in P if u[p])).count("1") + c) & 1
            for g in self.gorder:
                t, A, B = self.gates[g]
                val[g] = tbl_val(t, self.pin(A, val), self.pin(B, val))
            o = self.out[1] if self.out[0] == "c" else val[self.out[1]] ^ self.out[2]
            if o != self.orig_eval(x):
                raise Fail("semantics broken on Dom")

    @staticmethod
    def pin(p, val):
        return p[1] if p[0] == "c" else val[p[1]]

    # ---------- measures ----------
    def children(self, v):
        ch = set()
        for e, (terms, _) in self.eqs.items():
            if v in terms:
                ch.add(("e", e))
        for g, (t, A, B) in self.gates.items():
            if A == ("v", v) or B == ("v", v):
                ch.add(("g", g))
        return ch

    def W(self):
        s = len(self.eqs) + len(self.gates)
        q = len(self.partner) // 2
        I = sum(1 for p in self.params if p not in self.partner and self.children(p))
        return s + 5 * I + 3 * q

    def a(self):
        return len(self.params) - len(self.partner) // 2

    # ---------- primitive edits ----------
    def sub_value(self, v, lit):
        """Replace every use of value v by lit = ('c',b) or ('v',y,neg)."""
        for terms_c in self.eqs.values():
            terms = terms_c[0]
            if v in terms:
                terms.discard(v)
                if lit[0] == "c":
                    terms_c[1] ^= lit[1]
                else:
                    terms_c[1] ^= lit[2]
                    terms ^= {lit[1]}
        for g, gt in self.gates.items():
            for k in (1, 2):
                if gt[k] == ("v", v):
                    if lit[0] == "c":
                        gt[k] = ("c", lit[1])
                    else:
                        gt[k] = ("v", lit[1])
                        if lit[2]:
                            t = gt[0]
                            gt[0] = sum(tbl_val(t, a ^ (k == 1), b ^ (k == 2)) << (2 * a + b)
                                        for a in (0, 1) for b in (0, 1))
        if self.out[0] == "v" and self.out[1] == v:
            self.out = ("c", lit[1] ^ self.out[2]) if lit[0] == "c" else ("v", lit[1], lit[2] ^ self.out[2])

    def fix_param(self, x, lit):
        """x := constant or literal of another param (affine restriction)."""
        mk, c = (0, lit[1]) if lit[0] == "c" else (1 << lit[1], lit[2])
        self.phi = [((m ^ (1 << x) ^ mk, cc ^ c) if m >> x & 1 else (m, cc)) for m, cc in self.phi]
        self.params.discard(x)
        self.sub_value(x, lit)

    def unpair(self, p):
        r = self.partner.pop(p)
        self.partner.pop(r)
        self.excl.pop(p)
        self.excl.pop(r)
        return r

    def exchange(self, x, Q, F):
        mk, c = F[Q]
        assert mk >> x & 1 and x not in self.partner
        ell = mk ^ (1 << x)
        self.phi = [((m ^ (1 << x) ^ (1 << Q) ^ ell, cc ^ c) if m >> x & 1 else (m, cc)) for m, cc in self.phi]
        self.params.discard(x)
        self.auxes.discard(Q)
        self.params.add(Q)
        self.auxes.add(x)
        self.kind[x], self.kind[Q] = "aux", "param"

    def del_gate(self, g):
        del self.gates[g]
        self.gorder.remove(g)

    # ---------- deletions D1-D4, M ----------
    def gate_function(self, g, F):
        t, A, B = self.gates[g]
        # dependent affine parents (and constant-valued affine parents)
        for k, P in ((1, A), (2, B)):
            if P[0] == "v" and self.kind[P[1]] in ("param", "aux"):
                mk, c = F[P[1]] if self.kind[P[1]] == "aux" else (1 << P[1], 0)
                if mk == 0:
                    self.gates[g][k] = ("c", c)
        t, A, B = self.gates[g]
        if A[0] == "v" and B[0] == "v" and A[1] != B[1] and all(self.kind[P[1]] in ("param", "aux") for P in (A, B)):
            fa = F[A[1]] if self.kind[A[1]] == "aux" else (1 << A[1], 0)
            fb = F[B[1]] if self.kind[B[1]] == "aux" else (1 << B[1], 0)
            if fa[0] == fb[0]:
                neg = fa[1] ^ fb[1]
                self.gates[g][2] = ("v", A[1])
                if neg:
                    self.gates[g][0] = sum(tbl_val(t, a, b ^ 1) << (2 * a + b) for a in (0, 1) for b in (0, 1))
        t, A, B = self.gates[g]
        ins = sorted({P[1] for P in (A, B) if P[0] == "v"})
        outs = {}
        for bits in range(1 << len(ins)):
            env = {v: bits >> i & 1 for i, v in enumerate(ins)}
            av = A[1] if A[0] == "c" else env[A[1]]
            bv = B[1] if B[0] == "c" else env[B[1]]
            outs[bits] = tbl_val(t, av, bv)
        vals = set(outs.values())
        if len(vals) == 1:
            return ("c", vals.pop())
        for i, v in enumerate(ins):
            others = [j for j in range(len(ins)) if j != i]
            if all(outs[bits] == outs[bits & ~(sum(1 << j for j in others))] for bits in outs):
                if all(outs[bits] == outs[bits ^ (1 << j)] for bits in outs for j in others):
                    return ("v", v, outs[0])
        return None

    def closure(self, log):
        while True:
            F = self.forms()
            done = False
            for g in list(self.gorder):  # D1
                lit = self.gate_function(g, F)
                if lit is not None:
                    self.del_gate(g)
                    self.sub_value(g, lit)
                    log.append("D1")
                    done = True
                    break
            if done:
                continue
            for e, (terms, c) in list(self.eqs.items()):  # D3
                if len(terms) <= 2:
                    av = [v for v in terms if self.kind[v] == "aux"]
                    if not av:
                        raise Fail("equation without auxiliary term")
                    v = av[0]
                    rest = [u for u in terms if u != v]
                    del self.eqs[e]
                    self.auxes.discard(v)
                    self.sub_value(v, ("v", rest[0], c) if rest else ("c", c))
                    log.append("D3")
                    done = True
                    break
            if done:
                continue
            for g in list(self.gorder):  # D2
                if not self.children(g) and not (self.out[0] == "v" and self.out[1] == g):
                    self.del_gate(g)
                    log.append("D2")
                    done = True
                    break
            if done:
                continue
            for g in list(self.gorder):  # M
                t, A, B = self.gates[g]
                if t in XOR_TABLES and A[0] == B[0] == "v" and A[1] != B[1] and \
                        self.kind[A[1]] in ("param", "aux") and self.kind[B[1]] in ("param", "aux"):
                    self.del_gate(g)
                    self.kind[g] = "aux"
                    self.auxes.add(g)
                    self.eqs[self.neq] = [{g, A[1], B[1]}, tbl_val(t, 0, 0)]
                    self.neq += 1
                    log.append("M")
                    done = True
                    break
            if done:
                continue
            for x in sorted(self.params):  # D4
                if x in self.partner or (self.out[0] == "v" and self.out[1] == x):
                    continue  # an affine output is a Section 1 exclusion state
                ch = self.children(x)
                if len(ch) == 1 and next(iter(ch))[0] == "e":
                    e = next(iter(ch))[1]
                    Q = next((Q for Q, (mk, _) in sorted(F.items()) if mk >> x & 1), None)
                    if Q is None:
                        raise Fail("D4: no dependent auxiliary")
                    self.exchange(x, Q, F)
                    del self.eqs[e]
                    self.auxes.discard(x)
                    log.append("D4")
                    done = True
                    break
            if not done:
                return


def safe(st, p):
    return ("c", 1 - st.excl[p])


def step(st):
    """Perform one reducing row (with closure). Returns row name or None to stop."""
    log = []
    st.closure(log)
    if log:
        return "closure"
    if st.out[0] == "c" or st.kind[st.out[1]] in ("param", "aux"):
        return None  # Section 1 exclusion (affine/constant output)
    # paired rows
    for p in sorted(st.partner):
        if not st.children(p):
            r = st.partner[p]
            val = safe(st, r)
            st.unpair(p)
            st.fix_param(r, val)
            return "P0"
    for p in sorted(st.partner):
        if len(st.children(p)) >= 2:
            val = safe(st, p)
            st.unpair(p)
            st.fix_param(p, val)
            return "P2"
    for p in sorted(st.partner):
        r = st.partner[p]
        (hp,) = st.children(p)
        (hr,) = st.children(r)
        if hp == hr:
            pts = [(a, b) for a in (0, 1) for b in (0, 1) if (a, b) != (st.excl[p], st.excl[r])]
            if hp[0] == "g":
                t, A, B = st.gates[hp[1]]

                def hv(a, b, t=t, A=A, B=B):
                    env = {p: a, r: b}
                    return tbl_val(t, env[A[1]] if A[0] == "v" else A[1],
                                   env[B[1]] if B[0] == "v" else B[1])
            else:
                hv = lambda a, b: a ^ b
            pair = next((P1, P2) for i, P1 in enumerate(pts) for P2 in pts[i + 1:] if hv(*P1) == hv(*P2))
            (a1, b1), (a2, b2) = pair
            st.unpair(p)
            if a1 == a2:
                st.fix_param(p, ("c", a1))
            elif b1 == b2:
                st.fix_param(r, ("c", b1))
            else:
                st.fix_param(r, ("v", p, a1 ^ b1))
            return "PC"
    for p in sorted(st.partner):
        (hp,) = st.children(p)
        if hp[0] == "e":
            r = st.partner[p]
            val = safe(st, r)
            st.unpair(p)
            st.fix_param(r, val)
            st.closure([])
            return "PE"
    if not st.gorder:
        return None
    G = st.gorder[0]
    t, A, B = st.gates[G]
    if not is_and_type(t) or A[0] != "v" or B[0] != "v":
        raise Fail("first gate is not a nondegenerate AND")
    if st.out == ("v", G, 0) or st.out == ("v", G, 1):
        return None  # Section 1 exclusion
    F = st.forms()
    for k, P in ((1, A), (2, B)):
        if st.kind[P[1]] == "aux":
            other = (B if k == 1 else A)[1]
            mk = F[P[1]][0]
            cands = [x for x in sorted(st.params) if mk >> x & 1 and x != other]
            if not cands:
                raise Fail("AR: no exchange candidate")
            if any(x in st.partner for x in cands):
                raise Fail("AR: paired parameter in system")
            st.exchange(cands[0], P[1], F)
            return "AR"
    x, y = A[1], B[1]
    px, py = x in st.partner, y in st.partner
    if px and st.partner[x] == y:
        raise Fail("PC case reached main row")

    def kill_value(w):  # value of w forcing G constant
        for bw in (0, 1):
            if w == y and tbl_val(t, 0, bw) == tbl_val(t, 1, bw):
                return bw
            if w == x and tbl_val(t, bw, 0) == tbl_val(t, bw, 1):
                return bw
        raise Fail("no forcing value")

    if px != py:
        p, w = (x, y) if px else (y, x)
        r = st.partner[p]
        if st.children(r) and next(iter(st.children(r)))[0] != "g":
            raise Fail("AP1: partner child not a gate")
        vr = safe(st, r)
        st.unpair(p)
        st.fix_param(r, vr)
        st.fix_param(w, ("c", kill_value(w)))
        return "AP1"
    if px and py:
        r, h = st.partner[x], st.partner[y]
        vr, vh = safe(st, r), safe(st, h)
        st.unpair(x)
        st.unpair(y)
        st.fix_param(r, vr)
        st.fix_param(h, vh)
        st.fix_param(y, ("c", kill_value(y)))
        return "AP2"
    # GP
    corner = next((a, b) for a in (0, 1) for b in (0, 1)
                  if sum(tbl_val(t, a2, b2) == tbl_val(t, a, b) for a2 in (0, 1) for b2 in (0, 1)) == 1)
    const = tbl_val(t, 1 - corner[0], 1 - corner[1])
    st.partner[x], st.partner[y] = y, x
    st.excl[x], st.excl[y] = corner
    k = len(st.children(G))
    st.del_gate(G)
    st.sub_value(G, ("c", const))
    st.gp_k = k
    return "GP"


BOUND = {"P0": -3, "P2": 0, "PC": -4, "PE": 0, "AP1": -11, "AP2": -9, "AR": 0, "closure": 0}


AND_TABLES = [t for t in range(16) if is_and_type(t)]


def random_circuit(n, g, rng):
    gates = []
    for j in range(g):
        t = rng.randrange(16)
        a, b = rng.randrange(n + j), rng.randrange(n + j)
        gates.append((t, a, b))
    return gates, n + g - 1


def structured_circuit(n, extra, rng):
    """Random affine preprocessing, random shared AND-type layer, XOR of
    products of disjoint-ish affine forms (inner-product shape), plus random
    shared junk gates feeding the final XOR chain. Keeps the output
    nonlinear for many elimination steps."""
    gates = []

    def add(t, a, b):
        gates.append((t, a, b))
        return n + len(gates) - 1

    vals = list(range(n))
    for _ in range(rng.randrange(0, n)):  # affine mixing
        a, b = rng.sample(vals, 2)
        vals.append(add(rng.choice(XOR_TABLES), a, b))
    terms = []
    order = list(range(n))
    rng.shuffle(order)
    for i in range(0, n - 1, 2):
        u = order[i] if rng.random() < 0.6 else rng.choice(vals)
        v = order[i + 1] if rng.random() < 0.6 else rng.choice(vals)
        if u == v:
            continue
        terms.append(add(rng.choice(AND_TABLES), u, v))
    for _ in range(extra):
        pool = vals + terms
        a, b = rng.choice(pool), rng.choice(pool)
        terms.append(add(rng.choice(AND_TABLES + list(XOR_TABLES)), a, b))
    if not terms:
        terms.append(add(rng.choice(AND_TABLES), 0, 1))
    acc = terms[0]
    for tm in terms[1:]:
        acc = add(rng.choice(XOR_TABLES + tuple(AND_TABLES[:2])), acc, tm)
    return gates, acc


def ap2_circuit(n, extra, rng):
    """Two minimal ANDs whose GP pairs later meet at one AND (AP2 shape),
    with random tables, an optional shared child and random extra gates."""
    gates = []

    def add(t, a, b):
        gates.append((t, a, b))
        return n + len(gates) - 1

    x = list(range(n))
    rng.shuffle(x)
    A = add(rng.choice(AND_TABLES), x[0], x[1])
    B = add(rng.choice(AND_TABLES), x[2], x[3])
    C = add(rng.choice(AND_TABLES), x[0], x[2])
    D = add(rng.choice(AND_TABLES), x[1], x[4] if rng.random() < 0.7 else x[3])
    E = add(rng.choice(AND_TABLES), x[3], x[5])
    pool = [A, B, C, D, E] + x[4:]
    for _ in range(extra):
        pool.append(add(rng.choice(AND_TABLES + list(XOR_TABLES)), rng.choice(pool), rng.choice(pool)))
    acc = pool[0]
    for v in pool[1:]:
        if v >= n:
            acc = add(rng.choice(XOR_TABLES), acc, v)
    return gates, acc


def run(n, g, rng, stats, structured=False):
    if structured == "ap2":
        gates, out = ap2_circuit(n, g, rng)
    else:
        gates, out = structured_circuit(n, g, rng) if structured else random_circuit(n, g, rng)
    st = State(n, gates, out)
    W0, a0 = st.W(), st.a()
    st.check()
    while True:
        W1, a1 = st.W(), st.a()
        s1, I1, q1 = len(st.eqs) + len(st.gates), (st.W() - len(st.eqs) - len(st.gates) - 3 * (len(st.partner) // 2)) // 5, len(st.partner) // 2
        row = step(st)
        if row is None:
            break
        if row != "closure":
            st.closure([])
        st.check()
        dW, da = st.W() - W1, st.a() - a1
        stats[row] = stats.get(row, 0) + 1
        bound = -8 - st.gp_k if row == "GP" else BOUND[row]
        if da not in (0, -1):
            raise Fail(f"{row}: Delta a = {da}")
        if dW > 9 * da:
            raise Fail(f"{row}: Delta W = {dW} > 9 Delta a = {9 * da}")
        if dW > bound:
            key = f"table_overshoot_{row}_dW{dW}"
            stats[key] = stats.get(key, 0) + 1
            s2, q2 = len(st.eqs) + len(st.gates), len(st.partner) // 2
            I2 = (st.W() - s2 - 3 * q2) // 5
            ex = f"{row}: ds={s2 - s1} dI={I2 - I1} dq={q2 - q1} da={da} output_constant_after={st.out[0] == 'c'}"
            stats.setdefault("examples", set()).add(ex)
    if W0 - st.W() < 9 * (a0 - st.a()):
        raise Fail("global inequality violated")
    stats["dims_removed"] = stats.get("dims_removed", 0) + (a0 - st.a())


def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    rng = random.Random(seed)
    stats = {}
    for i in range(trials):
        structured = (False, True, "ap2")[i % 3]
        n = rng.randrange(6, 12) if structured else rng.randrange(3, 9)
        g = rng.randrange(0, 2 * n) if structured else rng.randrange(n, 4 * n + 8)
        try:
            run(n, g, rng, stats, structured)
        except Fail as exc:
            print(f"FAIL trial {i} n={n} g={g}: {exc}")
            sys.exit(1)
    print("all passed", dict(sorted(stats.items())))


if __name__ == "__main__":
    main()
