# Finite-field facts for the norm-one family (supporting file for proof.md)

These facts are used by `proof.md` Section 1.1 (the family, the rectangle
bound and the counting lemma) as indicated there.

**Fact LINPOLY.** A nonzero polynomial of degree D over an integral
domain has at most D roots (induction on D via the factor theorem: a
root r gives g = (T - r)q with deg q = D - 1). A finite integral domain
is a field (injectivity of w -> wx for w != 0 on a finite set).

**Fact FROB.** In characteristic 2, (u + v)^(2^c) = u^(2^c) + v^(2^c)
for all c >= 0 (the binomial coefficients (2^c choose i) are even for
0 < i < 2^c), so w -> w^(2^c) is additive; it is multiplicative, and on
a field it is injective with kernel {0}, hence an automorphism of any
finite field (the *c-th Frobenius power*). If W is a field with 2^d
elements, then w^(2^d) = w for all w in W (Lagrange in the group W* of
order 2^d - 1, plus w = 0), and the map phi: w -> w^2 has order exactly
d: phi^d = id by the previous sentence, and no smaller power c works,
else T^(2^c) - T of degree 2^c < 2^d would have 2^d roots in W,
contradicting LINPOLY. The powers of phi form a cyclic group of order
d; phi^c has order d/gcd(c, d).

**Fact KER (divisibility by T^(2^c) - T).** Let h in F_2[T] be monic
irreducible of degree d >= 1. Then h divides T^(2^c) - T in F_2[T] if
and only if d | c.

*Proof.* W_1 := F_2[T]/(h) is a field (h generates a maximal ideal of
the Euclidean domain F_2[T]) with 2^d elements, and the evaluation map
ev: F_2[T] -> W_1, g -> g(beta) with beta := the class of T, has kernel
exactly (h): the kernel is a proper ideal (1 |-> 1) containing (h),
whence equal to it by maximality. If d | c then beta^(2^c) = beta by
Fact FROB, so T^(2^c) - T lies in the kernel. If h | T^(2^c) - T then
beta^(2^c) = beta, and since every element of W_1 is g(beta) for some
g, Fact FROB gives g(beta)^(2^c) = g(beta^(2^c)) = g(beta): every
element of W_1 is a root of T^(2^c) - T, so 2^d <= 2^c, and applying
the same to the automorphism phi^c (which fixes beta and F_2, hence
fixes W_1 = F_2[beta] pointwise) gives phi^c = id on W_1, so
d = order(phi) divides c by Fact FROB. qed

**Fact IRRED (irreducibility, test and existence).** (a) A monic h in
F_2[T] of degree d is irreducible iff T^(2^d) = T mod h and
gcd(h, T^(2^(d/p)) - T) = 1 for every prime p | d. (b) F_2 has a monic
irreducible of every degree d >= 1.

*Proof.* (a) Suppose the test passes. Factor h = h_1...h_r into monic
irreducibles (unique factorization in the Euclidean domain F_2[T]),
with degrees d_i summing to d. By Fact KER applied modulo h_1 (where
T^(2^d) = T holds), d_1 | d; similarly every d_i | d. If r >= 2, some
d_i < d properly; pick that i and a prime p | d/d_i, so d_i | d/p.
Fact KER (forwards direction) gives h_i | T^(2^(d/p)) - T, so h_i |
gcd(h, T^(2^(d/p)) - T) != 1, contradicting the test. So r = 1 and h
is irreducible. Conversely let h be monic irreducible of degree d.
Fact KER gives T^(2^d) = T mod h. The gcd is a monic common divisor of
h, hence 1 or h; if it were h, then h | T^(2^(d/p)) - T and Fact KER
would give d | d/p, absurd. So the gcd is 1. (b) Let mu(d) be the
number of monic irreducibles of degree d. T^(2^d) - T is squarefree
(its derivative is 1 in characteristic 2) of degree 2^d, and by Fact
KER its monic irreducible factors are exactly those of degree c | d,
whence 2^d = sum_{c | d} c mu(c). Mobius inversion on the divisor
lattice gives mu(d) = (1/d) sum_{c | d} mobius(d/c) 2^c >= (1/d)(2^d -
sum_{1 <= c <= d-1} 2^c) = (2^d - 2^d + 2)/d = 2/d > 0: only the term
c = d contributes 2^d, every other c contributes at most 2^c, and
sum_{c=1}^{d-1} 2^c = 2^d - 2. In particular mu(m) >= 1 and
mu(t) >= 1: valid Q-blocks and R-blocks exist. qed

**Fact UNIT (roots of unity).** In any field, a finite subgroup of the
multiplicative group is cyclic; in a cyclic group of order M the
equation w^e = 1 has exactly gcd(e, M) solutions. Applied to the field
A_N (with |A_N*| = 2^(mt) - 1) and its subgroup
H_{Q,R} = {w in A_N : w^e = 1}: since e divides 2^(mt) - 1 (indeed
(2^m - 1) e = 2^(mt) - 1), |H_{Q,R}| = e exactly. Moreover e is odd
(e = 1 + even terms), and e >= 2^(mt-m) because
(2^m - 1) 2^(mt-m) = 2^(mt) - 2^(mt-m) <= 2^(mt) - 1.

*Proof of cyclicity.* Let G be a finite abelian group of order n and
exponent r (the lcm of the element orders). Every element of G is a
root of T^r - 1, so n <= r by LINPOLY (G embeds in the multiplicative
group of its field, a domain). The exponent is attained: for a in G
write ord(a) = prod p^u_p and a_p := a^(ord(a)/p^u_p), an element of
order p^u_p; given a, b, define c := the product over p of that one of
a_p, b_p having the larger order; the factors have pairwise coprime
orders, so their product c has order the product of the orders
(standard: if ord(u), ord(v) are coprime then ord(uv) =
ord(u)ord(v)), which is lcm(ord(a), ord(b)). Iterating over a list of
generators shows some element c* has order divisible by every element
order, hence by r; so r <= ord(c*) <= n. Combining, n = r = ord(c*)
and G = <c*>. For the solution count: a cyclic group of order M with
generator g: w = g^j satisfies w^e = 1 iff M | e j, i.e. (M/gcd(e,M))
| j: exactly gcd(e, M) values j in Z/M. qed

**Fact FIELD (composite fields).** Let d_1, d_2 >= 1 with
gcd(d_1, d_2) = 1, h_1 in F_2[X] and h_2 in F_2[Y] monic irreducible of
degrees d_1, d_2. Then the ideal (h_1(X), h_2(Y)) of F_2[X, Y] is
maximal: F_2[X,Y]/(h_1, h_2) is a field with 2^(d_1 d_2) elements.

*Proof.* E := F_2[X]/(h_1) is a field with 2^(d_1) elements, and
F_2[X,Y]/(h_1, h_2) = E[Y]/(h_2). We show h_2 is irreducible over E.
First, h_2 is squarefree over E: its derivative h_2' is a nonzero
polynomial of degree < d_2 over F_2 (if h_2' = 0 then h_2 = H(T)^2 for
H in F_2[T] by doubling exponents, making h_2 a nontrivial square
unless d_2 = 0, absurd), and gcd(h_2, h_2') = 1 already over F_2
(h_2 irreducible of larger degree), so a Bezout identity
u h_2 + v h_2' = 1 with u, v in F_2[Y] persists in E[Y]. Factor
h_2 = g_1 ... g_r over E[Y] into distinct monic irreducibles of degrees
D_1, ..., D_r (squarefree), with D_1 + ... + D_r = d_2. By the Chinese
remainder theorem (the g_i are pairwise coprime in the Euclidean domain
E[Y]), S := E[Y]/(h_2) is isomorphic to the product of the fields
F_i := E[Y]/(g_i), each of 2^(d_1 D_i) elements. In F_i let gamma be
the class of Y; h_2(gamma) = 0 (as g_i(gamma) = 0 is a factor), so by
Fact KER (h_2 | T^(2^(d_2)) - T) we get gamma^(2^(d_2)) = gamma. The
minimal polynomial of gamma over F_2 is h_2 (gamma is a root of the
irreducible h_2), so the ring F_2[gamma] has 2^(d_2) elements and is a
field (LINPOLY: finite domain). Put g' := gcd(d_2, d_1 D_i). The
automorphisms phi^(d_2) and phi^(d_1 D_i) of F_i (Fact FROB) both fix
gamma; since the powers of phi form a cyclic group of order d_1 D_i,
also phi^(g') fixes gamma, hence fixes F_2[gamma] pointwise
(Fact FROB additivity). The fixed set of phi^(g') is contained in the
root set of T^(2^(g')) - T in F_i, of at most 2^(g') elements
(LINPOLY). So 2^(d_2) <= 2^(g') <= 2^(d_2), forcing g' = d_2, i.e.
d_2 | d_1 D_i, and since gcd(d_1, d_2) = 1, d_2 | D_i. As D_i <= d_2,
D_i = d_2 for every i, so r = 1: h_2 is irreducible over E, and
S = E[Y]/(h_2) is a field of 2^(d_1 d_2) elements. qed

