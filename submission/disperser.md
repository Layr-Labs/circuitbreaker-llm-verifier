# Part D. The affine disperser: complete proof (Ben-Sasson–Kopparty, binary case, prime block degree)

This file proves the disperser property used in `proof.md`. It re-derives, for p = 2 and prime block degree k, Sections 3, 6 and 7 of E. Ben-Sasson and S. Kopparty, "Affine dispersers from subspace polynomials" (author preprint 2010-03-12, https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf, SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf; journal version SIAM J. Comput. 41(4):880–914, 2012). Every step is written out; nothing is used by citation except textbook facts named in the text (finite field arithmetic, uniqueness of binary representation, Gaussian elimination, a nonzero polynomial of degree less than |F| has a nonzero value on a field F).

Throughout: k is a prime, k >= 3, and K = F_{2^k} is a field with 2^k elements (any concrete model, e.g. F_2[Z]/(q) with q irreducible of degree k). Z_k denotes the residues mod k; [a, b] denotes the integers from a to b, and "[a, b] mod k" the cyclic interval of residues a, a+1, ..., b taken mod k (empty if b < a). For z in K, z^{2^k} = z.

**Trace.** Tr : K -> K, Tr(z) = sum_{s=0}^{k-1} z^{2^s}. It is F_2-linear (Frobenius z -> z^2 is additive), and Tr(z)^2 = sum_s z^{2^{s+1}} = Tr(z) since z^{2^k} = z, so Tr(z) in {0,1} = F_2 for every z in K. Tr is not the zero function: the polynomial sum_s X^{2^s} has degree 2^{k-1} < |K| and is not the zero polynomial, so it has a nonzero value on K.

## D0. Statement

**Theorem D.** Let 1 <= t <= r be integers and let k be a prime with k >= r(r+1)/t. For i in [1, r] put e_i := r + 1 - i and alpha_i := sum_{l=0}^{e_i} 2^l = 2^{r+2-i} - 1. Define F : K^r -> F_2 by

    F(x_1, ..., x_r) = Tr( sum_{I subset [1,r], |I| = t} prod_{i in I} x_i^{alpha_i} ).

View K^r as an F_2-vector space of dimension rk. If A subset K^r is an F_2-affine subspace with dim A > kr/t + kt + r(r+1), then F is not constant on A.

The proof occupies D1–D6. D7 applies it to the concrete family of `proof.md`.

## D1. Polynomials over K, reduction and binary weight

**D1.1 (nonzero polynomials).** Let h in K[Z_1, ..., Z_r] have degree at most 2^k - 1 in each variable. If h is not a constant polynomial, then the function h : K^r -> K is not constant. Proof: h - h(0,...,0) is a nonzero polynomial with all individual degrees < |K|. A nonzero polynomial with individual degrees < |K| does not vanish identically on K^r (induction on r: write it as sum_j c_j(Z_1..Z_{r-1}) Z_r^j; some c_j is nonzero, hence by induction nonzero at some point of K^{r-1}, where the resulting univariate polynomial of degree < |K| is nonzero and so has a nonzero value). Hence h takes a value different from h(0,...,0).

**Reduction.** For an integer e >= 1 let rd(e) be the unique element of [1, 2^k - 1] congruent to e modulo 2^k - 1, and rd(0) = 0. For a polynomial h, red(h) is obtained by replacing every monomial exponent e of every variable by rd(e). Since z^{2^k} = z for z in K, z^e = z^{rd(e)} for all z in K and e >= 0, so h and red(h) are the same function on K^r, and red(h) has individual degrees <= 2^k - 1. red(h_1 h_2) = red(red(h_1) red(h_2)); red is K-linear; and for polynomials in disjoint sets of variables the coefficient of a product monomial in red(h_1 h_2) is the product of the coefficients of the factors' monomials in red(h_1), red(h_2).

**Binary weight.** wt(e) is the number of ones in the binary expansion of the integer e >= 0.

**D1.2.** (i) wt(e + e') <= wt(e) + wt(e'). (ii) wt(rd(e)) <= wt(e) for all e >= 0.
Proof. (i) Adding in binary, each carry replaces two ones by one, so wt(e+e') = wt(e) + wt(e') - (number of carries). (ii) For e < 2^k, rd(e) = e. For e >= 2^k write e = 2^k u + v with 0 <= v < 2^k and u >= 1; then u + v = e - (2^k - 1)u is congruent to e mod 2^k - 1, u + v < e, u + v >= 1, and wt(u+v) <= wt(u) + wt(v) = wt(e) by (i). Iterating, we reach a number e* in [1, 2^k - 1] congruent to e with wt(e*) <= wt(e); e* = rd(e).

**D1.3.** For a monomial and a variable Z_i, its Z_i-weight is wt of its Z_i-exponent. If h_1, ..., h_m are polynomials and every monomial of h_j has Z_i-weight at most w_j, then every monomial of red(h_1 ... h_m) has Z_i-weight at most w_1 + ... + w_m. (A monomial of the product has Z_i-exponent equal to a sum of exponents from the factors; apply D1.2(i) then D1.2(ii).)

## D2. Linearized and subspace polynomials

A **linearized polynomial** is P(X) = sum_{j=0}^{k-1} a_j X^{2^j} with a_j in K. The function P : K -> K is F_2-linear (each z -> z^{2^j} is additive and fixes F_2). An **affine linearized** polynomial is a linearized polynomial plus a constant. A product of linearized polynomials is not linearized in general, but P(X)^2 = sum a_j^2 X^{2^{j+1}} and X^{2^k} is identified with X on K.

**D2.1 (linear maps are linearized polynomials).** Every F_2-linear map T : K -> K equals, as a function on K, a unique linearized polynomial of degree < 2^k. Proof: linearized polynomials of degree < 2^k form a K-vector space of dimension k, hence an F_2-space of dimension k^2; Hom_{F_2}(K, K) also has F_2-dimension k^2; the evaluation map from polynomials to functions is F_2-linear and injective (a nonzero polynomial of degree < |K| is not the zero function, D1.1), hence bijective.

**D2.2 (kernel polynomial).** Let L subset K be an F_2-linear subspace of dimension d, 0 <= d <= k. Then P_L(X) := prod_{a in L} (X - a) is a monic polynomial of degree 2^d of the form sum_{j=0}^{d} c_j X^{2^j} (c_d = 1); the function P_L : K -> K is F_2-linear with kernel exactly L, so its image L' := P_L(K) is an F_2-linear subspace of dimension k - d.
Proof. Induction on d. d = 0: P = X. If L = L_0 + F_2 b with b not in L_0 and dim L_0 = d - 1, then P_L(X) = P_{L_0}(X) P_{L_0}(X + b) = P_{L_0}(X) (P_{L_0}(X) + P_{L_0}(b)) = P_{L_0}(X)^2 + P_{L_0}(b) P_{L_0}(X), using F_2-linearity of P_{L_0} (induction) for the middle equality. This is of the stated form, monic of degree 2^d. The roots of P_L are exactly the elements of L (each with multiplicity one), so the kernel of the linear map is L and the image has dimension k - d (rank–nullity).

**D2.3 (image polynomial).** Let L subset K be linear of dimension d with 1 <= d <= k, and let L' := P_L(K) (dimension k - d). Put Q_L := P_{L'}, a monic linearized polynomial of degree 2^{k-d} <= 2^{k-1} whose roots are exactly the 2^{k-d} elements of L'. Then Q_L(K) = L.
Proof. Q_L(P_L(X)) is monic of degree 2^k and vanishes at every z in K (P_L(z) in L' is a root of Q_L), so Q_L(P_L(X)) = X^{2^k} - X as polynomials. Hence P_L(Q_L(P_L(X))) = P_L(X^{2^k} - X) = P_L(X)^{2^k} - P_L(X) (P_L is a sum of monomials c_j X^{2^j} with c_j^{2^k} = c_j, so P_L(X^{2^k}) = P_L(X)^{2^k}, and P_L(X - Y) = P_L(X) - P_L(Y) as polynomials). Thus g(Y) := P_L(Q_L(Y)) - (Y^{2^k} - Y) satisfies g(P_L(X)) = 0 as a polynomial; since deg g(P_L(X)) = 2^d deg g, g = 0. So P_L(Q_L(y)) = 0 for every y in K, i.e. Q_L(K) subset ker P_L = L. Q_L is linear with kernel L' of dimension k - d, so dim Q_L(K) = d = dim L and Q_L(K) = L.

Write Q_L(X) = sum_{j=0}^{k-1} a_j X^{2^j}. Then a_{k-d} = 1 and a_j = 0 for j in [k-d+1, k-1] (degree 2^{k-d}, monic). For 0 <= s <= k-1 define the **Frobenius shift** Q^{(s)}(X) := sum_{j=0}^{k-1} a_{(j-s) mod k}^{2^s} X^{2^j}. For every z in K, Q^{(s)}(z) = Q(z)^{2^s}: indeed Q(z)^{2^s} = sum_j a_j^{2^s} z^{2^{j+s}} and z^{2^{j+s}} = z^{2^{(j+s) mod k}}.

**D2.4 (structure of the coefficients; BK Lemma 3.10 for prime k).** Let L be linear of dimension d, 1 <= d <= k-1, with Q = Q_L and coefficients a_j as above; write a'_j = 1 if a_j != 0 and a'_j = 0 otherwise.
(1) Every cyclic interval J = [rho, rho + d - 1] mod k of d consecutive residues contains some j with a'_j = 1. In particular a'_0 = 1.
(2) If d >= 2, every cyclic interval J = [rho, rho + d - 2] mod k of d - 1 consecutive residues with rho != (k - d + 1) mod k contains some j with a'_j = 1.
Proof. (1) Suppose a_j = 0 for all j in J. Let s := (k - d - rho) mod k, so that J + s = [k - d, k - 1] mod k. Then Q^{(s)} has zero coefficient at every index in [k-d, k-1], so deg Q^{(s)} <= 2^{k-d-1}. Q^{(s)} is not the zero polynomial: since a_{k-d} = 1 != 0, the index k - d is not in J, so (k - d + s) mod k is not in J + s = [k-d, k-1], and the coefficient of Q^{(s)} at that index is a_{k-d}^{2^s} = 1. But Q^{(s)}(z) = Q(z)^{2^s} = 0 for all z in L', which has 2^{k-d} elements: a nonzero polynomial of degree < 2^{k-d} with 2^{k-d} roots, contradiction. For a'_0 = 1: if a_0 = 0 then [k-d+1, k] mod k = {k-d+1, ..., k-1, 0} would be d consecutive zero indices.
(2) Suppose a_j = 0 for all j in J. Let s := (k - d + 1 - rho) mod k; since rho != k - d + 1 mod k, 1 <= s <= k - 1. Now J + s = [k-d+1, k-1] mod k, so Q^{(s)} has zero coefficients at indices k-d+1, ..., k-1 and deg Q^{(s)} <= 2^{k-d}. Let c be the coefficient of X^{2^{k-d}} in Q^{(s)}. Then Q^{(s)} - c Q has degree <= 2^{k-d-1} and vanishes on the 2^{k-d} elements of L' (Q^{(s)}(z) = Q(z)^{2^s} = 0 = Q(z) there), so Q^{(s)} = c Q as polynomials. Hence for every y in L = Q(K), say y = Q(z): y^{2^s} = Q(z)^{2^s} = Q^{(s)}(z) = c Q(z) = c y. Fix y_0 in L, y_0 != 0 (d >= 1). Then c = y_0^{2^s - 1} != 0, and for any y in L \ {0}, w := y / y_0 satisfies w^{2^s - 1} = 1, i.e. w^{2^s} = w. Also w^{2^k} = w. Let M := { m >= 1 : w^{2^m} = w }. M contains s and k, is closed under addition (w^{2^{m+m'}} = (w^{2^m})^{2^{m'}} = w^{2^{m'}} = w), and if m > m' are in M then w = w^{2^m} = (w^{2^{m'}})^{2^{m-m'}} = w^{2^{m-m'}}, so m - m' in M. By the Euclidean algorithm gcd(s, k) in M; k is prime and 1 <= s <= k-1, so 1 in M: w^2 = w, w in {0, 1}, hence y = y_0. So L \ {0} = {y_0}, contradicting d >= 2.

(The hypothesis "k prime" is used exactly here, in (2). For composite k the conclusion would only be L subset y_0 F_{2^{gcd(s,k)}}, which is BK's hypothesis "L not contained in an affine shift of a proper subfield"; for prime k that hypothesis is automatic once d >= 2.)

## D3. The associated matrix and the coefficient formula

For a linearized Q(X) = sum_{j=0}^{k-1} a_j X^{2^j} define the k x k matrix M_Q over K, rows and columns indexed by Z_k, by M_Q[u, v] := a_{(v-u) mod k}^{2^u}. Then Q^{(u)}(X) = sum_v M_Q[u, v] X^{2^v}. For index sets A, B subset Z_k with |A| = |B|, M_Q[A, B] is the square submatrix with rows A and columns B, and Perm(M_Q[A, B]) := sum over bijections sigma : A -> B of prod_{u in A} M_Q[u, sigma(u)] (independent of any ordering). For s in Z_k, A + s := {(u + s) mod k : u in A}, B - s likewise. The 0/1 indicator matrix M'_Q has M'_Q[u, v] = a'_{(v-u) mod k}.

**D3.1 (Frobenius shift of minors).** M_Q[u + s, v] = (M_Q[u, v - s])^{2^s}, hence Perm(M_Q[A + s, B]) = Perm(M_Q[A, B - s])^{2^s}. Proof: M_Q[u+s, v] = a_{(v-u-s)}^{2^{u+s}} = (a_{((v-s)-u)}^{2^u})^{2^s}; z -> z^{2^s} is a ring homomorphism of K, so it commutes with sums and products, hence with permanents.

**D3.2 (coefficient formula; BK Lemma 3.13).** Let mu in K \ {0}, let Q_1, ..., Q_t be linearized polynomials with matrices M_1, ..., M_t, and let A_i, B_i subset Z_k with |A_i| = |B_i| >= 1 (i = 1..t). Put alpha_i := sum_{j in A_i} 2^j and beta_i := sum_{b in B_i} 2^b, and

    R(X_1, ..., X_t) := red( Tr( mu prod_{i=1}^t Q_i(X_i)^{alpha_i} ) ),   Tr(h) := sum_{s=0}^{k-1} h^{2^s} (as polynomials).

Then the coefficient c of the monomial prod_i X_i^{beta_i} in R is
    c = sum_{s=0}^{k-1} mu^{2^s} prod_{i=1}^t Perm(M_i[A_i + s, B_i]) = sum_{s=0}^{k-1} mu^{2^s} prod_{i=1}^t Perm(M_i[A_i, B_i - s])^{2^s}.
Proof. The s-th summand of Tr is mu^{2^s} prod_i Q_i(X_i)^{alpha_i 2^s}, and Q_i^{alpha_i 2^s} = prod_{j in A_i} Q_i^{2^{j+s}}; as a function on K, Q_i(X)^{2^{j+s}} = Q_i^{((j+s) mod k)}(X) = sum_v M_i[j+s, v] X^{2^v}, and these two polynomials have the same reduction (both have individual degree < 2^k after reduction and agree as functions: by D1.1 applied to their difference they are equal after reduction). Expand prod_{j in A_i} sum_v M_i[j+s, v] X_i^{2^v}: choosing v_j for each j in A_i gives the monomial X_i^{sum_j 2^{v_j}} with coefficient prod_j M_i[j+s, v_j]. If the v_j are pairwise distinct, sum_j 2^{v_j} < 2^k is the integer with digit set {v_j}; it equals beta_i iff {v_j : j in A_i} = B_i (uniqueness of binary representation). If two v_j coincide, the exponent has weight < |A_i| = wt(beta_i) by D1.2(i), and after reduction still weight < wt(beta_i) by D1.2(ii), so it is not beta_i. Hence the coefficient of X_i^{beta_i} in red(Q_i(X_i)^{alpha_i 2^s}) is sum over bijections v : A_i -> B_i of prod_{j} M_i[j+s, v(j)] = Perm(M_i[A_i + s, B_i]). The variables X_i are distinct, so the coefficient of prod_i X_i^{beta_i} in the reduced s-th summand is mu^{2^s} times the product over i. Summing over s gives the first formula; D3.1 gives the second.

**D3.3 (nonvanishing criterion; BK Remark 3.14).** In the setting of D3.2, suppose
(a) for each i there are orderings of A_i and of B_i such that the matrix M'_i[A_i, B_i], written with rows in the order of A_i and columns in the order of B_i, is lower triangular with every diagonal entry equal to 1; and
(b) for each s in [1, k-1] there is an i such that some column of M'_i[A_i, B_i - s] is all zero.
Then c != 0.
Proof. By (b), for s >= 1 some M_i[A_i, B_i - s] has an all-zero column, so its permanent is 0 and the s-th term vanishes. For s = 0: with the orderings of (a), any bijection other than the order-preserving one uses an entry above the diagonal, which is 0; so Perm(M_i[A_i, B_i]) is the product of the diagonal entries, all nonzero. Hence c = mu prod_i Perm(M_i[A_i, B_i]) != 0.

## D4. Independent sources (BK Theorem 6.1, binary, prime k)

**Theorem D4.1.** Let Y_1, ..., Y_t subset K be F_2-linear subspaces of dimensions d_1, ..., d_t with d_1 = max_i d_i and d_i >= 3 for all i. Let Q_i := Q_{Y_i} with coefficients a_{i,j} (D2.3) and indicator bits a'_{i,j}. Let e_1, ..., e_t be integers with 1 <= e_i <= d_i - 2, alpha_i := sum_{l=0}^{e_i} 2^l, and mu in K \ {0}. If
    sum_{i=1}^t (d_i - e_i - 1) >= k - d_1 + 1,
then there are sets B_i subset Z_k with |B_i| = e_i + 1 such that, with beta_i := sum_{b in B_i} 2^b (so wt(beta_i) = e_i + 1), the monomial prod_i X_i^{beta_i} has a nonzero coefficient in R := red(Tr(mu prod_i Q_i(X_i)^{alpha_i})).

Proof. Put A_i := {0, 1, ..., e_i}, so alpha_i = sum_{j in A_i} 2^j. We construct the B_i and verify (a) and (b) of D3.3. Recall (D2.3) a'_{i, k - d_i} = 1 and a'_{i,j} = 0 for j in [k-d_i+1, k-1] (the "zero zone" of Q_i, a cyclic interval of length d_i - 1), and a'_{i,0} = 1 (D2.4(1) when d_i < k; when d_i = k, Q_i = X by D2.3, so a_{i,0} = 1, and D2.4(2) is never invoked because every chosen l_i is then zero).

**Positions.** Let l_1 := 0 and l_{i+1} := l_i + (d_i - e_i - 1) for 1 <= i <= t (ordinary integers, each increment >= 1). By hypothesis l_{t+1} >= k - d_1 + 1; let t_0 be the least index with l_{t_0 + 1} >= k - d_1 + 1. For i <= t_0, minimality gives l_i <= l_{t_0} <= k - d_1 <= k - d_i, and l_{t_0+1} = l_{t_0} + d_{t_0} - e_{t_0} - 1 <= (k - d_1) + d_{t_0} - 2 <= k - 2. For i > t_0 we redefine l_i := 0 (these indices are only needed for (a)).

**Definition of B_i.** Fix i and write d = d_i, e = e_i, l = l_i, a'_j = a'_{i,j}.
Case (alpha): a'_l = 1 (this includes l = 0). B_i := {l} union {k-d+1, k-d+2, ..., k-d+e}.
Case (beta): a'_l = 0. Then l >= 1 and (since a'_{k-d} = 1 and l <= k - d) l <= k - d - 1. Let j := min{ j' > l : a'_{j'} = 1 } (exists, j <= k - d), jh := max{ j' < l : a'_{j'} = 1 } (exists, jh >= 0), q := j - jh - 1 >= 1 (the indices jh+1, ..., j-1 all have a' = 0, and l is among them), s_i := min(q, e). B_i := {j} union {jh+1, ..., jh+s_i} union {k-d+s_i+1, ..., k-d+e}.
In both cases |B_i| = e + 1 and B_i subset [0, k-2] subset Z_k: in case (alpha), l is not in [k-d+1, k-d+e] because a'_l = 1 while that range lies in the zero zone (e <= d - 2 gives k-d+e <= k-2); in case (beta), the three parts are pairwise disjoint since jh + s_i <= jh + q = j - 1 < j <= k - d < k-d+s_i+1, and all lie in [0, k-2].

**Claim D4.2 (condition (a)).** Order A_i as 0, 1, ..., e and order B_i as listed above (case (alpha): l, k-d+1, ..., k-d+e; case (beta): j, jh+1, ..., jh+s_i, k-d+s_i+1, ..., k-d+e). Then M'_i[A_i, B_i] is lower triangular with diagonal entries 1. The entry in row u and column b is a'_{(b-u) mod k}.
Case (alpha). Column 0 is b = l: diagonal entry (u = 0) is a'_l = 1. Column c in [1, e] is b = k-d+c: for u < c the index is k-d+(c-u) in [k-d+1, k-d+e], in the zero zone, so entries above the diagonal are 0; the diagonal entry (u = c) is a'_{k-d} = 1.
Case (beta). Column 0 is b = j: diagonal entry a'_j = 1. Column c in [1, s_i] is b = jh + c: for u < c the index is jh + (c - u) in [jh+1, jh+s_i] subset [jh+1, j-1], where a' = 0; the diagonal entry (u = c) is a'_{jh} = 1. Column c in [s_i+1, e] is b = k-d+c: as in case (alpha), entries with u < c have index in [k-d+1, k-d+e] (zero zone) and the diagonal entry is a'_{k-d} = 1. (No index arithmetic wraps mod k: all indices involved lie in [0, k-1].)

**Claim D4.3 (condition (b)).** For every s in [1, k-1] there is i <= t_0 such that M'_i[A_i, B_i - s] has an all-zero column.
First, for b in B_i the column of M'_i[A_i, B_i - s] indexed by (b - s) mod k has entries a'_{i,(b-s-u) mod k}, u = 0, ..., e_i, i.e. the indices of the cyclic interval [b-s-e_i, b-s] mod k. All these lie in the zero zone [k-d_i+1, k-1] iff (b - s) mod k lies in [k-d_i+1+e_i, k-1], iff (s - b) mod k lies in [1, d_i - e_i - 1]. So:
(*) if s is in the cyclic interval [b+1, b+d_i-e_i-1] mod k for some b in B_i, then M'_i[A_i, B_i - s] has an all-zero column.
Let C_i := union_{b in B_i} [b+1, b+d_i-e_i-1] mod k. It suffices to show [1, k-1] subset union_{i <= t_0} C_i.
Block 1: l_1 = 0 and a'_{1,0} = 1, so case (alpha): B_1 = {0} union [k-d_1+1, k-d_1+e_1]. From b = 0: C_1 contains [1, d_1-e_1-1] = [l_1+1, l_2]. From b = k-d_1+c, c = 1..e_1: the intervals [k-d_1+c+1, k+c-e_1-1] mod k; consecutive ones overlap or abut (since d_1 >= e_1 + 2), and their union is [k-d_1+2, k-1] (the last interval, c = e_1, ends at k-1). So C_1 contains [1, l_2] union [k-d_1+2, k-1].
Blocks 1 < i <= t_0: we show C_i contains [l_i+1, l_{i+1}]. Note 1 <= l_i <= k - d_i.
 - Case (alpha): l_i in B_i and the interval for b = l_i is [l_i+1, l_i+d_i-e_i-1] = [l_i+1, l_{i+1}].
 - Case (beta), q <= e_i (so s_i = q): then {jh+1, ..., jh+q} union {j} = [jh+1, j] contains l_i, so l_i in B_i and we conclude as in case (alpha).
 - Case (beta), q > e_i (so s_i = e_i and B_i = {j} union [jh+1, jh+e_i]): the indices jh+1, ..., jh+q = j-1 all have a' = 0. Apply D2.4(2) to Q_i with rho = jh+1: rho <= l_i <= k - d_i < k - d_i + 1, so rho != k - d_i + 1 and every run of d_i - 1 consecutive zero indices starting at rho is impossible; hence q <= d_i - 2, i.e. j <= jh + d_i - 1. The intervals for b = jh+1, ..., jh+e_i are [b+1, b+d_i-e_i-1]; consecutive b give abutting intervals, so their union is [jh+2, jh+e_i+d_i-e_i-1] = [jh+2, jh+d_i-1], which contains j. Together with the interval [j+1, j+d_i-e_i-1] for b = j, C_i contains [jh+2, j+d_i-e_i-1]. Since jh < l_i < j, this contains [l_i+1, l_i+d_i-e_i-1] = [l_i+1, l_{i+1}].
Hence union_{i <= t_0} C_i contains [1, l_{t_0+1}] union [k-d_1+2, k-1], and l_{t_0+1} >= k - d_1 + 1, so it contains [1, k-1]. (All intervals here are ordinary integer intervals inside [1, k-1] since l_{t_0+1} <= k-2.)

By D3.3 the coefficient of prod_i X_i^{beta_i} in R is nonzero. QED D4.1.

## D5. Echelon decomposition of an affine subspace (BK Lemma 7.1)

**D5.1.** Fix an F_2-basis of K, identifying K^r with F_2^{rk} (coordinates grouped into r blocks of k). Let A subset K^r be an F_2-affine subspace of dimension d and gamma in A. There are F_2-linear subspaces Y_1, ..., Y_r subset K with sum_i dim Y_i = d, and F_2-linear maps sigma_{ij} : Y_j -> K for 1 <= j < i <= r, such that
    A = { (gamma_1 + y_1 + sum_{j<1} ..., ..., gamma_i + y_i + sum_{j<i} sigma_{ij}(y_j), ...) : y_i in Y_i for all i }.
Proof. A = gamma + V with V linear of dimension d. Take a d x rk matrix over F_2 whose rows form a basis of V and bring it to reduced row echelon form by Gaussian elimination (row operations keep the row space V). Each row v^{(s)} has a pivot column p_s (its first nonzero entry), the pivots are distinct, and every row has 0 in the pivot columns of the other rows. Group the rows by the block containing their pivot: G_i := { s : p_s in block i }. Define Y_i := the span of the block-i parts of the rows in G_i; these parts are linearly independent (they have distinct pivot positions inside block i with the echelon pattern), so dim Y_i = |G_i| and sum_i dim Y_i = d. Every v in V is uniquely v = sum_s c_s v^{(s)}, and c_s is the coordinate of v at the pivot p_s. The block-i part of v is y_i + sum_{j < i} sum_{s in G_j} c_s (v^{(s)})_i where y_i := sum_{s in G_i} c_s (v^{(s)})_i in Y_i, because a row with pivot in block j has zeros in all blocks before j (its first nonzero entry is in block j). The coefficients (c_s)_{s in G_j} are determined linearly by y_j (they are the pivot coordinates of y_j), so sigma_{ij}(y_j) := sum_{s in G_j} c_s (v^{(s)})_i is a linear map Y_j -> K. Conversely every choice of (y_i in Y_i)_i arises from a unique (c_s). QED.

## D6. Proof of Theorem D

Let A subset K^r be affine with d := dim A > kr/t + kt + r(r+1); let gamma, Y_i, sigma_{ij} be as in D5.1 and d_i := dim Y_i, so sum_i d_i = d and 0 <= d_i <= k.

**Polynomial parametrization.** For each i with d_i >= 1 let Q_i := Q_{Y_i} (D2.3), a linearized polynomial with Q_i(K) = Y_i; if d_i = 0 put Q_i := 0 (then Q_i(K) = {0} = Y_i). For j < i, sigma_{ij} o Q_j is an F_2-linear map K -> K, so by D2.1 there is a linearized polynomial Q_{ij} with Q_{ij}(z) = sigma_{ij}(Q_j(z)) for all z in K. Put
    R_i(Z_1, ..., Z_r) := Q_i(Z_i) + sum_{j<i} Q_{ij}(Z_j) + gamma_i,
an affine linearized polynomial in Z_1, ..., Z_i (every monomial has total weight <= 1 and involves only one variable). The map H := (R_1, ..., R_r) : K^r -> K^r has image exactly A: as z_i ranges over K independently, y_i := Q_i(z_i) ranges over Y_i independently, and H(z) = (gamma_i + y_i + sum_{j<i} sigma_{ij}(y_j))_i. Therefore F(A) = (F o H)(K^r), and F o H is the function of the polynomial
    h := red( sum_{I subset [1,r], |I| = t} Tr( prod_{i in I} R_i^{alpha_i} ) ),
(with Tr applied to polynomials as in D3.2, mu = 1). By D1.1 it suffices to show that h has a nonconstant monomial with nonzero coefficient.

**Choice of the heavy blocks.** Let S := { i : d_i >= r + 2 }. Since d_i <= k for all i and d_i <= r + 1 for i not in S: |S| k + r(r+1) >= sum_i d_i = d > kt + r(r+1), so |S| > t, in particular |S| >= t. Also sum_{i in S} d_i >= sum_{i} (d_i - (r+1)) = d - r(r+1) > kr/t + kt. Let I_0 subset S consist of t indices with the largest d_i (ties arbitrary). Then sum_{i in I_0} d_i >= (t/|S|) sum_{i in S} d_i >= (t/r) sum_{i in S} d_i > (t/r)(kr/t + kt) = k + kt^2/r, so
    sum_{i in I_0} (d_i - (r+1)) > k + kt^2/r - (r+1)t >= k,
using kt >= r(r+1) (hypothesis k >= r(r+1)/t). For i in I_0: d_i >= r + 2 >= e_i + 2 (as e_i = r+1-i <= r), e_i >= 1 (i <= r), and d_i - e_i - 1 >= d_i - (r+1). Hence, with d_max := max_{i in I_0} d_i >= 1,
    sum_{i in I_0} (d_i - e_i - 1) > k >= k - d_max + 1.
So D4.1 applies to the subspaces (Y_i)_{i in I_0} (listed with a maximal-dimension one first), mu = 1, and the exponents e_i: there are B_i (i in I_0) with |B_i| = e_i + 1 such that the monomial
    M := prod_{i in I_0} Z_i^{beta_i},   beta_i = sum_{b in B_i} 2^b,  wt(beta_i) = e_i + 1,
has a nonzero coefficient in g' := red( Tr( prod_{i in I_0} Q_i(Z_i)^{alpha_i} ) ). M is nonconstant (each beta_i >= 1).

**Lemma D6.1 (BK Lemma 7.2).** The coefficient of M in g := red( Tr( prod_{i in I_0} R_i^{alpha_i} ) ) equals its coefficient in g'.
Proof. Write Ph_i := sum_{j<i} Q_{ij}(Z_j) + gamma_i, which involves only Z_j with j < i. For each s, the s-th trace summand of prod_{i in I_0} R_i^{alpha_i} is prod_{i in I_0} prod_{l=0}^{e_i} R_i^{2^{s+l}} and R_i^{2^m} = Q_i(Z_i)^{2^m} + Ph_i^{2^m}. Expanding prod_{l=0}^{e_i} (Q_i(Z_i)^{2^{s+l}} + Ph_i^{2^{s+l}}) = Q_i(Z_i)^{alpha_i 2^s} + G_i, every monomial of G_i (before reduction) has Z_i-exponent a sum of at most e_i powers of two, hence Z_i-weight <= e_i (D1.2(i)), and G_i involves only Z_j with j <= i. Now prod_{i in I_0} (Q_i^{alpha_i 2^s} + G_i) = prod_i Q_i(Z_i)^{alpha_i 2^s} + sum_{nonempty T subset I_0} prod_{i in T} G_i prod_{i in I_0 \ T} Q_i(Z_i)^{alpha_i 2^s}. For a nonempty T let i* := max T. In the corresponding product, the factor G_{i*} contributes Z_{i*}-weight <= e_{i*}; the factors G_i with i in T, i < i*, involve only Z_j with j <= i < i*; the factors Q_i(Z_i)^{...} with i != i* involve only Z_i. By D1.3 every monomial of the reduced product has Z_{i*}-weight <= e_{i*} < e_{i*} + 1 = wt(beta_{i*}), so M does not occur in it. Hence, after reduction, the coefficient of M in the s-th summand of g equals that in the s-th summand of g'; sum over s.

**Lemma D6.2 (BK Claim 7.4 and its use).** Let J subset [1, r] with |J| = t and J != I_0. Then the coefficient of M in red( Tr( prod_{j in J} R_j^{alpha_j} ) ) is 0.
Proof. The s-th trace summand is prod_{j in J} prod_{l=0}^{e_j} R_j^{2^{s+l}}; each factor R_j^{2^m} is a K-linear combination of the constant 1 and monomials Z_i^{2^b} with i <= j (R_j is affine linearized in Z_1..Z_j, and raising to 2^m maps a Z_i^{2^c} to Z_i^{2^{c+m}}). A monomial of the product (before reduction) is obtained by choosing, for each "slot" (j, l) with j in J and 0 <= l <= e_j = r+1-j, either the constant or a monomial in a variable Z_{lambda_j(l)} with lambda_j(l) <= j. Let n_i be the number of slots choosing Z_i. The Z_i-exponent of the chosen monomial is a sum of n_i powers of two, so its Z_i-weight is <= n_i, and reduction does not increase it (D1.2). For the result to be M we need n_i >= wt(beta_i) = e_i + 1 = r + 2 - i for every i in I_0. Claim D6.3 shows this forces J = I_0. So no monomial of any summand reduces to M.

**Claim D6.3.** Let I_0, J subset [1, r], |I_0| = |J| = t, and suppose partial functions lambda_j : {0, ..., r+1-j} -> {1, ..., j} (j in J) satisfy n_i := sum_{j in J} |lambda_j^{-1}(i)| >= r + 2 - i for all i in I_0. Then J = I_0.
Proof. Put P_m := sum_{i in I_0} i^m and S_m := sum_{j in J} j^m. For m >= 0:
    sum_{j in J} (r+2-j) j^m >= sum_{j in J} sum_{l in dom lambda_j} j^m >= sum_{j in J} sum_{l in dom lambda_j} lambda_j(l)^m = sum_{i=1}^r n_i i^m >= sum_{i in I_0} (r+2-i) i^m,
using |dom lambda_j| <= r+2-j, lambda_j(l) <= j, n_i >= 0, and the hypothesis. Rearranging: P_{m+1} - S_{m+1} >= (r+2)(P_m - S_m). For m = 0 the right side is (r+2)(t - t) = 0, so P_1 >= S_1, and by induction P_m >= S_m for all m >= 1. Let i_1 := max I_0 and j_1 := max J. If j_1 > i_1 then S_m > P_m for large m (the term j_1^m dominates), contradiction; so j_1 <= i_1. Also i_1 needs n_{i_1} >= r+2-i_1 >= 1 slots from blocks j >= i_1, so j_1 >= i_1. Hence j_1 = i_1; the only block that can serve i_1 is j_1, and it has exactly r+2-j_1 = r+2-i_1 slots, so all slots of block j_1 serve i_1 and no other element of I_0 receives anything from block j_1. Removing i_1 from I_0 and j_1 from J leaves a system satisfying the same hypotheses with t - 1; by induction on t (the case t = 0 is trivial) I_0 \ {i_1} = J \ {j_1}, so I_0 = J.

**Conclusion.** The coefficient of M in h is the sum over J of the coefficients of M in red(Tr(prod_{j in J} R_j^{alpha_j})); by D6.2 only J = I_0 contributes, and by D6.1 that contribution equals the coefficient of M in g', which is nonzero by D4.1. M is nonconstant, so h is a nonconstant polynomial with individual degrees <= 2^k - 1, hence a nonconstant function on K^r (D1.1), i.e. F is nonconstant on A = H(K^r). QED Theorem D.

## D7. The concrete family: parameters and the disperser dimension

For n >= 1 (as in `proof.md`): k = k(n) is the least prime with k > 2 n^{3/5} (equivalently k^5 > 32 n^3). When n > k, the payload length is m := n - k >= 1, r := ceil(m/k) >= 1, t := ceil(sqrt(r)), and for a monic irreducible Q of degree k the field is K = F_2[Z]/(Q) (any such Q; Theorem D holds for every field with 2^k elements). A payload y in F_2^m is padded with rk - m zeros, split into r blocks of k bits, and block i is mapped to K by the F_2-linear bijection phi(c_0 ... c_{k-1}) = sum_j c_j Z^j. The map iota_Q : F_2^m -> K^r, y -> (phi(block_1), ..., phi(block_r)), is F_2-linear and injective. g_{n,Q} := F o iota_Q with F as in Theorem D; this is the restriction of f_n to the inputs whose prefix encodes Q (`proof.md`, Section 1).

**D7.1.** For every n >= 10^7: (i) k >= r(r+1)/t; (ii) kr/t + kt + r(r+1) < 8 n^{4/5}; (iii) floor(8 n^{4/5}) + 1 < m = n - k.
Proof. Bertrand's postulate (a prime in (M, 2M] for every integer M >= 1) with M := ceil(2 n^{3/5}) gives 2 n^{3/5} < k <= 2M <= 4 n^{3/5} + 2 <= 5 n^{3/5} (as n^{3/5} >= 2). Hence r <= m/k + 1 <= n/k + 1 < n^{2/5}/2 + 1, sqrt(r) < n^{1/5}/sqrt(2) + 1, and sqrt(r) <= t < sqrt(r) + 1. Write N := n^{1/5}. For (iii): 8 N^4 + 1 + k <= 8 N^4 + 5 N^3 + 1 < 25 N^4 <= N^5 = n for N >= 25. For (i) and (ii):
 - r(r+1)/t <= sqrt(r)(r+1) < (N/sqrt2 + 1)(N^2/2 + 2) = N^3/(2 sqrt2) + N^2/2 + sqrt2 N + 2, which is < 2 N^3 < k for N >= 25 (at N = 25: 5524 + 313 + 36 + 2 = 5875 < 31250), proving (i);
 - kr/t <= k sqrt(r) < 5 N^3 (N/sqrt2 + 1) = (5/sqrt2) N^4 + 5 N^3;
 - kt < 5 N^3 (N/sqrt2 + 2) = (5/sqrt2) N^4 + 10 N^3;
 - r(r+1) < (N^2/2 + 1)(N^2/2 + 2) = N^4/4 + (3/2) N^2 + 2.
Summing, kr/t + kt + r(r+1) < (10/sqrt2 + 1/4) N^4 + 15 N^3 + 1.5 N^2 + 2 < 7.3211 N^4 + 15 N^3 + 1.5 N^2 + 2. For N >= 25 (i.e. n >= 25^5 = 9765625 <= 10^7): 15 N^3 + 1.5 N^2 + 2 <= (15/25 + 1.5/625 + 2/390625) N^4 < 0.6025 N^4, so the sum is < 7.9236 N^4 < 8 N^4 = 8 n^{4/5}, proving (ii).

**Corollary D7.2.** Let d(n) := floor(8 n^{4/5}) + 1. For every n >= 10^7 and every monic irreducible Q of degree k(n), the payload function g_{n,Q} : F_2^m -> F_2 (m = n - k) is nonconstant on every affine subspace of F_2^m of dimension >= d(n), and d(n) < m. In particular g_{n,Q} is an affine disperser for dimension d(n) in the sense of `elimination.md` (nonconstant on every affine subspace of dimension exactly d(n)), and d(n)/n -> 0.
Proof. d(n) < m is D7.1(iii). Let A subset F_2^m be affine with dim A >= d(n) > 8 n^{4/5}. Then iota_Q(A) is an affine subspace of K^r = F_2^{rk} of the same dimension (iota_Q is injective and linear), which exceeds kr/t + kt + r(r+1) by D7.1(ii); k is prime and k >= r(r+1)/t by D7.1(i); 1 <= t <= r since r >= 1 and ceil(sqrt(r)) <= r. By Theorem D (applied with the field K = F_2[Z]/(Q)), F is nonconstant on iota_Q(A), i.e. g_{n,Q} = F o iota_Q is nonconstant on A. A subspace of dimension > d(n) contains one of dimension d(n), so the two formulations agree.
