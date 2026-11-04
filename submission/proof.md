# A full-B2 family bound with explicit error constants 51/4 and 25/2, diagonal shared-partner elimination, and trajectory subspace extraction

## Main theorem

There is one fixed total language $L$ in $\mathrm{P}$ such that, writing $F_N : \{0,1\}^N \to \{0,1\}$ for its $N$-bit characteristic function, every circuit over all sixteen binary Boolean truth tables, with arbitrary sharing and fan-out and one designated output, satisfies

    C_B2(F_N) >= 4*N - (51/4)*N^(4/5)    for every integer N >= 2^20,
    C_B2(F_N) >= 4*N - (25/2)*N^(4/5)    for every integer N >= 2^25.

More precisely, for every integer $N \ge 2^{20}$ the proof establishes

    C_B2(F_N) >= 4*N - E(N),

    E(N) = (43849/3600)*N^(4/5) + (192/25)*N^(3/5)
           + (25/2)*N^(2/5) + 16.

The leading coefficient of $N$ is $4$, the polynomial exponent is $1$, and both logarithmic exponents are $0$. Both structured score-source claims in `submission.json` assert this identical complete family declaration and theorem. In addition, Section V proves two auxiliary unscored structural theorems (Unscored Theorems A and B, declared as string claims in `submission.json`) that establish the diagonal shared-partner elimination law (`DIAG`) and the arbitrary-depth level-1 rank / topological linearization / Dickson constant-subspace theorem for general XOR/AND circuits.

The proof is organized in five self-contained parts:
* **Part I:** Full-$B_2$ gate lower bound $C_{B_2}(f) \ge 4n - 8(d+1)$ via the ternary-pair potential $V = s + 4I + 2q$, including the complete 4-deletion proof for `AP2` and the diagonal shared-partner reduction `DIAG` (which achieves $\Delta V \le -10$ against the $-8$ requirement).
* **Part II:** Complete proof of the Ben-Sasson--Kopparty algebraic affine-disperser theorem over $\mathbb{F}_2$ with the sharper threshold $D > r(r+3)/2 + \ell\max(r/t, t-1)$.
* **Part III:** Elementary self-contained proofs of Bertrand's prime interval $(x, 2x]$ for $x \ge 512$, existence of irreducible polynomials of every prime degree $\ell$ over $\mathbb{F}_2$, and exactness of the polynomial-time prime-degree irreducibility test.
* **Part IV:** Complete actual-length analysis at every integer $N \ge 2^{20}$ using balanced prime selection $b = \lceil (24/25) N^{3/5} \rceil$ and injective zero-padded payload blocks ($N - m = \ell$), deriving the explicit rational error bounds $(51/4)N^{4/5}$ and $(25/2)N^{4/5}$ and all asymptotic $\varepsilon$-thresholds.
* **Part V:** Complete proofs of Unscored Theorem A (multi-gadget diagonal shared-partner elimination yielding $4n + 3t - 8(d+1)$ and $4n + 4t - 8(d+1)$) and Unscored Theorem B (arbitrary-depth admissible-rank + iterative topological linearization + quadratic subspace extraction).

## One fixed total evaluator

On an input $z \in \{0,1\}^N$ of actual bit length $N$:

1. If $N < 2^{20}$, output $0$.
2. Let $b$ be the least positive integer satisfying $25^5 b^5 \ge 24^5 N^3$. Thus $b = \lceil (24/25) N^{3/5} \rceil$.
3. Let $\ell$ be the first prime strictly larger than $b$, found by deterministic trial division up to $\lfloor\sqrt{\ell}\rfloor$.
4. Put $m = N - \ell$, $r = \lceil m / \ell \rceil$, and let $t = \lceil\sqrt{r}\rceil$ be the least integer with $t^2 \ge r$, using exact integer arithmetic.
5. Read the first $\ell$ bits of $z$ as the coefficients $p_0, \dots, p_{\ell-1} \in \mathbb{F}_2$ of the monic degree-$\ell$ polynomial
       P(X) = X^\ell + sum_(j=0)^(\ell-1) p_j X^j    in F_2[X].
6. Test $P$ for irreducibility over $\mathbb{F}_2$ by checking both
       X^(2^\ell) == X  (mod P)    and    gcd(P(X), X^2 - X) = 1.
   If either check fails, output $0$. Compute the first check by $\ell$ successive modular squarings starting from $X \bmod P$, and the second by polynomial Euclid.
7. Work in the quotient ring $K = \mathbb{F}_2[X]/(P)$. Append exactly $r\ell - m$ zero bits to the remaining $m = N - \ell$ payload bits of $z$, divide the resulting $r\ell$ bits into $r$ consecutive blocks of $\ell$ bits, and interpret block $i$ ($1 \le i \le r$) in the polynomial basis $1, X, \dots, X^{\ell-1}$ with the first bit of each block as the constant coefficient:
       x_i = sum_(j=0)^(\ell-1) x_(i,j) X^j    in K.
8. Put $w_i = r + 2 - i$ and $y_i = x_i^{2^{w_i} - 1}$ for $i = 1, \dots, r$. Compute each $y_i$ by the recurrence $u_1 = x_i$ and $u_{a+1} = u_a^2 \cdot x_i$ for $a = 1, \dots, w_i - 1$ in $K$. Then compute the order-$t$ elementary symmetric polynomial $\mathrm{Sym}_t(y_1, \dots, y_r)$ by the division-free dynamic program:
       E_0 = 1,    E_1 = ... = E_t = 0;
       for i = 1, ..., r:
           for j = min(i, t) down to 1:
               E_j <- E_j + y_i * E_(j-1).
9. Output the absolute field trace
       Tr_(K/F_2)(E_t) = sum_(j=0)^(\ell-1) E_t^(2^j)    in F_2.

The symmetric-polynomial recurrence is exact because the $j$-element subsets of $\{1, \dots, i\}$ partition into those excluding $i$ (sum $E_j$ from step $i-1$) and those including $i$ (sum $y_i E_{j-1}$ from step $i-1$); updating $j$ in descending order preserves the needed values from step $i-1$.

Every length $N \ge 0$ and every input string $z \in \{0,1\}^N$ is covered by this single deterministic algorithm. Part IV proves that for every $N \ge 2^{20}$ all parameter requirements of Part II hold and the evaluator runs in fixed polynomial time in $N$. Define the uniform NP witness bound $W(N) = 0$ and the polynomial-time verifier $V(z, w)$ that accepts if and only if $w$ is the empty string and the evaluator returns $1$ on $z$. Then $z \in L$ if and only if there exists $w$ with $|w| \le W(|z|)$ such that $V(z, w) = 1$. Both directions of NP membership hold unconditionally, with no per-length advice and no witness bits supplied as circuit inputs.

## Provenance and attribution

* **Algebraic disperser construction:** Eli Ben-Sasson and Swastik Kopparty, *Affine Dispersers from Subspace Polynomials*, author preprint dated 2010-03-12, Sections 3, 6, and 7 (Theorem 7.3): `https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf` (also at `https://www.math.toronto.edu/swastik/affine-disperser.pdf`). All needed binary-field lemmas, gap properties, and monomial-isolation steps are proved in full in Part II.
* **Fair affine system and parameter role exchanges:** Om Swostik Mishra and Rahul Santhanam, *An Explicit 3.5n - o(n) Size Lower Bound against General Boolean Circuits*, ECCC TR26-260, 2026-10-19, Sections 2--4: `https://eccc.weizmann.ac.il/report/2026/260/download/`.
* **Pinned full-$B_2$ circuit model conventions:** Jiatu Li and Tianqi Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4, and 3.3.
* **Prior Yukon submissions (declared in `dependencies.json`):**
  * `f4f2c464-3045-4c32-a6e4-32478c05f343` (mooselumph): certificate-prefix BK wrapper and zero-padded payload embedding.
  * `5e9fcf7f-8338-4729-9faa-a18994661af6` (ArthurOstapenko): disjoint three-point parameter pairs and the `GP` domain restriction.
  * `bf313f35-d97b-4836-b1e8-30afd1157ac2` / commit `0fe3ef48faab5aa682a459ff164eaf1ca04cf782` (GordoAR): potential $V = s + 4I + 2q$, 4-deletion `AP2` case analysis, and sharper BK entropy selection $\delta_i = \max(d_i - w_i, 0)$, which we combine here with balanced prime selection $b = \lceil (24/25) N^{3/5} \rceil$, zero-padded payload blocks ($N - m = \ell$), and the new structural theorems in Part V.
  No prior review verdict is used as a mathematical premise; every step needed for the scored theorem and unscored lemmas is proved in this file.


# Part I. Full-B2 gate bound from ordinary affine dispersal

**Theorem I.1.** *Let $d \ge 0$ be an integer. If $f : \mathbb{F}_2^n \to \mathbb{F}_2$ is nonconstant on every affine subspace of $\mathbb{F}_2^n$ of dimension strictly greater than $d$, then every full-$B_2$ circuit $C$ computing $f$ has at least*

    C_B2(f) >= 4*n - 8*(d + 1)

*gates. This holds in the pinned full-$B_2$ model with all sixteen binary truth tables, unrestricted fan-out and sharing, repeated pins, constants, unary gates, and one designated output.*

## 1. State representation, product domain, and potential

We maintain a state consisting of:
* $M$ independent binary parameters $u = (u_1, \dots, u_M) \in \mathbb{F}_2^M$ and an affine injective map $\Phi : \mathbb{F}_2^M \to \mathbb{F}_2^n$, defined and injective on the entire cube $\mathbb{F}_2^M$.
* A collection $P$ of $q$ pairwise disjoint unordered pairs of parameters $\{p, h\}$. Each pair $\{p, h\} \in P$ has a specified excluded corner $(e_p, e_h) \in \mathbb{F}_2^2$, so its allowed set is $T_{p,h} = \mathbb{F}_2^2 \setminus \{(e_p, e_h)\}$ of size $3$. Parameters not belonging to any pair in $P$ are called *free* (or *unpaired*). The allowed domain $D \subseteq \mathbb{F}_2^M$ is the Cartesian product of $T_{p,h}$ over all $\{p, h\} \in P$ and $\mathbb{F}_2$ over all free parameters.
* A fair semi-circuit comprising:
  1. A square invertible affine system $A v = N u + c$ over $\mathbb{F}_2$ in $r$ auxiliary variables $v = (v_1, \dots, v_r)$ and the $M$ parameters $u$, where $A \in \mathbb{F}_2^{r \times r}$ is invertible and every equation (row of $(A \mid N)$) has at most $3$ nonzero entries after cancellation over $\mathbb{F}_2$.
  2. An acyclic topologically ordered list of $B_2$ gates, each reading two pins that may be constants, parameters, auxiliary variables, or earlier gates in the list; no gate output appears in the affine system.
  3. A designated output, which is a constant or a literal ($t$ or $t \oplus 1$) of a parameter, auxiliary variable, or gate.
  On every point $u \in D$, the uniquely determined output of the semi-circuit equals $f(\Phi(u))$.

Let $s$ denote the total number of retained equations plus retained acyclic gates (one unit of charge each). The *children* of a parameter or auxiliary variable $w$ are the distinct equations in which $w$ appears with coefficient $1$ and the distinct acyclic gates that read $w$ on at least one pin (a repeated pin on the same gate counts as one child; the output reference is not a child). The *degree* $\deg(w)$ is its number of children. Let $I$ denote the number of free (unpaired) parameters $x$ with $\deg(x) \ge 1$. Set

    a = M - q,        V = s + 4*I + 2*q.

At initialization, $M = a = n$, $q = 0$, $\Phi = \mathrm{id}_{\mathbb{F}_2^n}$, $D = \mathbb{F}_2^n$, the affine system is empty, the acyclic gate list is $C$, $s = |C|$, and $V_0 \le |C| + 4n$. Throughout the process, $V \ge 0$, and unused parameters ($\deg(x) = 0$) remain in $\mathbb{F}_2^M$ and $\Phi$.

**Safe affine subspaces and terminal exclusions.** For each pair $\{p, h\} \in P$ with excluded corner $(e_p, e_h)$, fixing $p$ to its *safe value* $1 \oplus e_p$ leaves both values $0, 1$ of $h$ allowed in $T_{p,h}$. Fixing one member of each of the $q$ disjoint pairs to its safe value defines an affine subspace $U_{\mathrm{safe}} \subseteq D$ of dimension $M - q = a$. Because $\Phi$ is affine and injective on $\mathbb{F}_2^M$, $\Phi(U_{\mathrm{safe}})$ is an affine subspace of $\mathbb{F}_2^n$ of dimension $a$. Since $f$ is nonconstant on every affine subspace of dimension $> d$, $f \circ \Phi$ is nonconstant on $U_{\mathrm{safe}}$ whenever $a \ge d + 1$.

We execute the elimination loop while

    a >= d + 2.

At any state with $a \ge d + 2$:
* The designated output cannot be a constant, nor an affine function of the parameters (hence cannot be a literal of a parameter or auxiliary variable): on $U_{\mathrm{safe}}$ (dimension $a \ge d + 2$), any affine Boolean function is constant on an affine hyperplane of dimension $\ge a - 1 \ge d + 1 > d$, contradicting the disperser hypothesis on $f$.
* The designated output cannot be a literal of an AND-type gate $G = (U + \alpha)(W + \beta) + \gamma$ whose two inputs $U, W$ are affine functions of the parameters: on $U_{\mathrm{safe}}$, if either factor $U + \alpha$ or $W + \beta$ is nonconstant, setting that factor to $0$ restricts $U_{\mathrm{safe}}$ to an affine hyperplane of dimension $\ge a - 1 \ge d + 1 > d$ on which $G \equiv \gamma$ is constant; if both factors are constant on $U_{\mathrm{safe}}$, $G$ is already constant on $U_{\mathrm{safe}}$. Both cases contradict the disperser hypothesis.

## 2. Legal coordinate operations

Every coordinate operation below preserves the affine injectivity of $\Phi$ and ensures that the new product domain maps into the old domain $D$:
1. **Affine role exchange (`AR`):** Let $Q$ be an auxiliary variable whose solved affine expression in $v = A^{-1} N u + A^{-1} c$ has coefficient $1$ on a free parameter $x$, so $Q = x \oplus \ell(u_{-x})$ for an affine function $\ell$ of the remaining parameters. Substitute $x \leftarrow w \oplus \ell(u_{-x})$ in $\Phi$, rename $Q$ as the free parameter $w$, and rename $x$ as an auxiliary variable. Because $x$ was a free coordinate, for every fixed assignment to $u_{-x}$ the map $w \mapsto w \oplus \ell(u_{-x})$ is a bijection of $\mathbb{F}_2$, so the product domain $D$ is preserved (with $w$ replacing $x$) even if $\ell$ mentions paired coordinates. In the affine system, the column of $Q$ in $A$ and the column $N_x$ of $x$ in $N$ are swapped; no equation's terms change. Since $(A^{-1} N_x)_Q = 1$, replacing the $Q$-th column of $A$ by $N_x$ multiplies $\det(A)$ by $(A^{-1} N_x)_Q = 1$, so the new auxiliary matrix remains invertible over $\mathbb{F}_2$.
2. **Free parameter fixing (`UF`):** Fix a free parameter $x$ to a constant $c_0 \in \mathbb{F}_2$, or identify two distinct free parameters $x \leftarrow y \oplus c_0$. This composes $\Phi$ with an injective affine map on the remaining $M - 1$ coordinates, leaves $q$ unchanged, and gives $\Delta M = \Delta a = -1$.
3. **Safe fixing of a paired parameter (`PF-safe`):** For $\{p, h\} \in P$ with excluded corner $(e_p, e_h)$, fix $p = 1 \oplus e_p$ and remove the pair $\{p, h\}$ from $P$, making $h$ a free parameter. This gives $\Delta M = \Delta q = -1$ and $\Delta a = 0$.
4. **Pair-line restriction (`PF-line`):** Restrict a pair $\{p, h\} \in P$ to the affine line through any two of the three allowed points in $T_{p,h}$ (which either fixes $p = 1 \oplus e_p$, fixes $h = 1 \oplus e_h$, or imposes $p = h \oplus (e_p \oplus e_h \oplus 1)$) and remove $\{p, h\}$ from $P$, leaving one surviving free coordinate. This gives $\Delta M = \Delta q = -1$ and $\Delta a = 0$.
5. **Gate pairing (`GP`):** For two distinct free parameters $x, y$, declare $\{x, y\}$ to be a new pair in $P$ with a chosen excluded corner $(e_x, e_y) \in \mathbb{F}_2^2$, restricting $D$ accordingly while leaving $\Phi$ and $M$ unchanged. This gives $\Delta M = 0$, $\Delta q = +1$, and $\Delta a = -1$.

## 3. Affine movement and elementary deletions

Up to input and output negations, the sixteen binary truth tables consist of $2$ constant tables, $4$ unary tables (depending on exactly one input), $2$ XOR-type tables ($t_1 \oplus t_2 \oplus c$), and $8$ AND-type tables ($(t_1 \oplus \alpha)(t_2 \oplus \beta) \oplus \gamma$).

**Move (`M`):** If an acyclic gate $J$ computes an affine operation $b_1 t_1 \oplus b_2 t_2 \oplus c_0$ ($b_1, b_2, c_0 \in \mathbb{F}_2$) and its inputs $t_1, t_2$ are constants, parameters, or auxiliary variables, remove $J$ from the acyclic gate list, regard $J$ as a new auxiliary variable, and add the equation $J \oplus b_1 t_1 \oplus b_2 t_2 = c_0$ to the affine system. The equation has at most $3$ nonconstant terms, and the new auxiliary matrix is block lower-triangular with a new diagonal entry $1$, hence invertible. Because one gate is replaced by one equation and no parameter gains a new child, $\Delta M = \Delta q = \Delta a = \Delta s = 0$ and $\Delta I \le 0$, so $\Delta V \le 0$.

**Elementary deletions (`D1`--`D4` and `D2'`):** Each deletion below has $\Delta M = \Delta q = \Delta a = 0$, $\Delta s = -1$, and $\Delta I \le 0$, hence $\Delta V \le -1$.
* **`D1` (Degenerate or dependent gate):** If an acyclic gate $J$ is constant or unary, has a constant pin or two identical pins, or is an AND-type gate whose two inputs are affine values (parameters or auxiliaries) with linearly dependent linear parts over $\mathbb{F}_2$, then $J$ computes a constant or a literal of one of its inputs. Reroute all occurrences of $J$ to that constant or literal and delete $J$. Any parameter that gains a child from rerouting was already a parent of $J$ and thus already counted in $I$, so $\Delta I \le 0$.
* **`D2` (Childless gate):** Delete any acyclic gate $J$ that has no children and is not the designated output.
* **`D3` (Short equation):** Any equation $E$ with at most $2$ nonconstant terms must contain at least one auxiliary variable $v$ (otherwise the corresponding row of $A$ would be zero, contradicting invertibility of $A$). Write $E$ as $v = U \oplus c_0$, where $U$ is $0$ or one other variable (parameter or auxiliary). Substitute $U \oplus c_0$ for every occurrence of $v$ in the remaining equations, gates, and output, and delete $v$ and $E$. Adding row $E$ to the other rows containing $v$ is an elementary row operation on $A$ preserving $\det(A)$, after which column $v$ of $A$ has a single $1$ in row $E$; Laplace expansion along column $v$ shows the remaining $(r-1)\times(r-1)$ auxiliary matrix is invertible. Equation widths never increase ($U$ is at most one term), and any parameter gaining occurrences was already in $E$, so $\Delta I \le 0$.
* **`D4` (Isolated free parameter in an equation) and `D2'` (Isolated auxiliary in an equation):** If a free parameter $x$ has a single child and that child is an equation $E$, the column $N_x$ has a single $1$ (in row $E$), so $A^{-1} N_x \ne 0$. Choose an auxiliary $Q$ with $(A^{-1} N_x)_Q = 1$ and apply `AR` to exchange $x$ and $Q$. Now $x$ is an auxiliary variable whose only occurrence in the entire representation is in equation $E$ (it cannot be the designated output at $a \ge d+2$). This reduces `D4` to **`D2'`**: whenever an auxiliary variable $v$ has a single child and that child is an equation $E$ (and $v$ is not the output), column $v$ of $A$ has a single $1$ in row $E$. Delete $v$ and $E$ without substitution; Laplace expansion along column $v$ shows the remaining auxiliary matrix is invertible, and for every assignment to the surviving variables row $E$ uniquely determines the deleted $v$. In `D4`, the removed used free parameter $x$ is replaced by at most one used free parameter $Q$, so $\Delta I \le 0$ and $\Delta s = -1$.

**Batch deletion counting rule.** Whenever a parameter $w$ is fixed to a constant, every original gate child of $w$ acquires a constant pin (`D1`) and every original equation child of $w$ is shortened to at most $2$ nonconstant terms (`D3`). Subsequent `D1` and `D3` literal substitutions never lengthen an equation and never restore a nonconstant value to a constant pin, and invertibility of $A$ prevents any surviving equation from losing all auxiliary variables before it is deleted. Therefore, for any set $S_0$ of distinct original gates and equations that each acquire a constant pin or are shortened to $\le 2$ terms, all $|S_0|$ objects are guaranteed to be deleted; if one is removed earlier during cascading propagation, that removal is its single counted deletion.

## 4. Preliminary pair reductions and diagonal shared-partner reduction

At the start of each step (with $a \ge d+2$), exhaust `M`, `D1`--`D4`, and `D2'` until none applies. Then check the following reductions in order, restarting the loop after any reduction fires.

* **`P0` (Unused paired parameter):** If $p \in P$ has $\deg(p) = 0$, safely fix its partner $h$ (`PF-safe`) and remove $\{p, h\}$. The newly free $p$ still has $\deg(p) = 0$, so $\Delta M = \Delta q = -1$, $\Delta a = 0$, $\Delta I \le 0$, $\Delta s \le 0$, and $\Delta V \le -2$.
* **`P2` (Paired parameter of degree $\ge 2$):** If $p \in P$ has $\deg(p) \ge 2$, safely fix $p$ (`PF-safe`), remove $\{p, h\}$, and delete at least $2$ distinct children of $p$ by `D1`/`D3`. Only the partner $h$ becomes newly free and can enter $I$, so $\Delta M = \Delta q = -1$, $\Delta a = 0$, $\Delta s \le -2$, $\Delta I \le +1$, and $\Delta V \le -2 + 4 - 2 = 0$.
* **`PC` (Partners sharing their unique child):** After `P0` and `P2` are absent, every paired parameter has degree $1$. If $\{p, h\} \in P$ share their unique child $B$ (either a gate with pins $p, h$ or a $3$-term equation $v \oplus p \oplus h = c_0$), by pigeonhole at least two of the three allowed points of $T_{p,h}$ give the same Boolean value to the gate $B$ (or to $p \oplus h$). Restrict $\{p, h\}$ to the affine line through those two allowed points (`PF-line`) and delete $B$ (`D1` or `D3`), rerouting $B$ (or $v$) to its constant value. Because $B$ was the unique child of both $p$ and $h$, the surviving coordinate of the line has degree $0$, so $\Delta M = \Delta q = -1$, $\Delta a = 0$, $\Delta s \le -1$, $\Delta I \le 0$, and $\Delta V \le -1 - 2 = -3$.
* **`PE` (Paired parameter in an equation):** If $p \in P$ has unique child equation $E$, its partner $h$ has a distinct unique child $B \ne E$ (by `PC`). Safely fix $h$ (`PF-safe`) and delete $B$ (`D1` or `D3`). Because $\deg(p) = 1$ with child $E \ne B$, $B$ does not contain $p$, so rerouting $B$ cannot introduce $p$ elsewhere. If $E$ was already removed during the deletion of $B$, that supplies the second deletion; otherwise apply `D4` to $(p, E)$. In both cases at least $2$ distinct objects ($B$ and $E$) are deleted and at most one new free parameter enters $I$, giving $\Delta M = \Delta q = -1$, $\Delta a = 0$, $\Delta s \le -2$, $\Delta I \le +1$, and $\Delta V \le -2 + 4 - 2 = 0$.

When `M`, `D1`--`D4`, `D2'`, and `P0, P2, PC, PE` are all inapplicable, every paired parameter has degree $1$ and its unique child is an acyclic gate, so **no paired parameter appears in the affine system** ($N$ has zero columns on all paired parameters).

## 5. Minimal AND preparation and elimination (`GP`, `AP1`, `AP2`, and `DIAG`)

Since $a \ge d+2$ and the output cannot be affine, at least one acyclic gate remains. Let $G$ be the first gate in topological order. By `M` and `D1`, $G = (U \oplus \alpha)(W \oplus \beta) \oplus \gamma$ is an AND-type gate whose two inputs $U, W$ are affine values with linearly independent linear parts over $\mathbb{F}_2$, $G$ is not the output, and $G$ has $k = \deg(G) \ge 1$ distinct gate children.

**Optional Diagonal Shared-Partner Reduction (`DIAG`):** Suppose $U$ and $W$ are distinct auxiliary variables of degree $2$, appearing only in $G$ and in $3$-term equations $E_U : U \oplus s \oplus b_1 = c_1$ and $E_W : W \oplus s \oplus b_2 = c_2$ respectively, where $s, b_1, b_2$ are three distinct free parameters and $s$ has degree $2$, appearing only in $E_U$ and $E_W$ (if this configuration is not present, skip `DIAG` and proceed to standard `AR` preparation below). Impose the free-parameter identification (`UF`, $\Delta a = -1$)
    b_1 <- b_2 + (c_1 + c_2 + alpha + beta + 1).
Then in $E_U$ and $E_W$, $(U \oplus \alpha) \oplus (W \oplus \beta) = b_1 \oplus b_2 \oplus c_1 \oplus c_2 \oplus \alpha \oplus \beta = 1$, so $(U \oplus \alpha) = (W \oplus \beta) \oplus 1$ identically for all $s, b_2 \in \mathbb{F}_2$. Because $z(z \oplus 1) = 0$ over $\mathbb{F}_2$, the minimal AND gate $G = (U \oplus \alpha)(W \oplus \beta) \oplus \gamma \equiv \gamma$ computes the constant $\gamma$. Replace $G$ by $\gamma$ and delete $G$ (`D1`), then immediately apply `D2'` to delete the now-degree-$1$ auxiliaries $(U, E_U)$ and $(W, E_W)$, and finally delete all $k \ge 1$ gate children of $G$ by `D1` (as each acquired a constant pin from $G$). Since $E_U$ and $E_W$ were the only children of $s$, the free parameter $s$ becomes unused ($\deg'(s) = 0$) in addition to $b_1$ being eliminated by the identification. Thus $\Delta M = \Delta a = -1$, $\Delta q = 0$, $\Delta s \le -(3 + k) \le -4$, and two used free parameters ($b_1$ and $s$) leave $I$ ($\Delta I \le -2$), giving
    Delta V <= -4 - 8 = -12 <= -8 = 8 * Delta a.

**Standard `AR` preparation of $G$:** If $G$ has an auxiliary parent $Q$, its solved form $Q = \sum_x \lambda_x x \oplus c_0$ depends only on free parameters (since paired columns of $N$ are zero). Because the two parents of $G$ have independent linear parts, there exists a free parameter $x$ with $\lambda_x = 1$ (and $x$ distinct from the other parent of $G$ if that parent is already a free parameter). Exchange $Q$ and $x$ by `AR` ($\Delta V = 0$), recheck `M`, deletions, and pair rules (ending the step if any reduction fires), and repeat at most once for the second parent of $G$. Then both parents of $G$ are distinct parameters.

Now eliminate $G$ according to the number of paired parents of $G$:

1. **`GP` (Both parents $x, y$ are free parameters):**
   An AND-type gate $G = (x \oplus \alpha)(y \oplus \beta) \oplus \gamma$ takes the constant value $\gamma$ on the $3$-point set $T_{x,y} = \mathbb{F}_2^2 \setminus \{(1 \oplus \alpha, 1 \oplus \beta)\}$. Add $\{x, y\}$ to $P$ with excluded corner $(1 \oplus \alpha, 1 \oplus \beta)$. On the restricted domain, $G \equiv \gamma$ is constant; delete $G$ and all $k \ge 1$ distinct gate children of $G$ by `D1`. Both $x$ and $y$ were used free parameters and are now paired, and no unused parameter gains a child, so $\Delta M = 0$, $\Delta q = +1$, $\Delta a = -1$, $\Delta s \le -(1 + k)$, $\Delta I \le -2$, and
       Delta V <= -(1 + k) - 8 + 2 = -7 - k <= -8.

2. **`AP1` (One parent $p$ is paired with partner $h$, and one parent $x$ is free):**
   Let $B \ne G$ be the unique gate child of $h$ (`PC` ensures $B \ne G$). Safely fix $h$ (`PF-safe`, removing $\{p, h\}$) and fix the free parameter $x = \beta$ (`UF`) so that $G$ becomes constant. Now $p$ is free and unused ($\deg'(p) = 0$ after $G$ is deleted), $x$ is removed ($\Delta I \le -1$), and $G, B$ are deleted. For a third distinct gate deletion: if $G$ has a child $J \ne B$, $J$ acquires a constant pin from $G$ and is deleted by `D1`; otherwise $B$ is the sole child of $G$, so $B$ reads both $G$ (constant) and $h$ (constant), hence $B$ becomes constant, cannot be the output (since $a' = a - 1 \ge d + 1 > d$), and by `D2` closure has a child $J \notin \{G, B\}$ (by acyclicity) that is deleted by `D1`. Thus $\Delta M = -2$, $\Delta q = -1$, $\Delta a = -1$, $\Delta s \le -3$, $\Delta I \le -1$, and
       Delta V <= -3 - 4 - 2 = -9 <= -8.

3. **`AP2` (Both parents $p, r$ are paired, with partners $h, k$):**
   By `PC`, $\{p, h\}$ and $\{r, k\}$ are distinct pairs, and the unique gate children $B_h$ of $h$ and $B_k$ of $k$ satisfy $B_h \ne G$ and $B_k \ne G$. Safely fix $h$ and $k$ (`PF-safe`, removing both pairs), then fix the now-free $p = \alpha$ (`UF`) to make $G$ constant. The surviving free $r$ had only $G$ as child and becomes unused ($\deg'(r) = 0$), so $\Delta I \le 0$, $\Delta M = -3$, $\Delta q = -2$, and $\Delta a = -1$. We prove that at least **four** distinct gates are deleted:
   * **Case 1 ($B_h \ne B_k$):** The three gates $S_0 = \{G, B_h, B_k\}$ are distinct and deletable. Let $S_{\mathrm{const}} \subseteq S_0$ be the subset of gates in $S_0$ that become constant functions under the restriction; $G \in S_{\mathrm{const}}$, so $S_{\mathrm{const}} \ne \emptyset$. Let $M_*$ be the topologically last gate in $S_{\mathrm{const}}$. Since $a' = a - 1 \ge d + 1 > d$, $M_*$ cannot be the designated output, so by `D2` closure $M_*$ has a gate child $J$. Because $M_*$ is constant, $J$ acquires a constant pin and is deletable by `D1`. Could $J \in S_0$? By acyclicity and minimality of $G$, $J \ne G$ and $J \ne M_*$. If $J = B_h$ (or $B_k$), then since $B_h$ already has parameter pin $h$ (fixed to a constant) and its other pin is $M_*$ (constant), $B_h$ itself becomes constant, so $B_h \in S_{\mathrm{const}}$ and $B_h$ is topologically after $M_*$, contradicting the choice of $M_*$! Hence $J \notin S_0$, giving $4$ distinct deleted gates $\{G, B_h, B_k, J\}$.
   * **Case 2 ($B_h = B_k = B$):** Here $B$ has pins $h$ and $k$, both fixed to constants, so **both** $G$ and $B$ become constant ($G \ne B$, and neither is a child of the other since $G$ reads $p, r$ and $B$ reads $h, k$). Neither $G$ nor $B$ can be the output, so each has at least one child. If we can choose $K_G \in \mathrm{ch}(G)$ and $K_B \in \mathrm{ch}(B)$ with $K_G \ne K_B$, then $\{G, B, K_G, K_B\}$ are $4$ distinct deleted gates. Otherwise $\mathrm{ch}(G) = \mathrm{ch}(B) = \{K\}$ for a single gate $K$; then the two pins of $K$ are $G$ and $B$, both constant, so $K$ becomes constant, is not the output, and has a child $J \notin \{G, B, K\}$, giving $4$ distinct deleted gates $\{G, B, K, J\}$.
   Therefore $\Delta s \le -4$ in all cases, and
       Delta V <= -4 + 0 - 4 = -8.

## 6. Termination and telescoping

Every step preserves all state invariants, satisfies $\Delta a \in \{0, -1\}$ and $\Delta V \le 8\,\Delta a$, never increases $M$ or $s$, and strictly decreases $M + s$ after at most finitely many `M` and `AR` preparations (each `M` removes an acyclic gate, and at most two `AR` exchanges occur for the fixed minimal AND $G$ before a reducing rule fires). Thus the algorithm terminates after $T \le n + |C|$ steps at a state with $a_T \le d + 1$. Summing $\Delta V \le 8\,\Delta a$ gives
    |C| + 4*n >= V_0 >= V_T + 8*(n - a_T) >= 8*(n - d - 1),
and hence $|C| \ge 4n - 8(d + 1)$. This proves Theorem I.1.


# Part II. Algebraic disperser with a sharper threshold

## 1. Precise theorem

**Theorem II.1.** *Let $\ell$ be an odd prime, let $K$ be any field of size $2^\ell$, and let $r, t$ be integers satisfying $1 \le t \le r$ and $r + 1 < \ell$. Put*

    w_i = r + 2 - i,    alpha_i = 2^(w_i) - 1    (i = 1, ..., r),
    W = sum_(i=1)^r w_i = r*(r + 3)/2.

*Define the Boolean function $f : K^r \to \mathbb{F}_2$ by*

    f(x_1, ..., x_r) = Tr_(K/F_2)( sum_(I subset [r], |I|=t) product_(i in I) x_i^(alpha_i) ),

*where $\mathrm{Tr}_{K/\mathbb{F}_2}(z) = \sum_{j=0}^{\ell-1} z^{2^j}$. Under any $\mathbb{F}_2$-linear identification $K^r \cong \mathbb{F}_2^{r\ell}$, $f$ is nonconstant on every affine subspace $A \subseteq K^r$ of $\mathbb{F}_2$-dimension*

    D > W + ell * max(r/t, t - 1).                              (T)

*In particular, for $t = \lceil\sqrt{r}\rceil$, both $r/t \le \sqrt{r}$ and $t - 1 < \sqrt{r}$, so*

    D > r*(r + 3)/2 + ell * sqrt(r)                             (T')

*suffices.*

## 2. Elementary finite-field and binary-weight preliminaries

Let $q = 2^\ell$. Reducing each variable modulo $X_i^q - X_i$ yields a unique polynomial representative of degree $< q$ in each variable $X_i$ for every function $K^r \to K$ (by induction on $r$ using the univariate root bound). Consequently, a nonzero coefficient on a positive-degree monomial in this reduced representative certifies that the function is nonconstant.

Every linearized polynomial $\sum_{j=0}^{\ell-1} a_j X^{2^j}$ ($a_j \in K$) defines an $\mathbb{F}_2$-linear map $K \to K$. Since there are $q^\ell = 2^{\ell^2}$ such polynomials, the map from linearized polynomials to $\mathbb{F}_2$-linear maps $K \to K$ is injective by the root bound (degree $\le 2^{\ell-1} < q$), and $\dim_{\mathbb{F}_2}\mathrm{Hom}_{\mathbb{F}_2}(K, K) = \ell^2$, every $\mathbb{F}_2$-linear map $K \to K$ is uniquely represented by a linearized polynomial of degree $< q$.

The absolute trace $\mathrm{Tr}(z) = \sum_{j=0}^{\ell-1} z^{2^j}$ satisfies $\mathrm{Tr}(z)^2 = \mathrm{Tr}(z)$ (since $z^q = z$), hence maps $K \to \mathbb{F}_2$, and for odd $\ell$ satisfies $\mathrm{Tr}(1) = \ell \bmod 2 = 1$, so $\mathrm{Tr}$ is nonzero.

For an exponent $b \in \{0, 1, \dots, q-1\}$, let $\mathrm{wt}(b)$ denote the number of $1$s in its binary expansion. When multiplying powers of two and reducing modulo $X^q - X$ (which replaces $X^{2^\ell}$ by $X^{2^0}$), power positions are reduced modulo $\ell$, and any carry $2^j + 2^j = 2^{(j+1)\bmod \ell}$ merges two units into one unit. Thus a product of $h$ power-of-two monomials $X^{2^{j_1}} \cdots X^{2^{j_h}}$ reduces modulo $X^q - X$ to a monomial $X^b$ with $\mathrm{wt}(b) \le h$, and equality $\mathrm{wt}(b) = h$ holds (when $h < \ell$) if and only if the $h$ residues $j_1, \dots, j_h \bmod \ell$ are pairwise distinct.

## 3. Image-subspace polynomials and the prime-degree gap lemma

For the zero subspace $V = \{0\}$, take $Q_{\{0\}}(X) = 0$. For an $\mathbb{F}_2$-linear subspace $V \subseteq K$ of dimension $d \in \{1, \dots, \ell\}$, define the subspace vanishing polynomial
    P_V(X) = product_(v in V) (X + v).
By induction on $d$ (using $P_{U \oplus \mathbb{F}_2 b}(X) = P_U(X) P_U(X+b) = P_U(X)^2 + P_U(b) P_U(X)$), $P_V(X)$ is a monic linearized polynomial of degree $2^d$ with kernel $V$. Its image $W_0 = P_V(K)$ is an $\mathbb{F}_2$-linear subspace of $K$ of dimension $\ell - d$. Set $Q_V(X) = P_{W_0}(X)$. Then $Q_V(P_V(X))$ is monic of degree $2^\ell = q$ and vanishes on all of $K$, so $Q_V(P_V(X)) = X^q + X$. Moreover, composing $P_V(Q_V(X)) + (X^q + X)$ on the right with $P_V(X)$ gives $P_V(Q_V(P_V(X))) + P_V(X)^q + P_V(X) = P_V(X^q + X) + P_V(X)^q + P_V(X) = 0$, and since $P_V$ is nonconstant, $P_V(Q_V(X)) = X^q + X$ as a polynomial. In particular, $Q_V(K) \subseteq \ker(P_V) = V$, and since $\dim_{\mathbb{F}_2} \ker(Q_V) = \dim_{\mathbb{F}_2} W_0 = \ell - d$, we have $Q_V(K) = V$. Write
    Q_V(X) = sum_(j=0)^(\ell-d) a_j X^(2^j),
and extend the coefficients $a_j$ cyclically modulo $\ell$ by setting $a_j = 0$ for $\ell - d < j < \ell$. Notice $a_0 \ne 0$ (since $Q_V$ is squarefree, having $2^{\ell-d}$ distinct roots $W_0$) and $a_{\ell-d} = 1 \ne 0$.

**Lemma II.2 (Gap lemma for prime $\ell$).** *If $d \ge 2$ and $\ell$ is prime, the cyclic coefficient sequence $(a_j)_{j \in \mathbb{Z}/\ell\mathbb{Z}}$ of $Q_V$ has:*
1. *no $d$ consecutive zero entries; and*
2. *no $d - 1$ consecutive zero entries other than the designated tail $\ell - d + 1, \dots, \ell - 1$.*

*Proof.* (1) If $d$ consecutive coefficients were zero, raising $Q_V(X)$ to an appropriate Frobenius power $2^k$ and reducing modulo $X^q - X$ would cyclically shift those $d$ zeros into the top positions $\ell - d, \dots, \ell - 1$, producing a nonzero polynomial of degree $\le 2^{\ell-d-1}$ that still vanishes on all $2^{\ell-d}$ elements of $W_0$, contradicting the root bound.
(2) Suppose another cyclic run of $d - 1$ consecutive zeros existed; shifting it into the tail $\ell - d + 1, \dots, \ell - 1$ by a Frobenius power $2^k$ with $k \not\equiv 0 \pmod{\ell}$ yields $Q'(X) = Q_V(X)^{2^k} \bmod (X^q - X)$ of degree $\le 2^{\ell-d}$ with the exact same root set $W_0$ as $Q_V$. Comparing leading terms and applying the root bound to $Q' - a'_{\ell-d} Q_V$ gives $Q'(X) = c\, Q_V(X)$ for a nonzero constant $c \in K^\times$. Thus every $y \in Q_V(K) = V$ satisfies $y^{2^k} = c y$. Choosing a nonzero $y_0 \in V$, every ratio $z = y / y_0$ ($y \in V$) satisfies $z^{2^k} = z$ as well as $z^{2^\ell} = z$. Since $\ell$ is prime and $k \not\equiv 0 \pmod{\ell}$, $\gcd(k, \ell) = 1$; writing $1 = u k + v \ell$ ($u, v \in \mathbb{Z}$) shows that the $1$-st power of the Frobenius automorphism $z \mapsto z^2$ fixes $z$, so $z \in \mathbb{F}_2$ and $V \subseteq y_0 \mathbb{F}_2$ has dimension $\le 1$, contradicting $d \ge 2$. $\blacksquare$

## 4. Independent-product coefficient lemma

**Lemma II.3.** *Let $V_1, \dots, V_s \subseteq K$ be $\mathbb{F}_2$-linear subspaces of dimensions $d_i$ with $2 \le w_i < d_i \le \ell$, and let $Q_i = Q_{V_i}$. Put $e_i = w_i - 1$, $h_i = d_i - w_i \ge 1$, and $b_i = \ell - d_i$. If $\sum_{i=1}^s h_i \ge \ell$, then for every nonzero $\mu \in K$ the reduced polynomial*
    R(X_1, ..., X_s) = Tr_(K/F_2)( mu * product_(i=1)^s Q_i(X_i)^(2^(w_i) - 1) )
*has a nonzero coefficient on a monomial $\prod_{i=1}^s X_i^{\beta_i}$ with $\mathrm{wt}(\beta_i) = w_i$ for every $i = 1, \dots, s$.*

*Proof.* Write $Q_i(X) = \sum_{j=0}^{\ell-1} a_{i,j} X^{2^j}$ with indices modulo $\ell$. For a subset $B_i \subset \mathbb{Z}/\ell\mathbb{Z}$ of size $w_i$ and $\beta_i = \sum_{b \in B_i} 2^b$, define the $\ell \times \ell$ matrix $M_i(u, b) = a_{i, b-u}^{2^u}$. Expanding $Q_i(X_i)^{2^{w_i}-1} = \prod_{u=0}^{e_i} Q_i(X_i)^{2^u}$ and using the no-carry characterization of weight-$w_i$ exponents from Section 2, the coefficient of $\prod_{i=1}^s X_i^{\beta_i}$ in $R(X_1, \dots, X_s)$ is
    sum_(tau=0)^(\ell-1) mu^(2^tau) * product_(i=1)^s Perm( M_i([0, e_i], B_i - tau) )^(2^tau).     (C)
For each $i$ and any anchor $z \in \{0, \dots, b_i\}$:
* If $a_{i,z} \ne 0$, choose the ordered column tuple $B_i = (z, b_i + 1, \dots, b_i + e_i)$. The minor $M_i([0, e_i], B_i)$ is lower triangular (entries above the diagonal have coefficient indices in the zero tail $\{b_i + 1, \dots, \ell - 1\}$) with nonzero diagonal $(a_{i,z}, a_{i,b_i}^{2^1}, \dots, a_{i,b_i}^{2^{e_i}})$.
* If $a_{i,z} = 0$, let $u < z < v$ be the surrounding nonzero coefficient indices of $Q_i$ (which exist since $a_{i,0} \ne 0$ and $a_{i,b_i} \ne 0$). By Lemma II.2, the interior zero gap $g = v - u - 1$ satisfies $g \le d_i - 2$. Set $a = \min(g, e_i)$ and choose the ordered tuple of $w_i$ distinct columns
      B_i = (v, u + 1, ..., u + a, b_i + a + 1, ..., b_i + e_i).
  The minor $M_i([0, e_i], B_i)$ is again lower triangular with nonzero diagonal $(a_{i,v}, a_{i,u}^{2^1}, \dots, a_{i,u}^{2^a}, a_{i,b_i}^{2^{a+1}}, \dots, a_{i,b_i}^{2^{e_i}})$, as entries above the diagonal in the first $a + 1$ columns lie in the zero gap $(u, v)$ and those in the remaining columns lie in the zero tail $(b_i, \ell - 1)$.

Hence at $\tau = 0$, $\mathrm{Perm}(M_i([0, e_i], B_i)) \ne 0$ for every $i$. Moreover, any selected column $k \in B_i$ becomes an all-zero column in $M_i([0, e_i], B_i - \tau)$ whenever $\tau \in [k + 1, k + h_i] \pmod{\ell}$, because for $\tau = k + c$ ($1 \le c \le h_i$) and $0 \le u \le e_i$, the index $k - \tau - u = -(c + u) \equiv \ell - (c + u) \pmod{\ell}$ lies in $[\ell - d_i + 1, \ell - 1] = [b_i + 1, \ell - 1]$. For $B_i$ anchored at $z$, the union of $[k + 1, k + h_i]$ over $k \in B_i$ contains $[z + 1, z + h_i]$: if $a_{i,z} \ne 0$ or $g < e_i$, then $z \in B_i$; if $g \ge e_i$, the columns $u + 1, \dots, u + e_i$ cover $[u + 2, u + e_i + h_i] = [u + 2, u + d_i - 1]$, and since $v = u + g + 1 \le u + d_i - 1$ by Lemma II.2, adding $[v + 1, v + h_i]$ covers $[u + 2, v + h_i] \supseteq [z + 1, z + h_i]$.

Now reorder the subspaces so $d_1 = \max_i d_i$ ($b_1 = \ell - d_1$), let $s_0$ be the first index with $\sum_{i=1}^{s_0} h_i \ge b_1 + 1$, and set $z_i = \sum_{j < i} h_j \le b_1 \le b_i$ for $1 \le i \le s_0$ (and $z_i = 0$ for $i > s_0$). The intervals $[z_i + 1, z_i + h_i]$ for $1 \le i \le s_0$ cover $[1, b_1 + 1]$, while the tail columns $b_1 + 1, \dots, b_1 + e_1$ of $B_1$ additionally cover $[b_1 + 2, b_1 + e_1 + h_1] = [b_1 + 2, \ell - 1]$. Thus every nonzero shift $\tau \in \{1, \dots, \ell - 1\}$ zeroes at least one factor in (C), leaving only the nonzero $\tau = 0$ term. $\blacksquare$

## 5. Affine triangularization, weighted packet isolation, and completion of Theorem II.1

Let $A = \gamma + V \subseteq K^r$ be an affine subspace of $\mathbb{F}_2$-dimension $D$. Filtering $V$ by the subspaces whose first $0, 1, \dots, r$ block coordinates vanish yields $\mathbb{F}_2$-linear subspaces $Y_i \subseteq K$ of dimensions $d_i = \dim_{\mathbb{F}_2} Y_i$ with $\sum_{i=1}^r d_i = D$ and $\mathbb{F}_2$-linear maps $\sigma_{i,j} : Y_j \to K$ ($j < i$) such that $A = \{( \gamma_i + y_i + \sum_{j < i} \sigma_{i,j}(y_j) )_{i=1}^r : y_i \in Y_i\}$. Letting $Q_i = Q_{Y_i} : K \to Y_i$ be the image polynomial from Section 3 and $Q_{i,j}(X_j) = \sigma_{i,j}(Q_j(X_j))$, the polynomial map $R : K^r \to K^r$ with coordinates
    R_i(X_1, ..., X_r) = Q_i(X_i) + sum_(j < i) Q_(i,j)(X_j) + gamma_i
surjects onto $A$.

Put $\delta_i = \max(d_i - w_i, 0)$, so $\sum_{i=1}^r \delta_i \ge D - W$. Under hypothesis (T), $\sum_{i=1}^r \delta_i > \ell \max(r/t, t - 1)$. Since each $\delta_i < \ell$, strictly more than $t - 1$ indices have $\delta_i > 0$. Choose a $t$-element subset $I_0 \subset [r]$ of indices with the $t$ largest values of $\delta_i$; then $\delta_i \ge 1$ (so $d_i > w_i \ge 2$) for all $i \in I_0$, and by averaging,
    sum_(i in I_0) (d_i - w_i) = sum_(i in I_0) delta_i >= (t/r) * sum_(i=1)^r delta_i > ell.
By Lemma II.3 (with $\mu = 1$), $\mathrm{Tr}(\prod_{i \in I_0} Q_i(X_i)^{2^{w_i}-1})$ has a nonzero coefficient on a reduced monomial $M_* = \prod_{i \in I_0} X_i^{\beta_i}$ with $\mathrm{wt}(\beta_i) = w_i$ for each $i \in I_0$ (and degree $0$ in $X_i$ for $i \notin I_0$).

Now consider any $t$-element subset $J \subset [r]$ in the expansion of $f(R_1(X), \dots, R_r(X))$. Each factor $R_j(X)^{2^{w_j}-1}$ ($j \in J$) is a product of $w_j$ linearized factors $R_j(X)^{2^u}$, each of which is a sum of power-of-two terms in variables $X_i$ with $i \le j$ and a constant. For a product of choices from $J$ to reduce modulo $(X_i^q - X_i)$ to $M_*$, let $n_{j,i} \ge 0$ denote the number of factors from $j \in J$ that choose a term in $X_i$ ($i \in I_0$). Then $n_{j,i} > 0 \implies i \le j$ (hence $w_i \ge w_j$, with strict inequality $w_i > w_j$ whenever $i < j$), $\sum_{i \in I_0} n_{j,i} \le w_j$ for each $j \in J$, and by binary-weight subadditivity (Section 2), $\sum_{j \in J} n_{j,i} \ge \mathrm{wt}(\beta_i) = w_i$ for each $i \in I_0$. Summing $n_{j,i} / w_i$ gives
    t = sum_(i in I_0) 1 <= sum_(i in I_0, j in J) n_(j,i) / w_i <= sum_(i in I_0, j in J) n_(j,i) / w_j <= sum_(j in J) 1 = t.
Hence every inequality is an equality: $n_{j,i} = 0$ whenever $i < j$, so every factor for $X_i$ ($i \in I_0$) must come from $j = i$, forcing $I_0 \subseteq J$ and thus $J = I_0$, with all $w_i$ factors of $R_i$ choosing $Q_i(X_i)$ (never $\gamma_i$ or $Q_{i,j}(X_j)$ for $j < i$). Therefore, the coefficient of $M_*$ in $f(R_1(X), \dots, R_r(X))$ equals its nonzero coefficient in $\mathrm{Tr}(\prod_{i \in I_0} Q_i(X_i)^{2^{w_i}-1})$. Since $M_*$ has positive degree, $f$ is nonconstant on $A$. This proves Theorem II.1. $\blacksquare$


# Part III. Elementary prime and field support

## 1. Bertrand's interval $(x, 2x]$ for every integer $x \ge 512$

Let $P(y) = \prod_{p \le y} p$ denote the product of all primes $p \le y$. We first prove by strong induction on integers $N \ge 1$ that $P(N) \le 4^N$ (and hence $P(y) = P(\lfloor y\rfloor) \le 4^y$ for all real $y \ge 1$). The cases $N = 1, 2$ hold ($1 \le 4$, $2 \le 16$). For even $N = 2n \ge 4$, $2n$ is composite so $P(2n) = P(2n - 1) \le 4^{2n-1} < 4^{2n}$. For odd $N = 2n + 1 \ge 3$, every prime $p \in (n + 1, 2n + 1]$ divides $\binom{2n+1}{n} = \frac{(2n+1)!}{n!(n+1)!}$, and since $\binom{2n+1}{n} = \binom{2n+1}{n+1}$ are two equal terms in $\sum_{k=0}^{2n+1} \binom{2n+1}{k} = 2^{2n+1}$, we have $\binom{2n+1}{n} \le 4^n$, giving $P(2n + 1) \le P(n + 1)\binom{2n+1}{n} \le 4^{n+1} 4^n = 4^{2n+1}$.

Now let $x \ge 512$ be an integer and suppose for contradiction that there is no prime in $(x, 2x]$. Consider $B = \binom{2x}{x}$. Since $B$ is the largest of the $2x + 1$ binomial coefficients summing to $4^x$, $B \ge 4^x / (2x + 1)$. By Legendre's formula, $v_p(B) = \sum_{j \ge 1} (\lfloor 2x / p^j \rfloor - 2\lfloor x / p^j \rfloor)$, in which each summand is in $\{0, 1\}$ and vanishes for $p^j > 2x$, so $p^{v_p(B)} \le 2x$ for every prime $p$. Primes $p \le \sqrt{2x}$ contribute at most $(2x)^{\sqrt{2x}}$ to $B$. For $p > \sqrt{2x}$, $v_p(B) \le 1$; for $\frac{2}{3}x < p \le x$, we have $2 \le \frac{2x}{p} < 3$ and $1 \le \frac{x}{p} < \frac{3}{2}$, so $\lfloor 2x/p\rfloor - 2\lfloor x/p\rfloor = 2 - 2 = 0$; and by assumption there is no prime in $(x, 2x]$. Thus all primes $p > \sqrt{2x}$ dividing $B$ satisfy $p \le \frac{2}{3}x$ and contribute at most $P(\frac{2}{3}x) \le 4^{2x/3}$. Hence
    4^x / (2x + 1) <= B <= (2x)^(sqrt(2x)) * 4^(2x/3),
and taking base-$2$ logarithms and dividing by $x$ yields
    2/3 <= sqrt(2/x) * log_2(2x) + log_2(2x + 1) / x.
For $x \ge 512$, both terms on the right-hand side are strictly decreasing functions of $x$: since $\ln 2 > 3/5$, $\ln(2x) \ge 10\ln 2 > 6 > 2$, so $\frac{d}{dx}\frac{\log_2(2x)}{\sqrt{x}} = \frac{2 - \ln(2x)}{(2\ln 2)x^{3/2}} < 0$; similarly $\frac{d}{dx}\frac{\log_2(2x+1)}{x} = \frac{2x/(2x+1) - \ln(2x+1)}{(\ln 2)x^2} < 0$ since $\frac{2x}{2x+1} < 1 < \ln(2x+1)$. At $x = 512 = 2^9$, since $1025 < 2048 = 2^{11}$, the right-hand side equals $\frac{1}{16}\cdot 10 + \frac{\log_2(1025)}{512} < \frac{320 + 11}{512} = \frac{331}{512} < \frac{2}{3}$, a contradiction. Thus for every integer $x \ge 512$ there is a prime $\ell$ with $x < \ell \le 2x$.

## 2. Existence of irreducible polynomials and exactness of the prime-degree test

In a splitting field $K_0$ of $H(X) = X^{2^\ell} - X$ over $\mathbb{F}_2$, the derivative is $H'(X) = 1$, so $H$ is squarefree and its $2^\ell$ roots form a subfield $K_0$ of size $2^\ell$. The Frobenius automorphism $\sigma(a) = a^2$ satisfies $\sigma^\ell = \mathrm{id}_{K_0}$, so for prime $\ell$ every $\sigma$-orbit on $K_0$ has size $1$ or $\ell$. Only $0$ and $1$ have orbit size $1$, and $2^\ell > 2$, so at least one orbit $O$ has size $\ell$. For any orbit $O$ of size $\ell$, $P_O(X) = \prod_{\alpha \in O} (X - \alpha)$ has $\sigma$-invariant coefficients (hence in $\mathbb{F}_2$) and is irreducible over $\mathbb{F}_2$ (any nonconstant factor over $\mathbb{F}_2$ has a $\sigma$-invariant root set $\emptyset \subsetneq O' \subseteq O$, forcing $O' = O$). Conversely, every irreducible factor of $H(X)$ over $\mathbb{F}_2$ is the orbit polynomial of one $\sigma$-orbit, hence has degree $1$ (namely $X$ or $X + 1$) or degree $\ell$.

Now let $P(X) \in \mathbb{F}_2[X]$ be monic of prime degree $\ell$. If $X^{2^\ell} \equiv X \pmod{P}$ and $\gcd(P, X^2 - X) = 1$, then $P$ divides the squarefree polynomial $X^{2^\ell} - X$ and is coprime to $X(X+1)$, so $P$ is a product of distinct irreducible factors of degree $\ell$; since $\deg(P) = \ell$, $P$ is itself irreducible. Conversely, if $P$ is irreducible of degree $\ell$, then $\mathbb{F}_2[X]/(P)$ is a field of size $2^\ell$, so by Lagrange's theorem on its multiplicative group of order $2^\ell - 1$ every residue satisfies $z^{2^\ell} = z$, hence $X^{2^\ell} \equiv X \pmod{P}$, and since $\ell \ge 2$ is irreducible, $\gcd(P, X^2 - X) = 1$. Thus the evaluator's test in Step 6 is exact in both directions, and runs in $O(\ell^3)$ bit operations.


# Part IV. Complete actual-length analysis and error constants 51/4 and 25/2

## 1. Parameter validity and polynomial-time uniformity at every integer $N \ge 2^{20}$

Fix any integer $N \ge 2^{20}$ and put $u = N^{1/5} \ge 16$, so $N = u^5$ and $N^{3/5} = u^3 \ge 4096$. Step 2 of the evaluator finds the least positive integer $b$ with $25^5 b^5 \ge 24^5 N^3$, which by binary search takes $O(\log N)$ exact integer comparisons and satisfies
    (24/25) * u^3 <= b < (24/25) * u^3 + 1.
In particular, $b \ge \lceil (24/25) \cdot 4096 \rceil = 3933 \ge 512$, so Part III Section 1 guarantees a prime $\ell$ with $b < \ell \le 2b$ (odd since $\ell > 3933$). Writing $A = \ell / u^3$, we have
    24/25 < A < 48/25 + 2 / u^3.                                (A-bounds)
Since $u \ge 16$, $2/u^3 \le 1/2048$, so $24/25 < A < 98329/51200 < 2$.

In Step 4, $m = N - \ell$ and $r = \lceil m / \ell \rceil = \lceil (N - \ell)/\ell \rceil = \lceil N / \ell \rceil - 1$. Observe that for **every** integer $N$ and $\ell \ge 1$,
    N / ell - 1 <= r < N / ell.
Indeed, write $N = q\ell + \rho$ with $0 \le \rho < \ell$: if $\rho = 0$, then $(N - \ell)/\ell = q - 1$ is an integer so $r = q - 1 = N/\ell - 1 < N/\ell$; if $1 \le \rho < \ell$, then $q - 1 < (N - \ell)/\ell < q$ so $r = q = (N - \rho)/\ell < N/\ell$. Consequently,
    u^2 / A - 1 <= r < u^2 / A <= (25/24) * u^2.                (r-bounds)
For $u \ge 16$ and $A < 2$, $r \ge \lceil 256 / 2 \rceil - 1 = 127$, so $2 \le t = \lceil\sqrt{r}\rceil < r$, and $r + 1 < (25/24)u^2 + 1 < (24/25)u^3 < \ell$. Thus all hypotheses of Theorem II.1 hold for every integer $N \ge 2^{20}$. Trial division for $\ell$, modular squarings and Euclid for $P$, $O(r^2)$ field multiplications for $y_i$, and $O(rt)$ operations for $\mathrm{Sym}_t$ all run in deterministic time $O((N+1)^4)$.

## 2. Hard prefix restriction and exact remainder bound $E(N)$

Fix any integer $N \ge 2^{20}$. By Part III Section 2, choose one monic irreducible polynomial $P(X) \in \mathbb{F}_2[X]$ of degree $\ell$, and fix the first $\ell$ input bits of an arbitrary full-$B_2$ circuit $C$ for $F_N$ to the coefficients of $P$. This restriction fixes $\ell$ input wires to constants without adding any gates, producing a full-$B_2$ circuit of size $\le |C|$ for the restricted payload function $G_P : \mathbb{F}_2^m \to \mathbb{F}_2$ on $m = N - \ell$ actual input bits.

In Step 7, appending $r\ell - m$ zeros ($0 \le r\ell - m < \ell$) to the $m$-bit payload defines an injective $\mathbb{F}_2$-linear embedding $\iota : \mathbb{F}_2^m \hookrightarrow K^r$. Every affine subspace of $\mathbb{F}_2^m$ of dimension $D$ maps under $\iota$ to an affine subspace of $K^r$ of the same dimension $D$. Put
    T = r*(r + 3)/2 + ell * sqrt(r),        d_N = floor(T).
By Theorem II.1, $G_P$ is nonconstant on every affine subspace of $\mathbb{F}_2^m$ of dimension strictly greater than $d_N$. Applying Theorem I.1 to $G_P$ therefore yields
    C_B2(F_N) >= C_B2(G_P) >= 4*m - 8*(d_N + 1)
               >= 4*(N - ell) - 8*(r*(r + 3)/2 + ell*sqrt(r) + 1)
                = 4*N - [ 4*ell + 4*r^2 + 12*r + 8*ell*sqrt(r) + 8 ].    (1)

Notice that because the $m = N - \ell$ payload bits are injectively zero-padded into $r$ blocks rather than discarding tail bits, the linear header cost in (1) is $4\ell$ (with $r < N/\ell$), saving $4\ell$ gates.

Now substitute $\ell = A u^3$ and $r < u^2 / A$ from `(r-bounds)` into the leading $u^4 = N^{4/5}$ terms of (1):
    4*r^2 + 8*ell*sqrt(r) < u^4 * g(A),        where  g(A) = 4 / A^2 + 8 * sqrt(A).
On $(0, 2]$, $g''(A) = 24 A^{-4} - 2 A^{-3/2} = 2 A^{-4}(12 - A^{5/2}) > 0$ (as $A^{5/2} \le 2^{5/2} = \sqrt{32} < 6 < 12$), so $g$ is strictly convex on $[24/25, 48/25 + 2/u^3]$ and attains its maximum on that interval at one of the two endpoints:
* **Left endpoint $A_0 = 24/25$:** We have $g(24/25) = 4(25/24)^2 + 8\sqrt{24/25} = \frac{625}{144} + \frac{16}{5}\sqrt{6}$. Since $(49/20)^2 = 2401/400 > 6$, we have $\sqrt{6} < 49/20$, and hence
      g(24/25) < 625/144 + (16/5)*(49/20) = 625/144 + 196/25 = 43849/3600.
* **Right endpoint $A_1 = 48/25 + 2/u^3$:** On $[48/25, 2]$, since $A > 1$, $g'(A) = 4 A^{-1/2} - 8 A^{-3} < 4(1) - 0 = 4$, so by the mean value theorem $g(A_1) \le g(48/25) + 4(2/u^3) = \frac{625}{576} + \frac{32}{5}\sqrt{3} + \frac{8}{u^3}$. Since $(97/56)^2 = 9409/3136 > 9408/3136 = 3$, we have $\sqrt{3} < 97/56$, and for $u \ge 16$ ($8/u^3 \le 1/512$),
      g(A_1) < 625/576 + (32/5)*(97/56) + 1/512 = 1963219/161280 < 43849/3600.

Therefore $g(A) < 43849/3600$ for every $A \in [24/25, 48/25 + 2/u^3]$. For the lower-order terms in (1), `(A-bounds)` and `(r-bounds)` give $4\ell < 4(\frac{48}{25}u^3 + 2) = \frac{192}{25}u^3 + 8$ and $12r < 12(\frac{25}{24}u^2) = \frac{25}{2}u^2$. Substituting into (1) proves, for every integer $N \ge 2^{20}$,
    C_B2(F_N) >= 4*N - E(N),
    E(N) = (43849/3600)*N^(4/5) + (192/25)*N^(3/5) + (25/2)*N^(2/5) + 16.    (2)

## 3. Explicit rational constants 51/4 (at $N \ge 2^{20}$) and 25/2 (at $N \ge 2^{25}$)

Dividing (2) by $u^4 = N^{4/5}$, each negative power of $u$ is decreasing in $u$:
* For every integer $N \ge 2^{20}$ ($u \ge 16$):
      E(N) / N^(4/5) <= 43849/3600 + (192/25)/16 + (25/2)/256 + 16/65536
                      = 43849/3600 + 12/25 + 25/512 + 1/4096
                      = 11712937 / 921600 < 11750400 / 921600 = 51/4,
  which proves $C_{B_2}(F_N) \ge 4N - \frac{51}{4}N^{4/5}$ for every integer $N \ge 2^{20}$.
* For every integer $N \ge 2^{25}$ ($u \ge 32$):
      E(N) / N^(4/5) <= 43849/3600 + (192/25)/32 + (25/2)/1024 + 16/1048576
                      = 43849/3600 + 6/25 + 25/2048 + 1/65536
                      = 183324673 / 14745600 < 184320000 / 14745600 = 25/2,
  which proves $C_{B_2}(F_N) \ge 4N - \frac{25}{2}N^{4/5}$ for every integer $N \ge 2^{25}$.

For any real $\varepsilon > 0$, setting $N_\varepsilon = \max(2^{20}, \lceil (51 / (4\varepsilon))^5 \rceil)$ ensures $\frac{51}{4} N^{4/5} \le \varepsilon N$ and hence $C_{B_2}(F_N) \ge (4 - \varepsilon)N$ for every integer $N \ge N_\varepsilon$. Thus $\liminf_{N \to \infty} C_{B_2}(F_N)/N \ge 4$. Moreover, by (2), $\limsup_{N\to\infty} \max\{0, 4N - C_{B_2}(F_N)\} / N^{4/5} \le 43849/3600$.


# Part V. Unscored structural theorems: diagonal shared-partner elimination and arbitrary-depth admissible-rank subspace extraction

## 1. Unscored Theorem A: Diagonal shared-partner elimination in full-cube fair representations

**Unscored Theorem A.** *Let $d \ge 0$ be an integer, and let a full-cube fair representation on $n$ free parameters $u \in \mathbb{F}_2^n$ (with a square invertible auxiliary matrix $A v = N u + c$, rows of at most $3$ nonconstant terms, an acyclic full-$B_2$ gate list, and total size $S = \text{rows} + \text{gates}$) compute a function $f : \mathbb{F}_2^n \to \mathbb{F}_2$ that is nonconstant on every affine subspace of dimension $> d$. Suppose the representation contains $t$ minimal AND-type gates $G_j = (U_j \oplus \alpha_j)(W_j \oplus \beta_j) \oplus \gamma_j$ ($1 \le j \le t$) such that for each $j$, $U_j$ and $W_j$ are distinct auxiliary variables of degree $2$ appearing only in $G_j$ and in $3$-term equations $E_{j,1} : U_j \oplus s_j \oplus b_{j,1} = c_{j,1}$ and $E_{j,2} : W_j \oplus s_j \oplus b_{j,2} = c_{j,2}$ respectively, where $\{s_j, b_{j,1}, b_{j,2}\}_{j=1}^t$ are $3t$ distinct free parameters and each $s_j$ has degree $2$ (appearing only in $E_{j,1}$ and $E_{j,2}$). If $n - 2t \ge d + 2$, then*

    S >= 4*n + 3*t - 8*(d + 1),

*and if in addition $G_1, \dots, G_t$ have at least $t$ distinct gate children outside $\{G_1, \dots, G_t\}$, then*

    S >= 4*n + 4*t - 8*(d + 1).

*Proof.* Since $\{s_j, b_{j,1}, b_{j,2}\}_{j=1}^t$ are $3t$ distinct free parameters, we may simultaneously impose the $t$ disjoint free-parameter identifications (`UF`)
    b_(j,1) <- b_(j,2) + (c_(j,1) + c_(j,2) + alpha_j + beta_j + 1)        (1 <= j <= t).
This restricts the parameter space to an affine subspace of dimension $a' = M' = n - t \ge d + t + 2 \ge d + 2$ while keeping $q' = 0$. In each pair of equations $E_{j,1}, E_{j,2}$, we have $(U_j \oplus \alpha_j) \oplus (W_j \oplus \beta_j) = b_{j,1} \oplus b_{j,2} \oplus c_{j,1} \oplus c_{j,2} \oplus \alpha_j \oplus \beta_j = 1$, so $(U_j \oplus \alpha_j) = (W_j \oplus \beta_j) \oplus 1$ identically on the restricted cube. Hence every minimal AND gate $G_j = (U_j \oplus \alpha_j)(W_j \oplus \beta_j) \oplus \gamma_j \equiv \gamma_j$ computes the constant $\gamma_j$.

Replace every $G_j$ ($1 \le j \le t$) by the constant $\gamma_j$ and delete $G_1, \dots, G_t$ by `D1`. Because $U_j$ and $W_j$ had degree $2$ with children $G_j$ and $E_{j,1}$ (respectively $G_j$ and $E_{j,2}$), removing $G_1, \dots, G_t$ reduces the degree of every auxiliary variable $U_j$ and $W_j$ to $1$, with sole remaining occurrence in $E_{j,1}$ and $E_{j,2}$ respectively (since $a' \ge d + 2$, no auxiliary variable is the designated output). Furthermore, the equations $E_{j,1}, E_{j,2}$ ($1 \le j \le t$) are $2t$ distinct equations whose only auxiliary variables are $U_j$ and $W_j$ respectively. Applying `D2'` (isolated auxiliary variable in an equation) to $(U_j, E_{j,1})$ and $(W_j, E_{j,2})$ for each $j = 1, \dots, t$ deletes all $2t$ equations $E_{j,1}, E_{j,2}$ and their auxiliary variables $U_j, W_j$ without substitution.

At this point:
1. All $3t$ distinct objects $\{G_j, E_{j,1}, E_{j,2}\}_{j=1}^t$ have been deleted, so $s' \le S - 3t$.
2. For each $j = 1, \dots, t$, the free parameter $s_j$ appeared only in $E_{j,1}$ and $E_{j,2}$, both of which are now deleted; hence $\deg'(s_j) = 0$. Since the $t$ parameters $b_{1,1}, \dots, b_{t,1}$ were eliminated by the identifications and the $t$ surviving free parameters $s_1, \dots, s_t$ have degree $0$, at most $n - 2t$ free parameters have positive degree ($I' \le n - 2t$). Any further `D1`/`D3` simplifications caused by substituting $b_{j,1} \leftarrow b_{j,2} \oplus \text{const}$ into other equations can only decrease $s'$ and cannot increase $I'$.
Thus at the resulting valid state $(M', s', I', q') = (n - t, s' \le S - 3t, I' \le n - 2t, 0)$, the potential $V' = s' + 4I' + 2q'$ satisfies
    V' <= (S - 3*t) + 4*(n - 2*t) = S + 4*n - 11*t.
Running the elimination loop of Theorem I.1 from this state (starting at $a' = n - t$ and terminating at $a_T \le d + 1$ with $\Delta V \le 8\,\Delta a$) gives
    S + 4*n - 11*t >= V' >= 8*(a' - d - 1) = 8*(n - t - d - 1),
which rearranges to $S \ge 4n + 3t - 8(d + 1)$.

If in addition $G_1, \dots, G_t$ have at least $t$ distinct gate children $D_1, \dots, D_t \notin \{G_1, \dots, G_t\}$, each $D_j$ acquires a constant pin when $G_1, \dots, G_t$ are replaced by constants and is deleted by `D1`, giving $4t$ distinct deleted objects $\{G_j, D_j, E_{j,1}, E_{j,2}\}_{j=1}^t$. This improves $s' \le S - 4t$ and $V' \le S + 4n - 12t$, yielding $S \ge 4n + 4t - 8(d + 1)$. $\blacksquare$

## 2. Unscored Theorem B: Arbitrary-depth admissible-rank and topological linearization subspace bound

**Unscored Theorem B.** *Let $C$ be an acyclic XOR/AND circuit on $n$ input bits $x \in \mathbb{F}_2^n$ with AND-gate set $\mathcal{A} = \mathcal{G}_1 \cup \mathcal{G}_2$, where every gate is either an XOR-type gate $u \oplus v \oplus c$ or an AND-type gate $(u \oplus \alpha)(v \oplus \beta) \oplus \gamma$, $\mathcal{G}_1$ consists of all level-$1$ AND gates (those whose two input pins have no AND-gate ancestors, so both input pins compute affine forms $u_g, v_g$ in $x$), and $\mathcal{G}_2 = \mathcal{A} \setminus \mathcal{G}_1$. Let $\mathcal{G}_{2,\mathrm{top}} \subseteq \mathcal{G}_2$ be any subset of $\mathcal{G}_2$ with no directed path to any other AND gate in $\mathcal{A}$. Let $X_1 \subseteq \{x_1, \dots, x_n\}$ be the set of input variables appearing with nonzero linear coefficient in at least one input form of a gate in $\mathcal{G}_1$, and let $S_1 \subseteq X_1$ be any subset such that no gate in $\mathcal{G}_1$ has both of its affine input forms mentioning a variable in $S_1$. Then on an affine subspace $V_m \subseteq \mathbb{F}_2^n$ of dimension at least $|S_1| + (n - |X_1|) - |\mathcal{G}_2 \setminus \mathcal{G}_{2,\mathrm{top}}|$, every gate in $\mathcal{G}_1 \cup (\mathcal{G}_2 \setminus \mathcal{G}_{2,\mathrm{top}})$ restricts to an affine function, and if $r_{\mathrm{top}} \le 2|\mathcal{G}_{2,\mathrm{top}}|$ denotes the dimension of the span of the homogeneous linear parts (on $V_m$) of the input forms of $\mathcal{G}_{2,\mathrm{top}}$, then $C$ is constant on an affine subspace $V \subseteq V_m$ of dimension at least*

    dim(V) >= |S_1| + (n - |X_1|) - |G_2 \ G_(2,top)| - min(|G_(2,top)|, floor(r_top / 2)) - 1.

*Proof.* First fix all input variables in $X_1 \setminus S_1$ to $0$. Because $S_1 \subseteq X_1$, this imposes $|X_1 \setminus S_1| = |X_1| - |S_1|$ independent coordinate constraints and defines an affine subspace $V_0 \subseteq \mathbb{F}_2^n$ of dimension
    dim(V_0) = n - (|X_1| - |S_1|) = |S_1| + (n - |X_1|).
For every level-$1$ AND gate $g = (u_g \oplus \alpha_g)(v_g \oplus \beta_g) \oplus \gamma_g \in \mathcal{G}_1$, both $u_g$ and $v_g$ depend only on variables in $X_1$, and by hypothesis on $S_1$ at least one of $u_g, v_g$ has zero linear coefficient on every variable in $S_1$, hence depends only on $X_1 \setminus S_1$ and is constant on $V_0$. Since one factor of $g$ is constant on $V_0$, every $g \in \mathcal{G}_1$ restricts on $V_0$ to an affine function.

Next, order $\mathcal{G}_2 \setminus \mathcal{G}_{2,\mathrm{top}} = \{h_1, \dots, h_m\}$ ($m = |\mathcal{G}_2 \setminus \mathcal{G}_{2,\mathrm{top}}|$) topologically and construct affine subspaces $V_0 \supseteq V_1 \supseteq \dots \supseteq V_m$ with $\dim(V_k) \ge \dim(V_{k-1}) - 1$ such that all gates in $\mathcal{G}_1 \cup \{h_1, \dots, h_k\}$ restrict on $V_k$ to affine functions. At step $k$, every AND-gate ancestor of the two inputs $p_k, q_k$ of $h_k = (p_k \oplus \alpha)(q_k \oplus \beta) \oplus \gamma$ lies in $\mathcal{G}_1 \cup \{h_1, \dots, h_{k-1}\}$ (as no gate in $\mathcal{G}_{2,\mathrm{top}}$ has a path to another AND gate). By the induction hypothesis, all such ancestors are affine on $V_{k-1}$, and since only XOR-type gates intervene, both $p_k \oplus \alpha$ and $q_k \oplus \beta$ are affine on $V_{k-1}$. If $p_k \oplus \alpha$ is nonconstant on $V_{k-1}$, set $V_k = \{x \in V_{k-1} : p_k(x) \oplus \alpha = 0\}$ ($\dim(V_k) = \dim(V_{k-1}) - 1$), on which $h_k \equiv \gamma$; otherwise $p_k \oplus \alpha$ is already constant on $V_{k-1}$, so set $V_k = V_{k-1}$, on which $h_k$ is affine. After $m$ steps, $\dim(V_m) \ge |S_1| + (n - |X_1|) - |\mathcal{G}_2 \setminus \mathcal{G}_{2,\mathrm{top}}|$, and every gate in $\mathcal{G}_1 \cup (\mathcal{G}_2 \setminus \mathcal{G}_{2,\mathrm{top}})$ is affine on $V_m$.

Finally, on $V_m$, every gate $g \in \mathcal{G}_{2,\mathrm{top}}$ has only ancestors in $\mathcal{G}_1 \cup (\mathcal{G}_2 \setminus \mathcal{G}_{2,\mathrm{top}})$ and XOR gates, so its two shifted input pins compute affine functions $L_{g,1}(y)$ and $L_{g,2}(y)$ in coordinates $y \in \mathbb{F}_2^{\dim(V_m)}$ parameterizing $V_m$, and the circuit output on $V_m$ has the form
    out(y) = ell_0(y) + sum_(g in T_out) L_(g,1)(y) * L_(g,2)(y)
for an affine function $\ell_0(y)$ and a subset $T_{\mathrm{out}} \subseteq \mathcal{G}_{2,\mathrm{top}}$. We bound the codimension in $V_m$ needed to linearize the sum of products in two ways:
1. **Gate-by-gate bound ($\le |\mathcal{G}_{2,\mathrm{top}}|$ constraints):** For each $g \in \mathcal{G}_{2,\mathrm{top}}$ in turn, if $L_{g,1}$ is nonconstant on the current subspace, restrict to the hyperplane $L_{g,1} = 0$ (1 constraint); otherwise $L_{g,1} L_{g,2}$ is already affine. This linearizes $\mathrm{out}$ using at most $|\mathcal{G}_{2,\mathrm{top}}|$ constraints.
2. **Span-dimension bound ($\le \lfloor r_{\mathrm{top}} / 2 \rfloor$ constraints):** Let $z_1, \dots, z_r$ ($r = r_{\mathrm{top}}$) be a basis for the span of the homogeneous linear parts of $\{L_{g,1}, L_{g,2} : g \in \mathcal{G}_{2,\mathrm{top}}\}$ on $V_m$, so $\mathrm{out}(y) = \sum_{1 \le i < j \le r} c_{i,j} z_i z_j \oplus \ell_1(y)$ for an affine function $\ell_1(y)$. We prove by strong induction on $r \ge 0$ that any quadratic sum $\sum_{1 \le i < j \le r} c_{i,j} z_i z_j$ becomes affine after imposing at most $\lfloor r/2 \rfloor$ linear constraints on $z_1, \dots, z_r$. For $r \le 1$ or when all $c_{i,j} = 0$, $0$ constraints suffice. Otherwise, reorder $z_1, \dots, z_r$ so $c_{1,2} = 1$, and define the invertible linear change of variables $z'_1 = z_1 \oplus \sum_{j=3}^r c_{2,j} z_j$ and $z'_2 = z_2 \oplus \sum_{j=3}^r c_{1,j} z_j$. Then $\sum_{1 \le i < j \le r} c_{i,j} z_i z_j = z'_1 z'_2 \oplus Q'(z_3, \dots, z_r)$ for a quadratic form $Q'$ depending only on the $r - 2$ variables $z_3, \dots, z_r$. Imposing the single linear constraint $z'_1 = 0$ eliminates $z'_1 z'_2$, and by the induction hypothesis $Q'(z_3, \dots, z_r)$ becomes affine after at most $\lfloor (r - 2)/2 \rfloor$ additional constraints, giving $1 + \lfloor (r - 2)/2 \rfloor = \lfloor r/2 \rfloor$ constraints in total.

Taking the better of the two bounds yields an affine subspace $V' \subseteq V_m$ of codimension at most $\min(|\mathcal{G}_{2,\mathrm{top}}|, \lfloor r_{\mathrm{top}}/2 \rfloor)$ on which $\mathrm{out}$ is an affine function $\ell_*(y)$. If $\ell_*$ is nonconstant on $V'$, restricting to the hyperplane $\ell_*(y) = 0$ imposes $1$ final constraint and makes $C$ constant on an affine subspace $V \subseteq V_m$ of the stated dimension. $\blacksquare$
