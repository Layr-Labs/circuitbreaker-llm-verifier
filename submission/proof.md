# Coefficient four with a smaller explicit loss

For the total language below, we prove that, for every sufficiently large integer input length $n$,
$$
 C_{B_2}(f_n)\ \ge\ 4n-\frac{359}{25}n^{4/5}.
$$
Here $n$ includes the field-description bits. The full circuit can choose any of the sixteen two-pin truth tables at every gate. The proof first finds a hard restriction at each length, then proves a gate bound that applies to those restrictions.

The algebraic and circuit arguments below are adapted from the complete supplied proof of public Yukon submission `c733ccd7-71d9-4555-b156-100d49c75f54`, evaluated at commit `e1f9a7d30f29f3a79e17db1b31701ef7565fae2a`. Its original `proof.md` has SHA-256 `c7a5d52a3b04036c625debd2680160e8b1e63a70888e1d4705a5e4ca590d413c`. We reproduce the needed arguments here; no prior judgment is a premise. Our change is the prime cutoff and its all-length analysis, reducing the stated loss from $61/4$ to $359/25$. This does not increase the leading coefficient four.

## The function and its verifier

Set $f_0=0$. For $n\ge1$ let $s$ be the least positive integer with $s^5\ge4n^3$, and let $k$ be the least prime strictly above $s$. If $k\ge n$, return zero. Otherwise the first $N=n-k$ bits are data, and the last $k$ bits give the coefficients $h_0,\ldots,h_{k-1}$ of the monic polynomial $h(Z)=Z^k+\sum_{a=0}^{k-1}h_a Z^a$ over $\mathbb F_2$. Set $r=\lceil N/k\rceil$, $t=\lceil\sqrt r\rceil$ and $u_j=r+1-j$ for $1\le j\le r$. Pad the data at the end with zeros and divide the resulting $rk$ bits into $r$ blocks. In each block read coefficients from low to high degree, obtaining $X_j$ in the ring $R=\mathbb F_2[Z]/(h)$. Compute
$$
 f_n(x)=\pi\!\left(\sum_{I\subseteq\{1,\ldots,r\},\ |I|=t}
                         \ \prod_{j\in I} X_j^{2^{u_j}-1}\right),
$$
where $\pi$ takes the constant coefficient of the unique remainder of degree less than $k$. In particular this defines the answer on every reducible-modulus input.

One deterministic algorithm evaluates every length. It finds $s$ by integer fifth-power comparisons, then tests successive integers for primality by trial division. The prime estimate proved below gives $k=O(n^{3/5})$; the finitely many smaller lengths do not affect this bound or polynomial running time. Polynomial multiplication and reduction modulo a monic polynomial take polynomial bit time. Each block exponent has $u_j$ binary digits. Repeated squaring and multiplication uses $O(u_j)$ ring operations. Write the resulting powers as $b_j$. Starting with $E_0=1$ and $E_1=\cdots=E_t=0$, on block $j$ update $E_i\gets E_i+b_j E_{i-1}$ for $i=\min(t,j),\ldots,1$ in descending order. Splitting subsets by membership of $j$ shows inductively that $E_i$ is the degree-$i$ elementary symmetric sum on the processed blocks. Return $\pi(E_t)$. All loops and bit lengths are polynomial in $n$. Define $L=\{x:f_{|x|}(x)=1\}$. Take witness bound $W(n)=0$ and let $V(x,w)$ accept exactly if $w$ is empty and this algorithm returns one. A yes input has this witness and any accepting witness implies a yes input. Only the actual input $x$, not a witness, is supplied to the circuits.

We use the classical prime number theorem (Hadamard and de la Vallée Poussin, 1896) in the standard form
$$
 \Pi(x)\sim x/\log x\qquad(x\longrightarrow\infty),
$$
where $\Pi(x)$ counts primes at most the real number $x$ and this logarithm is natural. A precise reference is NIST DLMF, version 1.2.8 (2026-09-15), Section 27.12, equation 27.12.4: the first term of its asymptotic expansion is $\Pi(x)\sim x/\log x$ ([reference](https://dlmf.nist.gov/27.12#E4)). This familiar standard theorem is the only analytic number-theory input; it is not re-proved here. Its precise consequence needed in this proof follows directly. For every fixed $\eta>0$,
$$
 \Pi((1+\eta)x)-\Pi(x)
       =(\eta+o(1))\,x/\log x>0
$$
for every sufficiently large real $x$. Thus the least prime $p(x)>x$ obeys $1<p(x)/x\le1+\eta$ eventually. As $\eta$ is arbitrary, $p(x)/x\to1$. In our fixed construction $s/n^{3/5}\to2^{2/5}$, so
$$
 z_n:=k/n^{3/5}\longrightarrow 2^{2/5}. \tag{P}
$$
All of these limits hold over every sufficiently large integer input length, not just a chosen sequence of primes or lengths. In particular the scan terminates and has polynomial cost, and $k<n$ eventually.

## A restriction that varies on every large affine space

At each sufficiently large length fix the last $k$ input bits to one irreducible $h$ of degree $k$. Such a polynomial exists. In a splitting field the $2^k$ distinct roots of $Z^{2^k}-Z$ form a field by the Frobenius identities and closure under inversion of nonzero roots. The roots are simple since the formal derivative is one. This field has degree $k$ over $\mathbb F_2$. Since $k$ is prime, any element outside $\mathbb F_2$ has degree $k$ by the tower law. Its minimal polynomial works. This existential choice is a restriction on circuits for an already fixed function. No part of its evaluator takes advice.

Write $K=\mathbb F_2[Z]/(h)$ on the restricted input and take all dimensions below over $\mathbb F_2$. We prove that the restricted function takes both truth values on any affine subspace of the $N$ data bits of dimension
$$
 D>kt+\frac{r(r+1)}2. \tag{1}
$$

**Products of subspaces.** For nonzero subspaces $A,B\subseteq K$, let $P$ be the span of all $ab$ with $a\in A$, $b\in B$. We claim $\dim P\ge\min(k,\dim A+\dim B-1)$. Induct on $\dim B$. The one-dimensional case follows by multiplication by a nonzero scalar. The case $A=K$ follows as well. Otherwise, when $\dim B>1$, there is an $e\in K^*$ for which $0<\dim(A\cap eB)<\dim B$. Indeed, in the absence of such an $e$, fix $b\in B$ nonzero. For each nonzero $a\in A$, the space $(a/b)B$ meets $A$ at $a$ and so would lie inside $A$. Choose $z\in B/b$ outside $\mathbb F_2$. Then $A$ would be stable under $z$ and hence under the whole field $\mathbb F_2[z]=K$, contradicting $0\ne A\ne K$. For the indicated $e$, put $A'=A+eB$ and $B'=B\cap e^{-1}A$. Their dimensions sum to $\dim A+\dim B$, and $0<\dim B'<\dim B$. Products from $A'B'$ lie in $P$, since $B'\subseteq B$ and $eB'\subseteq A$. Induction proves the claim. Iterating it shows that products of $t$ nonzero subspaces $H_i$ span dimension at least $\min(k,1+\sum_i(\dim H_i-1))$.

**One block.** If $U\subseteq K$ has dimension $a$ and $1\le u\le a$, the values of
$$
 M(v_1,\ldots,v_u)=\det\big(v_\ell^{2^h}\big)_{0\le h<u,\,1\le\ell\le u},\qquad v_\ell\in U,
$$
span dimension at least $a-u+1$. Fix $u-1$ independent entries and vary the last one. The resulting map is binary-linear. Dependence of the entries makes the determinant zero. For independent entries, a zero determinant would supply a nonzero polynomial $\sum_{h=0}^{u-1} c_h T^{2^h}$ vanishing at the entries. Additivity would give $2^u$ roots, more than its ordinary degree. Hence the kernel is exactly the span of the fixed entries, proving the rank bound, including $u=1$.

**A coefficient on the given affine space.** The padded block encoding injects the affine space into $K^r$. In its direction space let $V_j$ be the subspace of vectors zero in the first $j-1$ blocks and let $U_j$ be the image of $V_j$ in block $j$. The kernel there is $V_{j+1}$. Choosing linear lifts of bases of the successive $U_j$ gives a basis of the whole direction space. With $a_j=\dim U_j$ one has $\sum a_j=D$ and $a_j\le k$. An affine point in block $j$ depends only on groups $1,\ldots,j$, with its own group's contribution precisely $y_j\in U_j$.

Put $q_j=\max(0,a_j-u_j)$. By (1), $\sum q_j>kt$ since $\sum u_j=r(r+1)/2$. Each $q_j\le k$, so at least $t$ of them are positive. Take the $t$ largest, indexed by $I$. Their sum is at least $(t/r)\sum q_j>kt^2/r\ge k$. For $i\in I$, form $H_i$ from the determinants in the preceding paragraph with $U=U_i$ and $u=u_i$. The product-span inequality now gives dimension $k$ because $\sum_{i\in I}(\dim H_i-1)\ge\sum_{i\in I}(a_i-u_i)>k$. By multilinearity the products of actual determinants span $K$. The functional $\pi$ is nonzero ($\pi(1)=1$), so we may choose $u_i$ directions $v_{i,\ell}$ in each selected group whose product of determinants has projection 1. In particular the chosen directions in each individual group are independent.

In the lifted basis take $y_i=\sum_\ell z_{i,\ell}v_{i,\ell}$ for $i\in I$, with independent Boolean bits $z_{i,\ell}$, and put all other groups to zero. These points give a cube inside the given affine space, with its original affine offset retained. Give each bit from group $i$ weight $1/u_i$. The product of all these bits is a squarefree monomial of weight $t$. In a summand indexed by a $t$-subset $J$, block $j\in J$ contributes the $u_j$ factors $X_j^{2^h}$ for $0\le h<u_j$. Each factor is affine in bits from groups $i\le j$ by Frobenius linearity. Any term picked from that factor has weight at most $1/u_j$. Equality requires a bit of group $j$, since the sequence $u_j$ strictly decreases and a constant has weight zero. The entire term therefore has weight at most $t$ even before reducing $z^2=z$. Any repeated bit strictly lowers its weight on reduction. To produce the target monomial, equality must hold in each factor and all chosen bits must be distinct. This forces $J=I$. The sum over permutations of the bits in each block is precisely its displayed determinant, since signs disappear in characteristic two. Applying $\pi$ to the product gives coefficient 1 on the target monomial. A Boolean function on a cube has a unique squarefree polynomial, as follows by induction using its values at the last bit 0 and 1. This nonempty monomial proves nonconstancy and establishes (1).

The starting point for this construction is the block polynomial of Eli Ben-Sasson and Swastik Kopparty, *Affine Dispersers from Subspace Polynomials*, 12 March 2010 author version, Section 7 and Theorem 7.3. The product-span and weighted-monomial arguments above give the binary instance used here. [Author preprint](https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf)

## From monochromatic affine spaces to a gate bound

We prove the following finite statement. If a Boolean function on $n$ bits takes both values on every affine space of dimension $d$, where $1\le d\le n$, then every exact single-output full-binary circuit for it has at least $4n-10d$ gates. Our direct source, submission `c733ccd7-71d9-4555-b156-100d49c75f54`, attributes the matched-parameter idea and finite bound to Om Swostik Mishra and Rahul Santhanam, ECCC TR26-260, revision 1 (3 November 2026), Sections 2–4, Theorem 4.1 ([source's citation](https://eccc.weizmann.ac.il/report/2026/260/revision/1/download/)). That paper was not retrievable during this preparation. The complete argument we use is supplied and checked below, rather than assuming its theorem from that citation. We use charges on singletons and pairs to give the accounting.

In this section $n$ denotes the number of inputs to the circuit under consideration. All the linear algebra below is over the field with two elements. Start with any circuit $D$ with $g$ gates. We seek an affine space on which its output is constant. At an intermediate stage let $u$ be $m$ independent coordinates and let $T$ be an injective affine map from their space to the original $n$-bit input space. Partition these coordinates into singletons and disjoint pairs. A pair $(p,r)$ imposes the promise $pr=0$. Let $Z$ be the assignments obeying all promises and $b$ the number of blocks (singletons plus pairs). Selecting one member of each pair and setting it to zero cuts out a $b$-dimensional affine space inside $Z$.

We represent $D(T(u))$ on $Z$ by a system of sparse affine equations followed by ordinary acyclic binary gates. With $k$ equations there are $k$ auxiliary bits $v$, uniquely determined by $A v=B u+c$ with $A$ invertible. Each individual equation mentions at most three variables after cancellation. The gates may read coordinates, auxiliary bits, constants and earlier gates. Equations never read gates. The output designates a bit (possibly complemented) or a constant. The initial map is the identity, the coordinates are all singletons, there are no equations, and the gates are those of $D$. Write $s=k+\#\text{gates}$, so initially $s=g$. These descriptions after the initial step just record values on promised inputs. In particular we need no assumption that solutions of a linear system can be implemented by free gates in the original model.

Call a coordinate used if some equation mentions it or some gate reads it. Its children are the distinct equations and gates in this count. Two pins on the same gate give just one child. Output designation is tracked separately. Give a used singleton charge 4, an unused singleton charge 10, and a pair charge 6 regardless of its uses. Define

$$
             R=(\text{sum of block charges})-s.
$$

Initially $R\ge 4n-g$ and always $R\le10b$. Our reductions never decrease $R$. Thus $g<4n-10d$ would imply $b>d$ throughout. We may terminate immediately if the output is affine in the coordinates or is a product of two affine forms plus a constant. This includes a designated coordinate, auxiliary or constant, or a front product gate on affine inputs. Solving the auxiliary system is used to reason about its affine values, not to introduce free gates into the original circuit. Indeed, take the $b$-dimensional space inside $Z$ described above. An affine output is constant on that space or has a constant fiber of dimension $b-1$. For a product, if both factors are constant on the space we are done there, and otherwise a nonconstant factor has a zero fiber of dimension $b-1$. A complemented output only changes the constant. The injection $T$ carries the resulting space to an affine space of the same dimension in the original domain.

We record the elementary changes of representation. When deleting a determined value, replace **every** occurrence, including any designated output. Signs on wires can be folded into the consuming truth tables or the constants in affine equations. A terminal sign is recorded in the description, where it preserves monochromaticity. All original gates count toward $g$, including any unary output gate. There are two constant operations, four essentially unary operations, two parity operations, and eight operations of the form $(x+i)(y+j)+h$. Indeed, this last group consists of the indicators of the four individual input pairs and their complements. These account for all sixteen possible truth tables.

* Delete a dangling nondesignated gate, or bypass a gate which on the current domain is a constant or a signed copy of one of its pins. In particular this handles constants, unary gates, repeated pins and a constant pin. For a gate on two affine bits whose linear parts have rank at most one, their pair takes at most two values. If it takes two, either nonconstant pin distinguishes them, so any output is constant or a signed copy of such a pin. A deletion never introduces use of any previously unused coordinate because any nonconstant replacement was already a pin of the removed gate.
* An affine gate on coordinates and auxiliaries can be written as a new equation for a fresh auxiliary, with coefficient one on that auxiliary. This substitutes an equation for a gate, keeps the size $s$, and preserves unique solvability because the previous system determines the new bit's inputs. A row of length at most two can always be deleted. It contains an auxiliary, since a row of the invertible $A$ cannot be zero. Solve for such an auxiliary as a constant or a signed copy and substitute everywhere, discarding its row and column. This leaves exactly one extension of every parameter assignment, and substitution does not widen other rows. Each remaining square system is therefore invertible. Neither operation can turn an unused singleton into a used one.
* An auxiliary $v_j$ with solved expression $x+L(u_{\ne x})$ can exchange roles with an unpaired coordinate $x$. The other coordinates together with $v_j$ are affine bijective coordinates, so substitute the inverse transformation only in $T$. Keep all equations and gates as written, treating $x$ as auxiliary and $v_j$ as coordinate. For each new parameter tuple the old $x$ is uniquely determined by $v_j$ and the other coordinates, hence there is exactly one solution of the new square system. Both the old $x$ and new $v_j$ occur in equations. Otherwise $v_j$ could not depend on $x$, or the old auxiliary column for $v_j$ would be zero. This exchange preserves the charge, row lengths and all pairs.
* In particular if a singleton $x$ occurs only in one row $E$, we can delete that row at no loss of coordinate count. When the circuit output is not directly $x$ or its complement, select an auxiliary depending on $x$. Such an auxiliary exists because the column for $x$ in $B$ has a single one, so its product by $A^{-1}$ is nonzero. Make the exchange just described. The old $x$ is now auxiliary and occurs only in $E$. Drop both. Every solution to the remaining equations extends in exactly one way by solving $E$ for $x$, so the smaller system is invertible. The old input values are already recorded by the new map $T$. If the output directly names $x$ before such an operation, it is affine and the terminal argument applies at a completed state.

In each of these operations $R$ stays fixed or grows, since a size decrease helps and a singleton changing from used to unused gains charge. A change fixing a coordinate to a literal of another coordinate or to a constant composes $T$ with the corresponding affine injection, retaining agreement on the surviving promise set. Such restrictions never widen the rows or alter the invertible auxiliary matrix before simplification. When the chosen restriction of a pair implies its promise, we remove that promise.

Here are the reductions. At the start of a round check for a terminal output. Perform elementary deletions when applicable and restart the round after each deletion. If an affine gate with affine input wires is at the front of the acyclic part (its pins are constants, coordinates or auxiliaries), move it into the affine system and restart. These moves are bounded by the number of gates until the next proper reduction, since each removes a gate from the acyclic part. Always cancel equal terms in equations. We next make one of the following changes, returning to the start after the whole change. During each grouped change, elementary substitutions may be performed in any topologically valid order, and only the charges at its two ends are compared.

1. If a paired $p$ has no children, set its partner to zero and keep $p$ as an unused singleton. Its block charge goes from 6 to 10. If a paired $p$ has at least two distinct children, set $p=0$ and eliminate two of those distinct children. Every gate among them now has a constant pin and every equation among them now has length at most two. Eliminating one can only substitute a literal or constant into the other, which therefore remains deletable. The surviving singleton contributes at least 4, so the loss from 6 is paid for by the two deleted objects.

2. All paired coordinates now have exactly one child. Suppose two partners have the same child. If it is a gate, its truth table agrees at two of $00,01,10$. Restrict to the affine line through those two points, given by one of $p=0,r=0,p+r=1$, and bypass the gate as a constant. If the child is a row, it must be $p+r+v=c$ with an auxiliary $v$. Restrict to $p+r=1$ and eliminate the row by substituting the resulting constant $v$. One free coordinate parameterizes the line and does not occur elsewhere. The pair's charge 6 becomes an unused singleton's 10.

3. Suppose a paired $p$ occurs in a row $E$. Its partner has a distinct unique child $H$. Set that partner to zero and delete $H$, reducing either a gate or an equation. Because $H$ did not contain $p$, this substitution introduces no further appearances of $p$. Its remaining unique row $E$ can be eliminated by the singleton row operation above or by short-row substitution. Even when $H$ is itself the designated output, its replacement is a constant or a literal of an input other than $p$, so the singleton row operation does not discard the output. Thus the surviving singleton (or its exchanged coordinate) has charge at least 4 and two different objects were removed. This pays for the change from 6 to at least 4.

After these changes are excluded, each paired coordinate occurs in one gate, never a row, and partners occur at distinct gates. Take the first remaining gate $G$ in a topological ordering. If none remains, the output is affine. The first gate must have product form, since gates with affine operation at the front have been moved or deleted. Its two affine inputs have independent linear parts by the rank-one bypass. We can make both of them actual coordinates. Keep any input which is already a coordinate. For an auxiliary input, select an unpaired coordinate on which it has coefficient one, other than the input already kept if any. This is possible because no paired coordinate appears in any row, and independence excludes dependence solely on a kept coordinate. Exchange roles as above and repeat once if necessary. Affine bijections preserve linear independence. The trades leave all paired symbols and their occurrences unchanged. We now have $G=(x+i)(y+j)+h$ for distinct coordinates $x,y$. If its literal is the output, the terminal argument applies. Otherwise $G$ has at least one child gate (dangling gates were deleted).

4. For two singleton inputs, complement each coordinate if needed, absorbing the sign in all its uses and in $T$, so that $G=xy+h$. Impose the new pair $xy=0$. Replace $G$ everywhere by the constant $h$ and eliminate one of its children. The two used singletons had combined charge 8. Their pair has charge 6. The two distinct eliminations compensate exactly, even with shared descendants.

5. Suppose just $y$ is paired. Fix its partner to zero and fix $x=i$. The gate $G$ is constant. Delete it and the (different) unique child of the partner, which has a constant pin. The surviving $y$ is an unused singleton as its only child was $G$. Thus the initial charge $6+4$ of these blocks is fully recovered by $y$ alone. If both $x,y$ are paired, their pairs must be different by item 2. Set both partners to zero, then set $x=i$. Again $G$ is constant; remove it and at least one partner's child. Such a child differs from $G$, even if both partners shared it. Only $y$ survives from these blocks, unused. Its charge 10 and the two deleted gates pay for the two initial charges of 6. In either case the partner's child cannot reroute anything to $y$, since $y$'s only child was $G$.

These cases remain valid when a child is also another child's descendant. A substitution is simultaneous at all uses, and the proofs charge only explicitly distinct gates or rows. No previously unused singleton acquires a use in these restrictions. A promise that was imposed in item 4 persists until one of its two parameters is restricted. Since paired parameters never take part in a coordinate exchange, the affine embedding respects all remaining pair promises. Each of items 1–3 and 5 decreases $m$ without increasing $s$. Item 4 and elementary deletions decrease $s$ without increasing $m$. The coordinate exchanges are bounded by two per selection of $G$. Thus $m+s$, together with the finite number of pure gate-to-row moves between decreases, ensures termination. With $g<4n-10d$, $b>d$ holds at the end of every round, as shown by the charge. Termination supplies a constant space of dimension at least $d$, which contains a $d$-dimensional constant space, a contradiction. This proves the finite circuit bound.

## Returning to actual input length

For the fixed irreducible slice set $d=kt+r(r+1)/2+1$, an integer. For sufficiently large $n$, $1\le d\le N$, as the estimates below also show. Any full $n$-input circuit with $g$ gates becomes a circuit on $N$ data bits by feeding constants into its selector pins. All occurrences of a selector receive the same constant. One may keep all $g$ original gates and apply the finite theorem directly with constant sources. Its normalization already handles unary gates and the signed terminal bit. Equation (1) and that theorem, with $N$ in place of the local dimension $n$, give
$$
 g\ \ge\ 4N-10d=4n-(4k+10d).
$$

Let $z=k/n^{3/5}$. By (P), $z\to z_0=2^{2/5}$. Since $r=\lceil(n-k)/k\rceil\le n/k$ and $t\le\sqrt r+1$, we have
$$
 kt\le\sqrt{nk}+k,\qquad
 \frac{r(r+1)}2\le\frac{n^2}{2k^2}+\frac n{2k}.
$$
In particular $k=O(n^{3/5})$, $d=O(n^{4/5})$, and $N/n\to1$. More precisely,
$$
 \frac{4k+10d}{n^{4/5}}
 \le 10\sqrt z+\frac{5}{z^2}
        +14z\,n^{-1/5}+\frac5z n^{-2/5}+10n^{-4/5}.
$$
Every term on the right has its asserted limit over all integer lengths. Because $z_0^{5/2}=2$,
$$
 10\sqrt{z_0}+\frac5{z_0^2}
    =\frac{25}{2}\,2^{1/5}=:A.
$$
The remaining three terms tend to zero; thus $\limsup(4k+10d)/n^{4/5}\le A$. The strict rational comparison needed for the declared loss is
$$
 718^5-2\cdot625^5=84{,}402{,}824{,}318>0,
 \qquad
 2^{1/5}<718/625.
$$
Consequently $A<(25/2)(718/625)=359/25$. Take half the positive gap between $359/25$ and $A$. The definition of convergence applied to the entire displayed upper bound gives one integer threshold after which $(4k+10d)/n^{4/5}<359/25$ at every length. Increase it if necessary so that $k<n$ and $1\le d\le N$. The finite circuit bound then yields
$$
 C_{B_2}(f_n)\ge4n-\frac{359}{25}n^{4/5}
$$
for every integer beyond that fixed threshold. No numerical estimate of the threshold or effective prime-gap rate is needed for this eventual theorem.

Finally $n^{4/5}/n=n^{-1/5}\to0$. For each $\epsilon>0$, taking $n$ beyond the preceding fixed threshold and with $(359/25)n^{-1/5}\le\epsilon$ gives $C_{B_2}(f_n)\ge(4-\epsilon)n$. The finite projection is $(4,1)$ for this same family.
