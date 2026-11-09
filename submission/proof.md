# One uniform family selected by the improved padding-aware certificate

## Claim and material sources

For one total uniform P language at actual input length N we prove

    C_B2(F_N)>=4N-L_N,
    L_N/N^(4/5) -> c_star=10*2^(1/5).

The leading guarantee and finite projection are jointly(4,1). The field
prefix is input and is fully counted. The global exact-selection strategy,
its full-block reference and leading-certificate optimization are credited
to8eef1799-c27f-46f3-8e61-3c4e0d28d145 (immutable source commit
b4a12ebc5267f63705e1c277565e2779c088c286). Here selection uses the improved
padding-aware certificate with recovered prefix gates. This defines a NEW
fixed language; its parameters can differ from earlier fixed-selector
families. Comparisons below are between proved losses, not optimal circuit
sizes, and do not identify the different languages.

Every unfamiliar premise is supplied in algebra.md, elimination.md and
degree.md. These reproduce the material BK algebra and unrestricted charged
full-B2 reduction, including its degree-sensitive connectivity refinement.
No earlier verdict, direct sum, restricted circuit basis, advice or
conjecture is a premise. The old separate universal+8 theorem is not used.

## 1. Exact parameters and one evaluator on every input

For each nonnegative integer N enumerate every prime integer2<=k<N.
Put m=N-k, r=ceil(m/k), and require r>=3 and k>r+1. For every integer
2<=t<r put ell=m-(r-1)k, p=k-ell and

    H=r(r+3)/2, T=(t-1)k-t(t+1)/2+1,
    U=rk-r(r+1)/2, B_full=t(t+3)/2,
    H_c=H-[ell=1], T_c=T-min(p,t-1),
    U_c=U-(k-1)[ell<=2], B_c=B_full+t[ell<=2],
    D_c=H_c+max(T_c,floor(U_c/(t+1)))+1.

Brackets denote integer indicators. Call(k,t) admissible if B_c<=D_c<m,
and set

    L(N;k,t)=3k+8D_c-4B_c-3.

If no pair is admissible, output0 on every input of length N. Otherwise
choose the lexicographically least(L(N;k,t),k,t); denote its loss L_N and
parameters k,m,r,t,D_N,B_N. All tests and tie breaks use exact integers.
This is one algorithm for all N, not a family chosen for a requested bound.

Read the first k bits as the lower coefficients of the monic polynomial
Q(Z)=Z^k+sum_(j=0)^(k-1)b_j Z^j. Require Z^(2^k)=Z modulo Q and
gcd(Q,Z^2+Z)=1; if either fails output0. For an accepted Q use
K=F_2[Z]/Q. Zero-pad the remaining m bits to rk bits and interpret them
as field blocks a_1,...,a_r, each first bit its constant coordinate.
Set w_i=r+2-i, y_i=a_i^(2^(w_i)-1). Output the constant field coordinate
of sum_(|I|=t) product_(i in I)y_i. Evaluate it by E_0=1,E_1=...=E_t=0,
updating E_j<-E_j+y_i E_(j-1) in decreasing j=min(i,t),...,1.
Splitting subsets by membership of i proves this recurrence. Exponentiation
uses repeated squaring, with O(w_i) multiplications per block.

## 2. Prefix correctness and uniform eligibility

In an algebraic closure of F_2, Z^(2^k)-Z has 2^k distinct roots because
its derivative is 1. The roots are closed under addition, multiplication
and nonzero inversion, hence form a field of size 2^k. Frobenius z->z^2
permutes them and its k-th iterate is the identity. A Frobenius orbit of
length d gives an irreducible degree-d polynomial over F_2: its product
of linear factors is Frobenius-fixed, while a proper factor over F_2 would
require a proper invariant subset of the orbit. Conversely, an irreducible
factor's roots form one orbit. Thus factor degrees divide k.

Since k is prime, the orbits have size 1 or k. Only 0 and 1 are fixed,
and 2^k>2, so an irreducible degree-k polynomial exists for every such k.
No algorithm selecting it is needed in the uniform evaluator.

The first prefix test says Q divides the squarefree polynomial Z^(2^k)-Z;
the second excludes the factors Z and Z+1. Remaining irreducible factors
all have degree k. Since Q has degree k, it must itself be irreducible.
Conversely, every irreducible degree-k Q passes both tests. The test is
therefore exact. Compute its first part by k squarings modulo Q and its
second by polynomial Euclid, both in polynomial bit time.

All parameter loops have at most N² iterations, each using polynomial-time
trial primality and polynomial-size integer arithmetic. Selected k,r,t are
at most N. Schoolbook field multiplication costs O(k²), all powers use
O(r²) multiplications, and the symmetric recurrence O(rt). Polynomial gcd
and the Frobenius test are polynomial as above. A coarse O((N+1)^6) bound
suffices, including rejected inputs. No practical runtime at every large
length or algorithm producing an irreducible prefix is assumed.

Let A be this evaluator, L={x:A(x)=1}, W(N)=0, and let V(x,w) accept exactly
when w is empty and A(x)=1. Membership gives an accepted empty witness;
conversely every accepted permitted witness implies A(x)=1. Thus one fixed
polynomial-time verifier proves both NP directions. Witness bits are not
circuit inputs. Only a valid-prefix restriction is a disperser, not the
whole prefix language.

## 3. The improved certificate for every admissible pair

Use widths L_i=k for i<r,L_r=ell and weights w_i=r+2-i. Let
S_c={i:L_i>w_i}. Because k>r+1, all full blocks are eligible, and only
the last can be excluded. Hence S_c has at least r-1>=t elements.
The displayed closed forms are exactly

    H_c=sum_i min(w_i,L_i),
    T_c=sum of the t-1 largest max(L_i-w_i,0),
    U_c=sum_(i in S_c)(k+1-w_i),
    B_c=sum of the t smallest eligible weights.

Indeed only weight2 has reduced width: H loses one exactly when ell=1;
its capacity k-2 drops by p until the next full weight t+1 replaces it,
losing min(p,t-1) from T. It is excluded exactly when ell<=2, removing
k-1 from U and replacing weight2 by t+2 in the t selected weights.
In particular B_c<=H_c<D_c and B_c>=5.

For any affine space A in the padded payload image, algebra.md Section4's
triangular decomposition has dim(V_i)=d_i<=L_i. The last bound holds
because the kernel of projection onto the first r-1 blocks lies in the
ell-dimensional field-coordinate padding subspace. Put delta_i=d_i-w_i,
S={i:delta_i>0}. Then S subset S_c and sum_S delta_i>=dim(A)-H_c.
If |S|<t this sum is at most T_c. If every t-subset failed
sum_I delta_i+max_I d_i>k+1, the same averaging as algebra.md Section5
gives(t+1)sum_S delta_i<=sum_S(k+1-w_i)<=U_c. Both are impossible
when dim(A)>=D_c. A t-subset satisfies the same coefficient-isolation
hypotheses, and algebra.md Sections3/6/7 supply its nonzero coefficient.
Its binary weight is at least B_c, so Section9's Boolean-degree bridge
proves degree at least B_c on every such payload space.

Apply degree.md at n=m,d=D_c,B=B_c. Its unrestricted B2 theorem gives
C_B2(P)>=4m-8D_c+4B_c+3 for every valid-prefix payload P. All affine
rows remain charged. Terminal degeneracies, cyclic fair rows, fanout,
sharing and all16binary tables are included in the supplied proofs.
Recovering the k prefix gates as proved next yields

    C_B2(F_(N;k,t))>=k+C_B2(P)>=4N-L(N;k,t).

This holds for every admissible pair, hence for the selected branch.

### 3.1 Recover k prefix gates without additivity

Every irreducible Q of degree k>1 has Q(1)=1. Flipping any coefficient
bit makes Q(1)=0, so the resulting prefix is rejected on all payloads.
After any subset of prefix bits is fixed to a valid prefix b, every
remaining prefix bit is still essential: set the other free bits to b
and take a payload with P=1, which exists because degree(P)>=B_c>=5.
Flipping that bit changes output to0.

If a full-B2 circuit depends on x and its restriction x=c is not a
constant or literal, that restriction has a circuit with one fewer gate.
Normalize internal constant/unary gates by substitution, absorbing
complements into consumer tables or the terminal driving gate. Dependence
ensures a gate using x. Fixing x makes that gate constant or unary in
its other pin, also with repeated pins; remove it by the same rules.
If a terminal complement has no surviving driving gate, the restricted
function is a constant or literal, excluded here. Shared consumer tables
are changed in place and require no new gates. Further simplifications
can only lower size. Each successive prefix fixing retains P as a further
restriction, so cannot become a constant or literal. Applying this step
k times to any original circuit proves C(F)>=k+C(P), not a direct sum.

## 4. Admissibility and the upper limit at all lengths

Let b be the least integer with b^5>4N³, k_plus the least prime at least b,
r_plus=ceil((N-k_plus)/k_plus), t_plus=ceil(sqrt(r_plus)). Bertrand's
theorem gives k_plus=Theta(N^(3/5)),r_plus=Theta(N^(2/5)),
t_plus=Theta(N^(1/5)). Thus r_plus>=3,k_plus>r_plus+1,2<=t_plus<r_plus
eventually. The full-block D_full=H+max(T,floor(U/(t+1)))+1 is
O(N^(4/5)); D_c<=D_full, B_c<=H_c<D_c and m=Theta(N). This reference
is admissible at EVERY sufficiently large length, and its loss is
O(N^(4/5)). Hence L_N exists and is O(N^(4/5)) over all those lengths.
The finite theorem already proves the joint4N-o(N) guarantee.

For its exact leading loss use the classical prime number theorem
pi(x)~x/ln(x), NIST DLMF27.2.3, https://dlmf.nist.gov/27.2.E3 .
This familiar standard theorem is not reproved. For each fixed eta>0,
pi((1+eta)x)-pi(x)~eta*x/ln(x)>0 eventually. Thus the least prime at
least x is x(1+o(1)), an all-length relative-gap consequence, not a
short-interval estimate. It gives k_plus/N^(3/5)->a_star=2^(2/5).

More generally k/N^(3/5)->a>0 and t=ceil(sqrt(r)) give
D_c=r²/2+kt+o(N^(4/5)), B_c=O(N^(2/5)); the width changes are lower
order by the displayed closed forms. The loss divided by N^(4/5)
tends to8sqrt(a)+4/a². Its derivative is4/sqrt(a)-8/a³, so its unique
global minimum is at a_star, with value c_star=10*2^(1/5). The reference
therefore proves limsup L_N/N^(4/5)<=c_star.

## 5. Matching lower limit for this explicit certificate class

Consider any admissible sequence with loss L<=C*N^(4/5). Other sequences
cannot yield a smaller finite liminf without such a subsequence. Put
M=max(T_c,floor(U_c/(t+1))). Since B_c<=H_c and H_c>=H-1,

    L=3k+8H_c+8M-4B_c+5>=3k+4H_c+8M+5.

It follows that k=O(N^(4/5))=o(N),r=O(N^(2/5)). Since rk=N+O(k)~N,
k=Omega(N^(3/5)) and k/r tends to infinity. Uniformly for2<=t<r,
T_c=(t-1)k-O(t²+t)>=kt/4 eventually, so kt=O(N^(4/5)). Also
U_c=rk-r(r+1)/2-O(k)=N+O(k+r²)~N. The bound
L>=8U_c/(t+1)-8 forces t=Omega(N^(1/5)). Thus k=O(N^(3/5)); the
kt bound and k lower bound give t=O(N^(1/5)), and rk~N supplies
the matching lower bound on r. Consequently

    a=k/N^(3/5), b=t/N^(1/5)

range in a compact subset of positive reals, depending only on C.
Uniformly there T_c=kt+O(N^(3/5)),
floor(U_c/(t+1))=N/t+O(N^(3/5)), H_c=r²/2+O(r), and
B_c=O(N^(2/5)). Hence

    L/N^(4/5)=G(a,b)+o(1),
    G(a,b)=4/a²+8max(ab,1/b)>=4/a²+8sqrt(a)>=c_star.

This proves the lower liminf for every admissible parameter sequence,
not merely sequences whose parameter ratios converge. Together with the
reference it proves L_N/N^(4/5)->c_star. The minimizer is unique:
a=a_star,b=a_star^(-1/2). Compactness and uniqueness also give the
selected limits k/N^(3/5)->a_star and t/N^(1/5)->2^(-1/5).
This optimizes the LEADING VALUE OF THIS CERTIFICATE CLASS, not actual
circuit sizes, other proof methods, or all algebraic constructions.

## 6. Comparison with the supplied previous selector

Define the previous full-block admissible set by the same prime k,
r>=3,k>r+1,2<=t<r and B_full<=D_full<m. Let its lexicographic minimum
of(4k+8D_full-4B_full-3,k,t) have loss L_old and degree k_old. This
exactly identifies the declared8eef source's selection rule as a
mathematical comparator, not as an imported accepted premise.

Every such pair is admissible here, since D_c<=D_full and B_c<=H_c<D_c.
At that same pair its new certificate is at most its old one minus k.
Therefore at every sufficiently large actual length,

    L_N<=L_old-k_old.

This is an unbounded saving Omega(N^(3/5)). Indeed the same reference
gives L_old=O(N^(4/5)); its positive full-block ledger
L_old>=4k_old+4H+8M+5 forces r_old=O(N^(2/5)), k_old=o(N), and
r_old*k_old~N, so k_old=Omega(N^(3/5)). The argument compares proved
losses for two different fixed languages, not their optimal complexities.
Likewise the earlier fixed evaluator's admissible parameters are eventually
in our search, so its padding-aware certificate is an upper bound on L_N.

Finally L_N=o(N): for every epsilon>0 choose a threshold after which
L_N<=4epsilon*N. Then C_B2(F_N)>=(1-epsilon)4N at every integer length
past it, jointly proving the full declaration and projection(4,1).


## 7. A sharp secondary lower envelope of this certificate

Let a_star=2^(2/5), c_star=10*2^(1/5). For the same selected loss,

    L_N>=c_star*N^(4/5)-6*a_star*N^(3/5)-O(N^(2/5))      (8)

at all sufficiently large lengths, and

    liminf_N (L_N-c_star*N^(4/5))/N^(3/5)=-6*a_star.     (9)

There is an explicit infinite sequence of lengths attaining(8) up to
O(N^(2/5)). Consequently on those lengths the SAME fixed language has

 C_B2(F_N)>=4N-c_star*N^(4/5)+6*a_star*N^(3/5)-O(N^(2/5)). (10)

The secondary correction in(10) is an infinite-length assertion, not
an all-length upper estimate on L_N. Equation(8) lower-bounds this
explicit certificate loss, not actual circuit complexity or all methods.
No stronger prime-gap result is assumed.

### 7.1 Uniform lower envelope

The selected parameters have the compact scales proved in Section5.
Enlarge a containing compact rectangle to include(a_star,b_star).
Write a=k/N^(3/5), b=t/N^(1/5). Uniformly on that compact set, expand
one further order in the exact capacities:

 T_c=ab*N^(4/5)-a*N^(3/5)+O(N^(2/5)),
 floor(U_c/(t+1))=N^(4/5)/b
       -(1/b²+1/(2a²b))*N^(3/5)+O(N^(2/5)),
 8H_c=4*N^(4/5)/a²+O(N^(2/5)), B_c=O(N^(2/5)).

For the second expansion use r=N/k+O(1),
U_c=N-N²/(2k²)+O(k+r), and
1/(t+1)=1/t-1/t²+O(1/t³). The resulting error terms
N/t³, N²/(k²t²), k/t and r/t are all O(N^(2/5)).
Indicators for the last block and the floor change no stated order.

Taking the average of the two arguments lower-bounds their maximum.
Thus for fixed uniform error constant,

 L_N>=N^(4/5)*J(a,b)+N^(3/5)*K(a,b)-O(N^(2/5)),
 J=4/a²+4(ab+1/b), K=-a-4/b²-2/(a²b).

AM-GM and Section4's calculus show J>=c_star, uniquely at
(a_star,b_star)=(2^(2/5),2^(-1/5)). There K=-6a_star.
The Hessian of J there has positive first entry24/a_star^4
and determinant80, so it is positive definite. Smoothness, uniqueness
and compactness imply J-c_star>=c*distance² on the parameter compact
set for some c>0. On that same set K is Lipschitz, so
K>=-6a_star-M*distance for some finite M. Completing the square gives

 c*N^(4/5)*distance²-M*N^(3/5)*distance>=-M²*N^(2/5)/(4c).

This proves(8), and in particular the lower inequality in(9).

### 7.2 Attainment at explicit prime-indexed lengths

Take ANY sequence of primes k tending to infinity. Set

 t=floor((k/2)^(1/3)), r=floor(t²+t/4+1/2), N=k(r+1).

Then k=2t³+O(t²), r=t²+t/4+O(1), m=kr and ell=k: there is no
padding. The reference pair(k,t) is admissible eventually because
k>r+1,2<=t<r and D_c=O(t^4)=o(m), while B_c<=H_c<D_c.
The selected family remains the original global-minimum family; this
construction selects lengths, not a new algorithm or per-length advice.
Infinitely many primes suffice; PNT rates or short intervals are unused.

Put x=sqrt(N/k)=sqrt(r+1)=t+1/8+O(1/t). Then
k=2x³+O(t²), t-x=-1/8+O(1/t). The exact full-width capacities give

 T_c=kx-9k/8+O(t²), U_c/(t+1)=kx-9k/8+O(t²),
 H_c=x^4/2+O(t²), B_c=O(t²).

For the U expansion, its numerator is
k(x²-1)-(x²-1)x²/2; divide by x+(t-x)+1 and use
x³=k/2+O(t²). The floor costs at most one. Hence the branch loss is

 L_ref=4x^4+8kx-6k+O(t²)
      =N^(4/5)*(4/a²+8sqrt(a))-6k+O(N^(2/5)),
 a=k/N^(3/5)=a_star+O(1/t).

The smooth function 4/a²+8sqrt(a) has derivative zero at a_star.
Its excess over c_star is O(1/t²). Since N~2t^5 and
k=a_star*N^(3/5)+O(N^(2/5)), this proves

 L_ref=c_star*N^(4/5)-6a_star*N^(3/5)+O(N^(2/5)).

The global minimum L_N<=L_ref. Combining with(8) proves two-sided
attainment up to O(N^(2/5)) on these lengths and(9). The finite circuit
inequality then gives(10). No optimal circuit upper bound is inferred.
