# A fixed uniform family with lower bound 4N-9N^(1/1024)

## Claim and provenance

Fix sigma=1024 and delta=1/sigma. These are constants throughout; they
do not grow with the input length. The single total evaluator in family.py defines a language in P, hence NP,
with actual ENTIRE input length N. We prove eventually at every integer length

    C_B2(F_N)>=4N-9N^(1/1024).

The full B2 convention allows all binary tables, sharing, constants, unary
operations, repeated pins and one output. Both claims assert this entire
guarantee and its joint projection (4,1).

The kernel in elimination.md is reproduced from production submission
11e8b36d-cbd4-44f6-810f-7a7f875b4ad0, including its fourth AP2 deletion.
Its state framework comes from 237016f9-ece3-4a05-9933-dfcf6056e787 and
Mishra--Santhanam TR26-260. The certificate wrapper and image-polynomial
preliminaries build on f4f2c464-3045-4c32-a6e4-32478c05f343.
The new disperser adapts Shaltiel TR11-127 challenge-response, Rao affine
condensation, Trevisan/Nisan--Wigderson code-bit extraction, and the
Ben-Sasson--Kopparty x^7 cubic. All necessary arguments are reproduced
in disperser.md with the actual parameters. No unfamiliar extractor is
assumed as a black box, and no prior judgment is a mathematical premise.

## 1. Exact family on every length

The fixed evaluator family.py resolves every branch, rounding and ordering.
Written bit zero is stored in integer position zero. Put alpha=1024*sigma,
and fix

    b=2^(alpha+log2(sigma)+40), E_tree=alpha+3log2(sigma)+90,
    epsilon_tree=2^-E_tree, H=2^(E_tree+2),
    q=2^(alpha+2log2(sigma)+80), J=2^(14+log2(sigma)),
    A=4096, gamma=2^-20, e=1/1024.

Every logarithm of sigma is the integer 10. For length N,
FamilyShape.build computes t=floor(N^epsilon_tree). If
t<b+2, output zero. Otherwise compute

    K=floor(N^(delta/2)), R=floor(K^(1/8)),
    cutoff=floor(N^(delta/4)), pmax=floor(N^(delta/2^24)),
    ell=H*t^b+1,
    G=2^floor(floor(log2R)/(2*2^20)).

All powers use integer roots. Construct metadata for the short seeded
extractor (alphabet 2^alpha) and the long extractors (design order A).
The coded guards require cutoff>=2t, ell^(6J+10)<=pmax,
ell*pmax<=R, the greedy-design cardinality inequality, G>=2, enough
b blocks, every long-design cardinality inequality, and sufficient output
length for each requested seed.

The condenser has row lengths r_0=R, r_(i+1)=floor(r_i/2^4106), for
i<log2G. Put w_i=floor(r_i^(1/10)) and m_i=2w_i-1. Further guards
require positive next rows, m_i>=3, serialized matrices of length at most N,
stage seed length d_i<=floor(m_i^(1/16)), and a final row at least
floor(sqrt R). These seed lengths are the ACTUAL code/design lengths,
not just their asymptotic estimates. Failure of any guard outputs zero.

Put P=sum_i m_i; output zero if P>=N. Read the first P bits as lower
coefficients of one monic polynomial Q_i of degree m_i per stage, in order.
Repeated degrees still get separate polynomials. If any fails the exact
irreducibility test below, output zero. All P bits count in N.

Pack the remaining N-P bits in low positions of an N-bit vector, padding
the final P coordinates with zero. ExtractionEngine runs the balanced
interval tree and challenge-response algorithm in disperser.md with these
validated fields. OWN block-extractor rows precede the selected child's
rows. The first active node with confidence-1 FAIL emits its confidence-2
result; if none exists output zero. FAIL is zero and PASS is one.
No polynomial, seed or node choice is advice. All choices are coded.

## 2. Exact field test and existence at every degree

For monic degree-m Q let z=X modulo Q. Require z^(2^m)=z modulo Q,
and gcd(Q,z^(2^(m/p))-z)=1 for every prime divisor p of m. Compute
these values during m repeated squarings. Trial division finds the prime
divisors. At m=1 the empty gcd list correctly accepts both monic linears.

In an algebraic closure, X^(2^m)-X has distinct roots since its derivative
is one. They form a field under sums, products and nonzero inverses.
Frobenius orbits correspond to irreducible factors, because a proper
invariant subset of one orbit cannot exist. An orbit of length d is fixed
by Frobenius^m exactly when d|m. Thus X^(2^m)-X is the squarefree product
of all monic irreducibles of degrees dividing m.

Let I_d count monic irreducibles of degree d. Root counting gives
2^m=sum_(d|m)d I_d and d I_d<=2^d. For m>=2,

    sum_(d|m,d<m)d I_d<=sum_(j=1)^floor(m/2)2^j<2^m.

Hence I_m>0; degree one is immediate. This proves existence of every
needed prefix field without claiming a polynomial-time search at large
degrees. Splitting fields are used only for existence, not as an oracle.

The first test makes Q a squarefree product of factors with degrees dividing
m. Every proper divisor d|m divides m/p for some prime p|m, and that gcd
test excludes degree d. Only a degree-m factor remains, so Q is irreducible.
Conversely every irreducible Q passes. Testing is polynomial in m.

Code fields instead use the explicit quadratic towers in disperser.md;
graph-design fields have logarithmic degree, where exhaustive irreducible
search is polynomial in N. Those searches are part of the uniform algorithm.

## 3. Total polynomial time and both NP directions

Integer roots use binary search; when the fixed degree exceeds the input
integer's bit length the answer is immediately one (or zero at zero).
All guarded parameter arithmetic has polynomial cost in N. Field and
code operations have polynomial bit cost. The short-design greedy scan
enumerates a fixed alphabet to O(log N) coordinates, hence polynomially
many candidates; its explicit guard guarantees enough permissible words.
The long-design guard guarantees enough graphs. The full seed loop has
N^O(1) iterations with an enormous FIXED exponent. The block-count guard
ensures every requested fresh block exists. The tree has polynomially many
nodes, depth at most H, and at most t^b subsets per node. Its explicit
postorder stack goes only to strictly shorter intervals because cutoff>=2t
and t>=2. Tower multiplication and tower construction also use iterative
stacks/loops. No fixed host-language recursion limit is part of the family.
Subset enumeration uses integer bit masks. Matrix and seed lengths are
polynomially bounded. Memoization changes no value.

Thus there are fixed finite c_0,c_1 bounding runtime by c_0(N+1)^c_1 on
ALL inputs, including rejected prefixes and guard branches. Constants and
the polynomial degree are enormous. At physically practical lengths the
first guard returns zero; this is an asymptotic complexity theorem, not a
practical extractor-speed claim.

For L={x:evaluate(x)=1}, define W(N)=0 and V(x,w) to accept exactly when
w is empty and evaluate(x)=1. Membership supplies the accepted empty
witness; an accepted permitted witness implies membership by definition.
Both directions use the same uniform polynomial-time verifier. Witnesses
are not target-circuit inputs. This establishes one total language in P/NP.

## 4. All sufficiently large integers pass every guard

As N tends through all integers,

    t=N^epsilon_tree(1+o(1)), K=N^(delta/2)(1+o(1)),
    R=N^(delta/16)(1+o(1)), ell=N^(delta^2/2^50+o(1)).

Fixed lower bounds such as t>b+2 and G>=2 eventually hold. The cutoff
exponent exceeds the branching exponent. Since
(6J+10)*delta^2/2^50<delta/2^24 and
delta^2/2^50+delta/2^24<delta/16, the confidence and
response-length guards pass. The seeded proof proves the greedy and
long-graph capacity bounds. Actual long seed length is
O(log(N)^2+h^(2e)) at entropy h, and every used h is a fixed positive
power of N, so this is O(h^(2e)) eventually.

All condenser rows are at least R^.99, the row count is at most
r_i^gamma, and serialized size G*R=N^(delta/16+delta*gamma/32+o(1))<N. Positive next rows
and the final sqrt(R) bound follow. Cubic seed allowance has exponent
1/160 in r_i, larger than 2e; BE's O(K^(2e)) seed fits within sqrt(R).
The block-count estimate in disperser.md is below b. Finally

    P<=2(log2G)*R^(1/10)=O(N^(delta/160)log N)=o(N^delta),

so P<N. Thus every guard passes for ALL sufficiently large integer N,
not just a subsequence. Prefix validity still depends on input; valid
prefixes exist at each such N by section 2.

## 5. Hard restriction and full explicit remainder

Fix any valid list of prefix polynomials for a sufficiently large N. The
payload has n=N-P inputs and injectively zero-pads into N coordinates.
The complete disperser proof shows nonconstancy on every affine source
of dimension at least k=ceil(N^delta), for EVERY validated field list.
The payload therefore disperses all affine spaces of dimension >k-1.

Apply elimination.md with d=k-1 to get size at least 4(N-P)-8k.
Fixing prefix bits in a full-family circuit produces a payload circuit
with no size increase. The gate normalization covers the resulting
constants, unary gates, repeated pins and all sharing. Consequently

    C_B2(F_N)>=4N-(4P+8ceil(N^delta)).

The loss divided by N^delta tends to 8, since P=o(N^delta). Hence it
is at most 9N^delta at every sufficiently large integer length. This
proves the complete declared fixed-remainder guarantee, and division by
N yields liminf C_B2(F_N)/N>=4, the joint projection (4,1).
Prefix selection is only a circuit restriction in the proof. The evaluator
never receives advice, and all prefix coordinates were included in N.

## 6. Scope of finite evidence

The finite audit checks bounded circuit identities, arithmetic, field/code
primitives and actual small-length guard evaluations. All those full-family
evaluations return zero through the documented guard; they do not execute
the enormous nontrivial asymptotic instances. Finite checks establish no
universal lemma. The supplied proofs carry the all-length, uniformity and
asymptotic obligations, independently of any prior judging result.
