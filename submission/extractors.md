# Locally evaluable linear binary code with a polynomial list bound

This is the standard polynomial-evaluation/trace-concatenation idea, with only a list SIZE bound required; efficient list decoding is not used by the extractor evaluator.

Let n>=1 and let0<rho<=1/4 be rational. Choose a prime integer s with q=2^s>1024(n+1)/rho^6. Use the deterministic prime-degree binary field construction to form F_q. Pack the n input bits, followed by at most s-1 zero bits, into coefficients of a polynomial p in F_q[T]. Its degree is at most D=max(1,ceil(n/s)-1). In particular D<=n+1 and D<q. Define the binary word indexed by(a,b) in F_q^2 as

 C_p(a,b)=Tr_(F_q/F2)(b*p(a)).

This encoding is binary linear and injective. Evaluation at any index takes polynomial time in n and s, by Horner evaluation and Frobenius trace. No array of q^2 entries is constructed. Prime search and field construction are polynomial in s for the chosen prime degree, assuming the supplied coarse prime-interval and constructor arguments. One may choose the least prime s exceeding the displayed binary-log threshold; then s=O(log(n+1)+log(1/rho)), eventually and globally after enlarging a fixed finite constant. Thus an index uses2s bits.

**List bound.** Any binary received word has agreement at least1/2+rho with at most q distinct encoded messages.

For a received word r(a,b), put h_a(b)=(-1)^r(a,b), and for v in F_q define c_a(v)=E_b h_a(b)*(-1)^Tr(bv). Agreement implies E_a c_a(p(a))>=2rho. Since each summand is at most1, the fraction of a for which c_a(p(a))>=rho is at least rho: if a fraction theta had this property, the average would be at most theta+(1-theta)*rho, and theta<rho would make that smaller than2rho. For each a, define V_a={v:c_a(v)>=rho}. Orthogonality of the q trace characters gives sum_v c_a(v)^2=1, so |V_a|<=rho^-2. Every good message therefore satisfies p(a) in V_a on at least rho*q field elements.

Set T=floor(rho*q/2). Choose a nonzero bivariate polynomial Q(A,Y) of weighted degree at most T, giving A weight1 and Y weightD, that vanishes on all pairs(a,v) with v in V_a. There are at most q/rho^2 homogeneous linear constraints on its coefficients. The monomials A^iY^j with 0<=i<=floor(T/2), 0<=j<=floor(T/(2D)) all have weighted degree at most T, so their count is at least T^2/(4D). Since rho*q is much larger than1, T>=rho*q/4. The number of available monomials is consequently at least rho^2*q^2/(64D)>q/rho^2, by q>1024(n+1)/rho^6 and D<=n+1. Linear algebra therefore gives such a nonzero Q.

For any good message, Q(A,p(A)) has degree at most T and vanishes on at least rho*q>T points. It is the zero polynomial. Polynomial division by the monic polynomial Y-p(A) over F_q[A] shows that Y-p(A) divides Q(A,Y). Distinct messages give distinct monic factors; viewing Q over F_q(A)[Y], the number of such factors is at most deg_Y Q<=T/D<q. If deg_YQ=0 there are no good messages. This proves the stated (looser) bound q in every case. No decoding algorithm, factor-extraction complexity or probabilistic construction is assumed.

A numerical detail to retain in later use: agreement and correlation differ by a factor2. The threshold at each row is rho, while global correlation is2rho. The list size costs at most s bits of advice. The code's evaluation complexity depends on the FIELD BIT LENGTH s, not on its exponentially larger number q of elements.


# Linear strong extraction from a design and a list-bounded code

This section proves the standard Trevisan reconstruction principle. The encoder lemma is the corresponding section of extractors.md. The following sections supply the explicit designs, parameter regimes and iteration.

Let C be an injective binary linear encoding of n-bit messages whose code positions are indexed by ell bits. Suppose at most L codewords agree with any received word on a fraction at least1/2+rho, where rho=epsilon/(2m), 0<epsilon<=1/4 and m>=1. Let S1,...,Sm be specified ell-element subsets of a d-element seed universe, ordered internally in a specified way. Put

 B=max_(1<=i<=m) sum_(j<i) 2^|Si intersect Sj|.

Define Ext(x,y)_i=C(x)[y restricted to Si]. If a source X on n bits has min-entropy at least

 k >= d+B+log2 L+log2(4m/epsilon),

then (Y,Ext(X,Y)), for independent uniform d-bit Y, is within total variation epsilon of uniform(Y,U_m). For each fixed seed the map x->Ext(x,y) is binary linear.

Proof. Suppose the total variation exceeds epsilon. A Boolean test D(y,z) distinguishes the two distributions with probability gap exceeding epsilon. Insert the m hybrids that replace suffix coordinates of Ext by independent uniform bits. One adjacent pair has signed gap exceeding epsilon/m; complement D if necessary to make it positive. Fix this coordinate i and let d_b(y,z_<i) be the conditional acceptance probability after assigning coordinate i the bit b and randomizing the later output bits. The difference between the adjacent hybrids is exactly

 (1/2)*E[(-1)^Ext_i * (d_0(Y,Ext_<i)-d_1(Y,Ext_<i))].

Consider a randomized predictor that outputs0 with probability(1+d_0-d_1)/2 on input(y,z_<i). Its prediction advantage over1/2 is precisely the displayed difference. Averaging its finite table of random choices permits fixing one deterministic predictor P with advantage alpha>epsilon/m over X,Y. No efficiency is required of P; it is fixed once and does not depend on the individual source sample.

Write a(x)=Pr_Y[P(Y,Ext(x,Y)_<i)=Ext(x,Y)_i]-1/2. It lies in[-1/2,1/2] and has expectation alpha. The set G of messages with a(x)>=alpha/2 has X-probability at least alpha/2, hence greater than epsilon/(2m). Indeed if its probability were below alpha/2, the expectation would be at most (alpha/2)*(1/2)+(1-alpha/2)*(alpha/2)<alpha. This intentionally uses a loose bound.

For each x in G, fix the seed coordinates outside Si to an assignment z for which the conditional prediction advantage remains at least alpha/2. Supply z as at most d advice bits. For each j<i, with these outside coordinates fixed, the previous code bit C(x)[Y|Sj] is a function only of the coordinates in Si intersect Sj. Supply its whole truth table as at most2^|Si intersect Sj| bits, with a fixed lexicographic table order. All tables together use at most B bits. Tables may be padded to a common fixed advice length.

Given any such advice, each t in{0,1}^ell determines the entire seed by combining t with z, and the tables determine every preceding output bit. Running the fixed predictor P therefore defines one received word r(t). For the particular advice belonging to x, this word agrees with C(x) on at least1/2+alpha/2>1/2+epsilon/(2m)=1/2+rho of positions. By the code hypothesis at most L different messages have this property for any fixed advice. There are at most2^(d+B) advice strings. Thus

 |G|<=2^(d+B)*L,
 Pr[X in G]<=2^(d+B)*L*2^-k<=epsilon/(4m),

contradicting the previously established strict lower bound epsilon/(2m). This proves strong extraction, with the seed included in the tested distribution. The seed was not assumed independent after fixing message-specific advice; independence is used only in the original expectation and hybrid argument.

Finally each output code bit is linear in x. Fixing a seed fixes its index, so all output bits jointly form a linear transformation of x. All explicit evaluators need only the locally computable encoder and the specified design; the hypothetical predictor, advice, received words and list are proof devices and are never computed by the family.


# A supplied logarithmic-seed linear strong extractor

This lemma uses Shaltiel--Umans (2005) interleaved-curve reconstruction with the seed retained. It strengthens its finite field bound by multiplicity-weighted interpolation, the Guruswami--Sudan/Koetter--Vardy soft-decision technique (Koetter--Vardy, IEEE TIT49(11),2003, DOI10.1109/TIT.2003.819332). Every required multiplicity and reconstruction argument is proved below; no decoding theorem is assumed. Keeping the seed establishes the strong version needed here. This is the constant-error linear strong extractor used by family.py.

All logarithms below have base2. Let n,K,m,d,h,s be positive integers, q=2^s, Q=q^d, rho=1/(32m), and assume

    h>=2, h<q-1, 2m<2Q,
    q>2^30*h*d*s/rho,
    K>=400*m*h*d*s^2+ceil(log(1/rho))+3.

Let F be the explicitly constructed field with q elements. Encode x in {0,1}^n injectively and F2-linearly as a d-variate polynomial xhat over F, of total degree less than h. For example assign the bits as coefficients0 or1 of the first n monomials in increasing total degree and then lexicographic exponent order, provided binom(h+d-1,d)>=n. Let A be multiplication by a primitive element of the degree-d extension of F, on its standard d-dimensional F basis. Thus A is invertible and its powers act transitively on F^d minus0. The finite field, extension, primitive element and monomial choices are deterministic and specified below.

Define E_q(x,v)_i=xhat(A^i v), for v in F^d and1<=i<=m. Define the binary map E(x;(v,b))_i=Tr_F/F2(b*E_q(x,v)_i), for independent uniform b in F. It has seed length (d+1)s and output length m. We prove that it is a strong extractor for min-entropy K and error1/4, and that it is linear in x for each fixed seed. Neither proof predictors nor their decoders are computed by its evaluator.

## Seed-aware q-ary prediction is impossible

Fix1<=i<=m. For any deterministic real weights omega(v,z_<i;c)
in[0,1] with sum_c omega^2<=1 at every prefix, put
S_x=E_v omega(v,E_q(x,v)_<i;E_q(x,v)_i) and
C_x=E_v sum_c omega(v,E_q(x,v)_<i;c)^2. We claim

 E_X(S_x-C_x/2)<=2rho.                              (SC)

Otherwise call x good if S_x-C_x/2>=rho. Its X-mass exceeds rho:
this quantity is<=1, while if good mass<=rho its mean is at most
rho+(1-rho)*rho<2rho. Every good x has S_x>=rho,C_x<=2S_x.
Weights are fixed finite proof devices, possibly depending on X but
not on an individual message; they are never evaluated by the family.

Fix a good x and put r=100ds and B=(2r-1)(h-1), so B+1<=2rh. The field-size condition ensures2r<q. Choose2r distinct uniform parameters t1,...,t_(2r) in F, and independent uniform y1,...,y_(2r) in F^d. Interpolate vector polynomials p1,p2 of degree less than2r, with

    p1(tj)=yj for every j;
    p2(tj)=A*yj for j<=r, and p2(tj)=yj for j>r.

For0<=j<Q put P_(2j+1)=A^j p1 and P_(2j+2)=A^j p2. Each P_a is a uniformly random vector polynomial of degree less than2r, independent of the entire tuple of parameters. Indeed, for each fixed distinct parameter tuple, interpolation is a bijective linear map from the uniform vector values to the coefficients; the same is true after applying A to some values or to all coefficients. Its values at any2r distinct arguments are independent uniform vectors, by the invertibility of the corresponding evaluation map.

Successive curves P_(a-1),P_a agree on r distinct uniform parameters S_a independent of P_a. If a is even take the last r parameters; if a is odd and at least3 take the first r. The identities follow directly from the two interpolation conditions. Independence follows from the preceding coefficient argument. These assertions are marginal assertions for each a; no independence between different curves is required.

Give a reconstruction procedure the coefficients of xhat composed with the INITIAL2m curves. Each composition has degree at most B, so this advice has at most

    a0=2m*(B+1)*s<=4mhrs=400mhds^2

bits. Pad to the integer upper bound a=400mhds^2. Giving all2m startup curves avoids any boundary problem when the weight-prefix index is i=1. Shared random curves are random coins of the reconstruction, not advice bits.

Suppose a>2m and all earlier compositions have been learned. For u=P_a(t), its previous q-ary output values at seed v=A^-i u are

    xhat(A^(j-i)u)=(xhat composed with P_(a-2(i-j)))(t), 1<=j<i.

All these curve indices are positive and earlier. Hence omega defines weights omega_a(t,c) from that known prefix AND the seed A^-i P_a(t). With ideal true prefix polynomials, these weights depend only on P_a, fixed x and omega, not on the overlap parameter set S_a.

Fix a good x and abbreviate S=S_x>=rho,C=C_x<=2S. On a random
P_a, the true-symbol weights W_t and norm costs
C_t=sum_c omega_a(t,c)^2 are2r-wise independent variables in[0,1],
with respective means S,C. The independence follows from the2r-wise
uniform curve evaluations; each is a function of that evaluation.

For2r-wise independent Y_t in[0,1] of mean p, expand the centered
2r-th moment. Singleton indices vanish, and E|Y_t-p|^e<=p for e>=2
because |Y_t-p|<=1 and Var(Y_t)<=p. For j distinct indices the absolute
product expectation is<=p^j. Choosing representative positions and
assigning remaining positions gives at most
binom(2r,j)*j^(2r-j)<=4^r*(2r)^(2r-j) partitions. Hence, with mu=qp,

 E[(sum Y_t-mu)^(2r)]<=sum_(j=1)^r4^r*(2r)^(2r-j)*mu^j.

If mu<=mu_star and mu_star>=2r, this is at most
r*4^r*(2r)^r*mu_star^r<=(16r*mu_star)^r. Apply first to W with
mu_star=Sq to obtain Pr[sum W<Sq/2]<=(64r/(Sq))^r.
Apply to C_t with mu_star=2Sq, even when qC<2r. Since qC<=2Sq,
Pr[sum C_t>4Sq]<=(8r/(Sq))^r. The field guard gives Sq>=rho q>2r
and64r/(rho q)<6400/(2^30 h)<2^-16. Each failure is<Q^-1600.
Thus except with probability2Q^-1600, score>=Sq/2 and cost<=4Sq.

Let M=ceil(4/rho)=128m, a fixed integer known to the decoder. Assign
integer multiplicities mu_tc=floor(M*omega_a(t,c)). Vanishing to
multiplicity mu at(t,c) means the shifted polynomial has no monomial
of TOTAL degree<mu. This imposes mu(mu+1)/2 homogeneous LINEAR
constraints in any characteristic. Their ACTUAL count Kcost is known
from the weights; it is<=sum mu_tc^2<=M^2 sum C_t<=4M^2Sq on the
good event. Integer score on the true polynomial is at least
M Sq/2-q>=M Sq/4, since MS>=4.

Choose T=1+ceil_sqrt(4B*Kcost), COMPUTED FROM ACTUAL Kcost, not S.
Monomials with exponents<=floor(T/2),floor(T/(2B)) have weighted
degree<=T for weights1,B; their number>T^2/(4B)>Kcost. Linear algebra
gives nonzero bivariate R meeting the multiplicities. On the good event,

 T<=2sqrt(B*Kcost)+2<=4M sqrt(BSq)+2<M Sq/8.

Indeed q>2^30hds/rho and B<200hds imply q>4096B/S and q>8;
the two terms on the preceding right side are each<M Sq/16.
The true z(t)=xhat(P_a(t)), degree<=B, gives a composition R(t,z(t))
of degree<=T with root multiplicities summing>=M Sq/4>T, so it is zero.
For completeness shift at(t0,z(t0)): every monomial has total degree
at least mu, and z(t)-z(t0) is divisible by t-t0, so the composition
is divisible by(t-t0)^mu. Distinct t give coprime factors; this is valid
in characteristic2. Monic division then gives Y-z(Tvar) as a factor.
There are at most deg_Y R<=T/B<Mq/8<q^2 candidate z, since M<q.
The decoder can reject lists larger than q^2, without losing a good
case. Choose the first nonzero interpolation solution deterministically;
no efficient interpolation or factorization is needed or called.

The candidate list depends only on P_a and fixed x and omega. Test these candidates using the r known values at the overlap with P_(a-1). A false candidate differs from the true polynomial on all but at most B field points. Since S_a is a uniform distinct r-tuple independent of P_a, the probability that it passes every test is at most (B/q)^r: sampling without replacement gives a product of factors at most B/q. Taking a union over fewer than q^2 candidates bounds the disambiguation failure by

    q^2*(B/q)^r < 2^(2s)*(2^-22)^r <=2^(-2198ds).

Union bounds over at most2Q steps show that all learning steps succeed with probability at least1-2Q*(2Q^-1600+2^(-2198ds)), exceeding3/4 under the displayed guards. This is an unconditional union bound on the ideal learning steps; it does not condition the uniform-curve assertions on earlier success.

Also p1(0) is uniform in F^d and is nonzero except with probability1/Q<1/4 (the field guard makes Q much larger than4). If it is nonzero, its A orbit, occurring among the odd curves, covers every nonzero vector. Hence with probability greater than1/2 reconstruction obtains xhat on all nonzero vectors. These evaluations uniquely determine a polynomial of total degree less than h<q-1. Indeed their difference vanishes at all nonzero points; on any line through0 there are q-1 roots and degree less than q-1, giving the value0 at the origin as well. A polynomial with individual degrees less than q vanishing on all of F^d is zero, by iterated univariate root bounds. The injective encoding then recovers x.

For each good message there is the specified a-bit advice, and reconstruction with that advice succeeds with probability at least1/2 over shared coins. Therefore the X-probability of messages reconstructible from SOME a-bit advice, averaged over those coins, exceeds rho/2. For any fixed coins the deterministic procedure has at most2^a possible outputs. Since Hinf(X)>=K, that probability is at most2^(a-K)<=rho/8 by K>=a+ceil(log(1/rho))+3. This proves(SC). Finite real weight tables cause no oracle issue: their floored multiplicities and finite interpolation choices define one fixed unbounded decoder, independent of the individual advice/message.

## Binary conversion with the seed retained

Suppose the joint output with uniform seed Y has total variation
>epsilon=1/4 from(Y,U_m). Use natural-log relative entropy D. The
elementary inequality D>=2TV^2 follows by grouping a distinguishing
event and its complement (log-sum/Jensen), then observing that binary
relative entropy as a function of its first probability has minimum0
at equality and second derivative1/[p(1-p)]>=4. This also covers
endpoints by continuity. Thus D>2epsilon^2.

The chain rule, obtained by factoring joint conditional probabilities,
gives D=sum_i E D(next bit|Y,prefix || U_1). Put
b_i=E[(-1)^nextbit|Y,prefix]. Since log(1+x)<=x,
D(bit||U_1)=E log(1+b_i*(-1)^bit)<=b_i^2. Therefore some fixed
coordinate i has c=E b_i^2>2epsilon^2/m=4rho.

Let H(v,b,zbin_<i) be this CONDITIONAL MEAN sign, with value0 on
unreachable prefixes, rather than rounding it to a sign. This is a
fixed function in[-1,1]. For a FIELD prefix z and c_symbol in F set

 gamma_(v,z)(c_symbol)=E_b H(v,b,(Tr(bz_j))_(j<i))*(-1)^Tr(b*c_symbol).

Trace characters are orthonormal: s is odd,Tr(1)=1, and multiplication
by a nonzero field element preserves uniformity. Parseval gives
sum gamma^2=E_b H^2<=1. Define omega=max(0,gamma). Its mean true-symbol
weight is>=E[H*(-1)^nextbit]=E H^2=c by conditional expectation.
Its average norm cost is<=E H^2=c by Parseval, averaging actual field
prefixes over X,v. Hence E_X(S_x-C_x/2)>=c/2>2rho, contradicting(SC).
This proves strong error<=1/4. The real conditional means and weights
are proof devices; no posterior, decoder or distribution is an input
to the explicit extractor evaluator.

Every output is an evaluation linear in the message coefficients, followed by multiplication by the fixed b and binary trace. Thus fixing the complete seed fixes a binary linear map on x.

## Explicit uniform parameters and deterministic fields

Fix one public rational beta=pRate/qRate in[3/20,53/100], with
positive pRate,qRate<=100ceil(log2(n+1)). Every invocation, helper and
cache key carries this context; it never depends on source bits or epsilon.
Check its size before variable powers. Invalid numerical parameters take
the specified fallback. Fractions are reduced by gcd, and r^(u/v) is
implemented by the exact integer root of r^u, not floating arithmetic.

Put epsilon0=1/10, bSU=beta, aSU=1/(1-beta),
seedguard=1+aSU*beta+epsilon0, c=(9999/10000)*(1-beta)/seedguard,
and group=1-beta. The regime labels zeta=beta/2 and slice=beta
are limiting logarithmic rates, NOT finite pure-power error or width bounds.
The actual widths below include logarithmic losses. No fixed positive
power loss is imposed on either rate. The SU output remains
floor(floor(K^beta)/(L+1)^12); its advice has no positive exponent gap.

cap(beta) is2,3,4,5,6 on the consecutive bins ending at
19/100,3/10,2/5,12/25,53/100. Substitution at the endpoints gives
seedguard<cap-1/20. The full-seed coefficient is at most2.127660,
its guard at most2.227660, with gap epsilon0. Also
c*seedguard/(1-beta)=.9999 and .1<c<.8.
All rational parameters and field searches are deterministic.

Let L=ceil(log2(n+1)), W=log2 K, w=floor W>=1,
 m=floor(floor(K^beta)/(L+1)^12), d=ceil(aSU*L/w), h=d*ceil_root(n,d)+1.
There are binom(h+d-1,d)>=((h-1)/d)^d>=n monomials of degree<h,
so the bit-coefficient encoding used above exists. Put
C0=2^30*h*d*(32m). Starting at the least odd prime>=bit_length(C0),
visit odd primes in order until2^s>C0*s, or reject when s>(L+1)^2.
This returns the least admissible odd prime: smaller starting degrees
have2^s<=C0, and2^r/r is increasing for integer r>=2. Set q=2^s,Q=q^d.
Require exactly the finite guards

 s<=(L+1)^2, q>C0*s, Q<=(n+1)^cap, h^d<=(n+1)^cap,
 h<q-1, m<Q, K>=400*m*h*d*s^2+ceil(log2(32m))+3.

Require m>=1. The finite predictor proof applies to every returned tuple.
A rejected tuple is not used as a primitive. No eventual estimate
replaces any finite check.

Here is the complete uniform eventual calculation for
L^(2/3)<=W<=L^(19/20). Uniformly, d=Theta(L/W), d tends to infinity,
W/log L tends to infinity, d log L=o(L), and

 log2 h=L/d+O(log L)=(1/aSU+o(1))*W,
 log2 m=beta*W-12log2(L+1)+o(1)=(beta+o(1))*W,
 dW=(aSU+o(1))*L.

For T=log2(C0*(L+1)^2)+O(1), the relative-prime-interval theorem gives
its least odd prime s_old=(1+o(1))*T, uniformly since T>=c0*W tends
to infinity for one positive fixed c0. Thus
s_old=O(W)<=(L+1)^2 and2^s_old>C0*(L+1)^2>=C0*s_old eventually.
The returned least degree exists with s<=s_old. Conversely q>C0*s
implies s>log2 C0. These two bounds give

 s=(1/aSU+beta+o(1))*W=(1+o(1))*W,
 (d+1)s=(aSU+o(1))*L=(1+aSU*beta+o(1))*L,
 d log2 h=L+O(d log L).

All errors are uniform: floor/ceil errors are O(1/W)+O(1/d), and
log L/W tends uniformly to zero; the leading prime threshold coefficient
stays in a fixed positive bounded interval. Each fixed relative-prime
error bound applies to EVERY sufficiently large real threshold.
The cap margin proves Q,h^d guards. For the advice use the EXACT bound
n^(1/d)<=2^(L/d)<=K^(1-beta), hence h<=3d*K^(1-beta).
Also d<=4(L+1) and every accepted s<=(L+1)^2. Consequently

 400*m*h*d*s^2<=19200*K/(L+1)^6=o(K).

The remaining ceil(log2(32m))+3<=L+8 is o(K) in this window.
Thus the advice guard eventually passes uniformly, and m>=1 because
beta*W-12log2(L+1) tends uniformly to infinity. No fixed exponent
slack or unjustified neglect of h is used. q>C0*s gives h<q-1
and m<q<=Q. This proves every guard eventually, with one threshold
for all allowed contexts in this entropy window. For the tree, clipping
Lpad=ceil(Wmax^(20/19)) and Wmin^3>=(globalL+1)^2 put each actual
W_j in that window for its component ambient; the current proof gives
the detailed clipping application.

In particular for L0=ceil(log2 n), the actual FULL seed satisfies
(d+1)s<seedguard*L0 eventually: L=L0+O(1) and the coefficient gap is
exactly epsilon0. At every accepted finite tuple,
(d+1)s<=2ds<=2cap*log2(n+1)<=12log2(n+1).
The degree search visits at most(L+1)^2 candidates. All integer tests
and rejected paths are polynomial-time. Context fractions have fixed-degree
rational expressions in pRate,qRate, so their numerators/denominators
and the bit lengths of all pre-root powers are polynomial in L. The
pre-rejection q^d has O(L^3) bits. The accepted Q,h^d guards bound every
subsequent enumeration by a fixed power of n, independently of context.

For completeness, the actual field algorithms need no oracle. Construct
F_(2^s) by the supplied prime-degree binary constructor. For its degree-d
extension enumerate monic polynomials in coefficient order, testing
Z^(q^d)=Z modulo f and gcd(Z^(q^(d/p))-Z,f)=1 for each prime p|d.
For d=1 use f=Z. The finite-field root/orbit facts proved in cubic.md
show that the first test forces distinct irreducible factors of degrees
dividing d; every proper divisor divides d/p for some prime p|d and is
excluded by the gcd tests. Conversely an irreducible degree-d polynomial
passes. Such a polynomial exists: the union of proper subfields in F_(q^d)
has at most sum_(e=1)^floor(d/2)q^e<q^d elements; an element outside
has orbit length d and its minimal polynomial is suitable. Existence
of the finite fields themselves follows from the root-set construction
in cubic.md. Trial division supplies the prime divisors.

Factor Q-1 by trial division and choose the first nonzero g satisfying
g^((Q-1)/p)!=1 for every prime p|Q-1. The elementary cyclicity proof
in cubic.md Appendix A.2 proves existence and this exact order criterion.
Multiplication by g supplies the transitive nonzero Singer cycle A.
Searching all q^d polynomials and Q elements, modular gcd/powers and
trial division are polynomial in Q,d,s; Q<=(n+1)^6 bounds them uniformly.
Monomial enumeration and evaluating the m trace-polynomials at successive
A^i v are polynomial as well. No list decoder, predictor table or
long information-theoretic curve is enumerated by the evaluator.

# Small-error linear extraction with logarithmic output trimming

For input length n and integer entropy parameter 1<=r<=n, let
L=ceil_log2(n+1), m=floor(r/(L+1)^12), and
U=floor(beta*floor(log2 r)/2)-24ceil_log2(L+1).
Take p the largest odd prime at most U, Q=2^p, z=floor(Q/(L+1)^2).
Reject if no such prime exists or m<1 or z<2. Use the preceding local
code at (n,m,z), with code field exponent s and index length ell=2s.
Let A=max(1,ceil(log_Q m)), d=Q^2. Require ell<=Q, A<=ell and

 2[d+m*2^(A-1)+s+z+ceil_log2(4m)]<=r.                 (NWG)

These are the exact constructor guards. There is no additional ambient
window restriction for this finite theorem. Use the first ell elements
of F_Q as evaluation positions. The ith base-Q coefficient polynomial
of degree<A supplies a set {(x,P_i(x)):x among those positions}.
Every set has ell points; different polynomials intersect in at most A-1.
The reconstruction lemma's advice is at most m*2^(A-1). Thus (NWG)
proves strong error epsilon=2^-z for entropy r/2. Fixed-seed output is
linear. No field table is enumerated: both binary fields and code bits
are evaluated by their supplied deterministic algorithms.

For an affine source of entropy at least r/2+M, M=floor(r/4),
iterate J=ceil(M/m) copies with independent seeds and trim the last.
Conditioning on previous linear outputs costs at most M bits. Every
possible history therefore retains entropy r/2. The hybrid argument
gives near-output error eta=J*2^-z and full seed length J*d. The
extra M entropy requirement is essential and is charged below.

Here are UNIFORM eventual estimates for beta in the entire public
interval and L^(2/3)<=W=log2 r<=L^(19/20). The supplied relative
prime-interval theorem implies p/U tends to1: apply it to
[U/(1+delta),U] for any fixed delta>0. No quantitative prime-gap bound
is needed. Since log L=o(W), log2 Q=(beta/2+o(1))*W and
log2 z=(beta/2+o(1))*W, uniformly. Also z dominates L.
The code field has s=O(L+z+log m)=O(z), by the supplied prime
existence bound; hence ell=O(z)<=Q eventually. A stays bounded by
one fixed constant uniformly. The exact upper bound
Q^2<=r^beta/(L+1)^48 follows from p<=U. Consequently d=o(r),
m*2^(A-1)=O(r/(L+1)^12)=o(r), and every (NWG) cost is o(r).
All guards eventually succeed. Moreover J=O((L+1)^12), so
J*d=O(r^beta/(L+1)^36). These estimates also prove deterministic
polynomial runtime: p=O(log r), Q<=r^(.265), and the numeric code
field exponent s is polynomial in n. Prime searches by trial division
and the supplied field constructor therefore have fixed polynomial cost.
The exponentially large code field itself is never listed.

# Augmented affine row map

Put R=floor(floor(K^beta)/(L+1)^12), g=floor(K^(1-beta)),
w=floor(floor(K^beta)/(L+1)^24), M=floor(K/4).
For 1<=t<=g and 0<=take<=min(R,M), retain the exact row evaluator:
pair consecutive w-bit prefixes, duplicating the last when odd;
apply the half-entropy cubic extractor to each pair and trim its
output to the full near-extractor seed at ORIGINAL input length n
and parameter K; apply that near extractor to the original raw X.
Each new row has M bits. Repeat until one row remains and return
its requested prefix. If t=1, return its prefix directly.

With t_j=ceil(t/2^j), Jround=ceil(log2 t), total exposed slice bits
are at most (2t+Jround)*w. The exact constructor requires this at
most floor(K/4), w<=min(R,M), take<=min(R,M), and all near/cubic
guards. Eventually R/w is asymptotic to (L+1)^12, while
(2g+ceil_log2 g)*w<=3K/(L+1)^24=o(K). The near full seed is
O(K^beta/(L+1)^36); cubic output is Theta(w), larger by a growing
logarithmic factor. Thus every guard succeeds uniformly in the window.
For any fixed delta>0, the eventual error is at most
2^(-K^(beta/2-delta)); the finite bound is the explicit mixture
error below, not a claimed exact pure-power exponent beta/2.

## Conditional affine invariant and seed independence

Maintain a mixture of joint affine distributions of(X,current rows) in which the raw X has entropy at least K minus the total number of previously exposed slice bits, current rows are affine functions of X, and at least one row is uniform at its full current length. Initially this is exactly the hypothesis. The mixture approximation incurs the stated error; it is not an assumption that the unconditional nonlinear algorithm always has affine output.

Fix one such mixture component at a round. Let Z collect all t_j*w current slice bits. These are affine linear observations of X. Choose a basis for the kernel of the observation map on the direction of X and extend it. This gives X=A(Z)+B with B an affine source independent of Z. On every possible fiber its entropy is at least

 K-(sum_(h<j) t_h*w)-t_j*w >=K-floor(K/4) >=K/2+M.

Consequently the supplied near-linear-output strong extractor applies to B at PARAMETER K. This is the point of retaining raw X: its entropy is not limited to the current good row length, and the output is always M rather than one fourth of the previous row length.

A pair containing the current uniform row has slice-pair entropy at least w. Its cubic seed S is etaCubic-close to uniform, etaCubic=2*2^(-w/2000). Since S is a function of Z, it is independent of B before conditioning on Z. For every fixed seed s, the output on B is affine, so its variation distance from uniform is either zero or at least1/2. Strong extraction with a uniform seed and error eta implies that at most2eta of seeds are bad. Under S, that probability is at most2eta+etaCubic. A translation A(Z) changes only the fixed affine offset of the output, so this same bad-seed set applies to the actual X after conditioning on Z.

For every good Z value, all seeds and A(Z) are fixed. Every next row is affine in B, the specified good pair gives a uniform M-bit row, and the raw X still has the displayed entropy. Thus the conditional joint distribution satisfies the invariant for the next round, with the additional slice-bit budget charged. On exceptional Z values, replace the joint distribution by any affine distribution with the same affine X marginal and one uniform M-bit row given as an affine projection of X. Such a projection exists because H(X)>=K/2+M>M. Other rows can be arbitrary affine functions. This replacement preserves the invariant and changes total variation by at most2eta+etaCubic. Future computation is a fixed deterministic function of the joint pair(X,rows), so data processing and convexity justify iterating the bound even on replacement components. No independence from the whole raw X is claimed for S.

After Jround=O(log K) rounds the total error is at most Jround*(2eta+etaCubic). Since log2 z=(beta/2+o(1))*log2 K and w=K^(beta+o(1)), it is eventually at most2^(-K^(beta/2-delta)) for every fixed delta>0. A final prefix cannot increase variation. The t=1 case is exact. This proves the augmented row lemma with all history conditioning and raw-entropy losses charged. It uses no external somewhere-random extractor theorem.


# A block-source extractor using the augmented row map

Let the joint source(X1,...,Xb) be affine, each block n bits, and H(Xi|X1,...,X_(i-1))>=K on every possible full-history fiber. Require L^(2/3)<=log K<=L^(19/20), L=ceil(log2(n+1)); the numerical entropy guard of the small-error helper holds eventually in this window. Use the supplied constant-error strong linear extractor on X1 and enumerate every full seed. Use its complete R=floor(floor(K^(bSU))/(L+1)^12) bits, forming a v-by-R affine matrix W1, v=2^d0. At least one row is uniform because each fixed-seed affine output is either uniform or at distance at least1/2, while the average strong error is at most1/4.

Let g=floor(K^((1-bSU))), eventually at least2. Let E be the SMALLER-output small-error extractor at parameter K and input length n, trimmed to R output bits. Its initial output m0=floor(K/(L+1)^12) is at least R exactly; its full seed length is d1=Q^2<=K^beta/(L+1)^48, and strong error is epsilon1=2^-z from the actual small_parameters tuple. It needs only entropy K/2. Do not substitute the near-linear-output routine for this transfer: that longer seed is used only INSIDE the augmented row map.

To turn Wi into Wi+1, partition its rows into consecutive groups of g, padding the last group with zero R-bit rows. For each group compute the augmented row map on the PAIR(Xi,group), trimming to d1 bits. This trim is available eventually because d1=o(R) and d1=o(M), M=floor(K/4). Use each resulting seed in E on the NEXT block Xi+1, giving the next R-bit matrix rows. The raw block supplied to the group map is Xi, not Xi+1. Code retains it as previous while advancing to the next block.

Choose b=1+min{j:g^j>=v}. The row-count recurrence gives one final row after b-1 transfers. Since log g=Theta(log K) and d0=O(log n), b=O(log n/log K)=O((log n)^(1/3)) throughout the entropy window. Output that final R-bit row or any shorter requested prefix.

## Full-history conditional invariant

There is a set of full histories(X1,...,X_(i-1)) of probability at least1-(i-1)*tau such that, conditioned on each such history, Wi is an affine function of Xi and has one uniform row. Here tau=2epsilon1+epsilon2, epsilon2=ceil_log2(g)*(2eta+etaCubic). Initially this holds exactly.

Fix a good history. The raw Xi has entropy at least K by the block-source hypothesis. Every current matrix row is an affine function of Xi because fixing the previous input blocks fixes all previous transfer seeds. At least one group has a uniform R-bit row. The augmented row lemma therefore shows that THIS group's seed, a function of Xi, is epsilon2-close to uniform. All other group outputs may be correlated and arbitrary functions of Xi.

The joint pair(Xi,Xi+1) on this history is affine. Decompose Xi+1=A(Xi)+B by a kernel/complement split, where B is affine, independent of Xi and H(B)=H(Xi+1|Xi,history)>=K. The good group's seed is independent of B, being a function of Xi. Strong linear extraction on B implies that at most2epsilon1 of uniform seeds fail to yield a uniform R-bit output. Thus its actual seed is bad with probability at most tau. For each Xi value where it is good, conditioning on Xi fixes all group seeds and A(Xi); every next row is affine in B, the specified next row is uniform, and the full block-source hypothesis remains available at the next step. This proves the invariant after discarding probability at most tau from each previously good history. No group output is assumed to reveal the whole Xi; full original blocks are conditioned.

The final output is within(b-1)*tau of uniform. Since log2 z=(beta/2+o(1))*log2 K and b and the history factors have polynomial size in n, for every fixed delta>0 the final error is at most2^(-K^(beta/2-delta)) eventually in the window. More directly, along the global schedule log2 z_j=(a/2+o(1))*Z, while all multiplicative history factors have logarithm polynomial in global log N, negligible compared with2^(cZ) for any fixed c>0. Any fixed p_conf<a/2 therefore gives EBmax<<2^(-cP) for every fixed c, where P=2^ceil(p_conf Z). Output length is R=floor(floor(K^(bSU))/(L+1)^12).

## Uniform runtime and scope

Initial seed enumeration is polynomial in n by the supplied SU field/seed guard. Each augmented group uses O(log K) rounds, explicitly supplied cubic field construction and locally evaluated linear extractors; its raw input length is n. The group count decreases by g up to a ceiling on every transfer. Numerical ranks and guards depend only on n,K and row counts, not on target bit values. Every step therefore has fixed polynomial runtime and is total under the stated finite fallbacks. Padding a short subblock to ambient n bits preserves affinity and entropy and never adds ignored target inputs.

# Response maps specialized to their actual public output

Default optimal_parameters/constant_parameters use
m=floor(floor(K^beta)/(L+1)^12), L=ceil_log2(n+1).
Optional public take requires1<=take<=this m, then sets m=take
BEFORE field degree selection. It retains d,h,context and ALL actual
finite SU guards. The generic proof uses actual m, not m=K^beta, so
its fixed-seed linearity/error1/4 apply. The wrapper below first retunes
K; BLOCK defaults are unchanged. All cache keys include public take.

For a response of requested m bits at ordinary tuning K and context
beta=p/q, response_parameters instead sets
Kresp=max(2,ceil_root((m*(L+1)^12)^q,p)), requires Kresp<=K, then calls
optimal_parameters(n,Kresp,mode,take=m). Thus floor(floor(Kresp^beta)/(L+1)^12)>=m
exactly. Its returned map uses its ACTUAL new d,h,field and complete
output; it is not a prefix of a map tuned at ordinary K. Every finite
SU guard is checked directly and T(c)<=Kresp<=K. Invalid cases reject.

For each fixed feasible tuple, log m=(p_conf+o(1))Z since m=r_jP and
log r_j=o(Z). Hence log Kresp=(p_conf/beta_j+o(1))Z. It lies in the
same uniform entropy window eventually: it is a fixed positive multiple
of Z>>L^(2/3), while Kresp<=K gives the upper window bound. The mode
context is unchanged. Also m=(1+o(1))*Kresp^beta_j/(L_j+1)^12. The supplied uniform
SU analysis therefore gives the SAME seedguard(beta_j)*log n bound,
not an unjustified monotonicity comparison at the old K. The response
entropy obeys log T(c)<=log Kresp<=((p_conf*w0/a)+o(1))Z. Since
p_conf<a/2 and wm>(w0+X)/2, this exponent is STRICTLY below wm.
All constructor/root powers and seed enumeration retain their checked
fixed polynomial bounds. Default BLOCK maps retain ordinary K.

# Finite entropy thresholds for the SAME tuned maps

The public tuning parameter K is distinct from a map's required source
entropy. Keep every constructor, field, output bit and error unchanged.
For any accepted SU tuple c=(m,s,d,h,seedlen) define

 T(c)=400*m*h*d*s²+ceil_log2(32m)+3.

The supplied generic SU reconstruction proves strong error1/4 at source
min-entropy T(c): its other hypotheses concern only the fixed map.
Consequently its fixed-seed linear maps have a uniform affine output at
some seed whenever the affine input has entropy>=T(c). The constructor
already checks T(c)<=K, including for the specialized response output.

For an accepted small tuple p=(m0,z,s,ell,pdeg,Q,A,D), define, for1<=r<=m0,

 P_p(r)=D+r*2^(A-1)+s+z+ceil_log2(4r).                (PREFIX)

The SAME small_extract map's first r bits are a strong extractor at
entropy P_p(r), with error eps=2^-z. Indeed, the code was built with
rho0=eps/(2m0). Its list bound at agreement1/2+rho0 also holds at the
larger agreement1/2+eps/(2r). Restrict the supplied reconstruction to
its first r design sets: its overlap advice is<r*2^(A-1), list advice
is s, and the remaining seed/hybrid costs are exactly bounded by(PREFIX).
All code positions and sets remain unchanged. Fixed-seed linearity is
unchanged. Moreover P_p(r)<=P_p(m0)<=K/2 by the existing finite guard.

The near-linear map computes M=floor(K/4) bits in successive m0-sized
blocks. Its first r<=min(m0,M) bits are EXACTLY one small_extract prefix,
using the first D seed bits. Additional unused seed bits are independent
when the full seed is uniform, so the strong error remains eps with the
FULL near seed included. This is a prefix theorem, not a claim that the
entire M-bit output is uniform at this smaller entropy.

At a block-extractor context let c be its default SU tuple, p its small
tuple, R its default row length, g its grouping size, w the reservoir's
slice width, and take=D its transfer seed length. Set

 e=max(w,D), S=(2g+ceil_log2 g)*w,
 Hsrc=max(T(c),S+P_p(e),P_p(R)).                      (BLOCK)

The existing guards give e<=min(R,M),R<=m0 and S<=K/4. Thus e fits in
the FIRST near block; Hsrc<=K. Also P_p(e)>=e and T(c)>=R.
The SAME block extractor works for affine blocks of conditional entropy
>=Hsrc with its original finite error bound.

Here is the full modification to the supplied row proof. Initially SU
at T(c) gives one uniform R-bit row. In reservoir rounds retain one
uniform e-bit PREFIX, rather than a uniform full M-bit row. All complete
rows remain affine functions of the raw source for fixed exposed slices.
At most S slice bits are exposed in total; the independent raw kernel
therefore retains entropy>=Hsrc-S>=P_p(e)>=e. A pair containing the good
w-bit slice still has joint affine entropy>=w, so the same cubic map's
FULL seed is etaCubic-close to uniform and independent of that kernel.
By(PREFIX) the near map's e-prefix has strong error eps; at most2eps of
uniform full seeds fail to give a uniform affine e-prefix. Translation
by the exposed slice values does not change this bad-seed set. Hence
the good pair fails with probability<=2eps+etaCubic. On bad slice fibers
replace rows by affine rows with one uniform e-prefix and the SAME raw
marginal; such a projection exists since raw rank>=e. This changes joint
total variation by at most that bad probability and preserves the next
round's invariant. Thus the supplied mixture/data-processing induction
applies, and the final take<=e seed has the old conservative group error:
its old near error Jnear*eps is>=eps. The one-row case is exact.

Finally(PREFIX) at r=R replaces the transfer extractor's old K/2 source
requirement by P_p(R). The supplied full-original-history block induction
now uses entropy Hsrc for each block. Its independent next-block kernel
has that entropy; the good-group seed is independent of it; except the
same conservative bad-seed probability the next row is uniform R bits.
Every hypothesis and approximation step of that induction is retained.
The original EB, path and output error inequalities therefore suffice.
