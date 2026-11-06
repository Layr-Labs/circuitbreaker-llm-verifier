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

This lemma follows the interleaved-curve reconstruction of Shaltiel--Umans (2005), with larger constants, an explicit elementary list-size proof, and the seed allowed as an argument to the predictor. That last detail establishes the strong version needed here. This is the constant-error linear strong extractor used by family.py.

All logarithms below have base2. Let n,K,m,d,h,s be positive integers, q=2^s, Q=q^d, rho=1/(32m), and assume

    h>=2, h<q-1, 2m<2Q,
    q>2^30*h*d*s/rho^4,
    K>=400*m*h*d*s^2+ceil(log(1/rho))+3.

Let F be the explicitly constructed field with q elements. Encode x in {0,1}^n injectively and F2-linearly as a d-variate polynomial xhat over F, of total degree less than h. For example assign the bits as coefficients0 or1 of the first n monomials in increasing total degree and then lexicographic exponent order, provided binom(h+d-1,d)>=n. Let A be multiplication by a primitive element of the degree-d extension of F, on its standard d-dimensional F basis. Thus A is invertible and its powers act transitively on F^d minus0. The finite field, extension, primitive element and monomial choices are deterministic and specified below.

Define E_q(x,v)_i=xhat(A^i v), for v in F^d and1<=i<=m. Define the binary map E(x;(v,b))_i=Tr_F/F2(b*E_q(x,v)_i), for independent uniform b in F. It has seed length (d+1)s and output length m. We prove that it is a strong extractor for min-entropy K and error1/4, and that it is linear in x for each fixed seed. Neither proof predictors nor their decoders are computed by its evaluator.

## Seed-aware q-ary prediction is impossible

For every1<=i<=m and every deterministic function f taking (v,z1,...,z_(i-1)) to a list of at most rho^-2 field elements, we claim

    Pr_(x~X,v uniform)[E_q(x,v)_i in f(v,E_q(x,v)_<i)]<=rho

whenever Hinf(X)>=K. Lists and predictors can be computationally unbounded. Suppose contrariwise that this probability exceeds rho. Call a message x good if its success probability over v is at least rho/2. Good messages have X-probability greater than rho/2: if the mass were at most rho/2, the total success would be at most rho/2+(1-rho/2)*rho/2<rho.

Fix a good x and put r=100ds and B=(2r-1)(h-1), so B+1<=2rh. The field-size condition ensures2r<q. Choose2r distinct uniform parameters t1,...,t_(2r) in F, and independent uniform y1,...,y_(2r) in F^d. Interpolate vector polynomials p1,p2 of degree less than2r, with

    p1(tj)=yj for every j;
    p2(tj)=A*yj for j<=r, and p2(tj)=yj for j>r.

For0<=j<Q put P_(2j+1)=A^j p1 and P_(2j+2)=A^j p2. Each P_a is a uniformly random vector polynomial of degree less than2r, independent of the entire tuple of parameters. Indeed, for each fixed distinct parameter tuple, interpolation is a bijective linear map from the uniform vector values to the coefficients; the same is true after applying A to some values or to all coefficients. Its values at any2r distinct arguments are independent uniform vectors, by the invertibility of the corresponding evaluation map.

Successive curves P_(a-1),P_a agree on r distinct uniform parameters S_a independent of P_a. If a is even take the last r parameters; if a is odd and at least3 take the first r. The identities follow directly from the two interpolation conditions. Independence follows from the preceding coefficient argument. These assertions are marginal assertions for each a; no independence between different curves is required.

Give a reconstruction procedure the coefficients of xhat composed with the INITIAL2m curves. Each composition has degree at most B, so this advice has at most

    a0=2m*(B+1)*s<=4mhrs=400mhds^2

bits. Pad to the integer upper bound a=400mhds^2. Giving all2m startup curves avoids any boundary problem when the predictor index is i=1. Shared random curves are random coins of the reconstruction, not advice bits.

Suppose a>2m and all earlier compositions have been learned. For u=P_a(t), its previous q-ary output values at seed v=A^-i u are

    xhat(A^(j-i)u)=(xhat composed with P_(a-2(i-j)))(t), 1<=j<i.

All these curve indices are positive and earlier. Hence the predictor defines a list L_a(t) by taking as input that known prefix AND the known seed A^-i P_a(t). This last argument is why the proof works for a seed-aware predictor. With the true previously computed polynomials, L_a depends only on P_a, x and the fixed f, not on the overlap parameter set S_a.

For uniform u in F^d the chance of correct prediction is at least rho/2: A^-i preserves uniform measure and x is good. On a random P_a, therefore, the q indicators of correct prediction are2r-wise independent Bernoulli variables with common mean p>=rho/2. Put mu=qp, which exceeds2r by the field guard.

Here is a self-contained tail estimate. For2r-wise independent Bernoulli variables Y_t of mean p, expand the2r-th centered moment of their sum. Any term in which some index appears exactly once vanishes. With j distinct indices,1<=j<=r, the absolute expected product is at most p^j, since for every e>=2, E|Y_t-p|^e<=p. Choose a representative position for each block and map the remaining positions to representatives. There are at most binom(2r,j)*j^(2r-j)<=4^r*(2r)^(2r-j) possible partitions, and at most q^j assignments of distinct field indices. Thus

    E[(sum Y_t-mu)^(2r)]
      <=sum_(j=1)^r 4^r*(2r)^(2r-j)*mu^j
      <=r*4^r*(2r)^r*mu^r <=(16r mu)^r.

The second inequality uses mu>=2r, and the last r<=2^r. Markov's inequality gives

    Pr[sum Y_t<mu/2]<=(64r/mu)^r
       <=(128r/(rho q))^r <2^(-16r)=Q^-1600.

Consequently, except for this probability, the true degree-B polynomial xhat composed with P_a agrees with the lists L_a(t) on at least rho*q/4 arguments.

To bound the possible such polynomials, let T=floor(rho*q/8). There are at most q/rho^2 input pairs(t,z), z in L_a(t). Form a nonzero bivariate polynomial R(Tvar,Y) of weighted degree at most T, assigning weights1 and B to its variables, which vanishes on all pairs. This is possible by homogeneous linear algebra: the rectangular monomial set with exponents up to floor(T/2), floor(T/(2B)) has more than T^2/(4B) coefficients, and, since T>=rho*q/16,

    T^2/(4B)>=rho^2*q^2/(1024B)>q/rho^2.

The final inequality follows from q>2^30 hds/rho^4 and B<200hds. If z(t) has degree at most B and agrees at rho*q/4 points, R(t,z(t)) has degree at most T and more than T roots, so vanishes identically. Thus Y-z(Tvar) divides R. Over F(Tvar), different such z give different linear factors. There are at most deg_Y R<=T/B<q candidates. This proves only a finite list-size bound; no decoding algorithm is needed.

The candidate list depends only on P_a and the fixed x,f. Test these candidates using the r known values at the overlap with P_(a-1). A false candidate differs from the true polynomial on all but at most B field points. Since S_a is a uniform distinct r-tuple independent of P_a, the probability that it passes every test is at most (B/q)^r: sampling without replacement gives a product of factors at most B/q. Taking a union over fewer than q candidates bounds the disambiguation failure by

    q*(B/q)^r < 2^s*(2^-22)^r <=2^(-2199ds).

Union bounds over at most2Q steps show that all learning steps succeed with probability at least1-2Q*(Q^-1600+2^(-2199ds)), exceeding3/4 under the displayed guards. This is an unconditional union bound on the ideal learning steps; it does not condition the uniform-curve assertions on earlier success.

Also p1(0) is uniform in F^d and is nonzero except with probability1/Q<1/4 (the field guard makes Q much larger than4). If it is nonzero, its A orbit, occurring among the odd curves, covers every nonzero vector. Hence with probability greater than1/2 reconstruction obtains xhat on all nonzero vectors. These evaluations uniquely determine a polynomial of total degree less than h<q-1. Indeed their difference vanishes at all nonzero points; on any line through0 there are q-1 roots and degree less than q-1, giving the value0 at the origin as well. A polynomial with individual degrees less than q vanishing on all of F^d is zero, by iterated univariate root bounds. The injective encoding then recovers x.

For each good message there is the specified a-bit advice, and reconstruction with that advice succeeds with probability at least1/2 over shared coins. Therefore the X-probability of messages reconstructible from SOME a-bit advice, averaged over those coins, exceeds rho/4. For any fixed coins the deterministic procedure has at most2^a possible outputs. Since Hinf(X)>=K, that probability is at most2^(a-K)<=rho/8 by K>=a+ceil(log(1/rho))+3. This contradiction proves the seed-aware q-ary prediction claim.

## Binary conversion with the seed retained

Suppose that (v,b,E(X;(v,b))) has total variation more than epsilon=1/4 from (v,b,U_m). The ordinary hybrid argument yields a fixed coordinate i and a deterministic next-bit predictor P(v,b,previous binary outputs) with advantage greater than epsilon/m. For completeness, telescope hybrids that successively replace output bits by independent uniform bits. A distinguishing test has adjacent gap greater than epsilon/m. Its difference is one half of the expected product of the actual signed bit and the test's difference of conditional acceptance probabilities after setting that bit to0 or1. Predict0 with probability(1+d0-d1)/2. This randomized predictor has that advantage; fix its finite table of coins by averaging to obtain a deterministic predictor.

For a seed v and a q-ary prefix z_<i, define a sign-valued function H_(v,z)(b) by running P on the trace prefix (Tr(b*z1),...,Tr(b*z_(i-1))). For c in F set

    gamma_(v,z)(c)=E_b H_(v,z)(b)*(-1)^Tr(bc).

The q trace characters are an orthonormal basis: for nonzero c, multiplication by c permutes F and the nonzero linear trace functional is balanced. Nondegeneracy of trace follows from the supplied prime-degree field construction, since s is odd and Tr(1)=1. Parseval gives sum_c gamma(c)^2=1. Thus the list of symbols with gamma(c)>=2rho has size at most1/(4rho^2)<=rho^-2, and it depends only on v and the FIELD prefix, as permitted in the preceding lemma.

The mean gamma at the actual next symbol, averaged over x,v, is twice the binary prediction advantage, strictly larger than2epsilon/m=16rho. If the chance theta that it is at least2rho were at most rho, this mean would be at most theta+(1-theta)*2rho<=3rho. Therefore theta>rho. That violates the seed-aware q-ary prediction claim. We obtain strong error at most1/4.

Every output is an evaluation linear in the message coefficients, followed by multiplication by the fixed b and binary trace. Thus fixing the complete seed fixes a binary linear map on x.

## Explicit polynomial-time choices in the subpolynomial-entropy regime

Let L=ceil(log(n+1)), w=floor(log K), m=floor(K^(333/500)), d=ceil(3L/w), and h=d*ceil_root(n,d)+1. Use these only when w>=1. The number of monomials of total degree less than h is binom(h+d-1,d)>=((h-1)/d)^d>=n, so the stated bit-coefficient encoding is available. Choose s as the least odd prime not smaller than the bit length of

    2^40*h*d*(L+1)^2*(32m)^4,

and q=2^s. Require the finite guards

    s<=(L+1)^2, q>2^30*h*d*s*(32m)^4,
    q^d<=(n+1)^10, h^d<=(n+1)^8,
    h<q-1, m<Q,
    K>=400mhds^2+ceil(log(32m))+3.

If any guard fails this component returns a specified fallback and is not used in the eventual theorem. Guards use only integer arithmetic and roots. For L^(2/3)<=log K<=L^(19/20) every guard holds eventually. Indeed w=Theta(log K), d=Theta(L/log K), d*log L=o(L), and log h=L/d+O(log L)=(1/3+o(1))*log K. The FULL supplied relative-prime-interval theorem, applied to the explicit bit-length threshold as detailed below, gives s=(1124/375+o(1))*log K uniformly over this entropy window. Consequently

 ds=(1124/125+o(1))*L<10L,
 d log h=L+O(d log L)<8L.

The earlier coarse factor2 prime interval is insufficient for the q^d guard at these longer outputs and is not used. The advice logarithm is

 log(400mhds^2)=(333/500)log K+log h+O(log L)
              <=(1499/1500+o(1))*log K,

strictly below log K eventually. The field-size guard holds because s=O(log K+log L)=O(L), while (L+1)^2/s tends to infinity; the defining threshold has the additional factor(L+1)^2. All other guards follow from these estimates. They are uniform since log K dominates log L, d tends to infinity uniformly, and d*log L=o(L). The tree uses log K=(148+o(1))*L^(3/4)*(log L)^(1/4), lying in this window eventually. The proof below of the full prime estimate uses only the supplied analytic prime theorem and the displayed integer threshold, not acceptance of any extractor guard.

Here are deterministic field and Singer-cycle constructions without oracles. The base field F_(2^s) is the supplied polynomial-time prime-degree binary field constructor. For the degree-d extension enumerate monic degree-d polynomials in coefficient order and use the following irreducibility criterion: Z^(q^d)=Z modulo f, and gcd(Z^(q^(d/p))-Z,f)=1 for every distinct prime divisor p of d. For d=1 just use f=Z. To justify the criterion, the roots of Z^(q^j)-Z in an algebraic closure form the field with q^j elements (they are distinct, closed under addition/multiplication/inversion, and there are q^j roots). Its irreducible factors have degrees dividing j: the Frobenius orbit length of a root is exactly its minimal-polynomial degree. The stated congruence forces all irreducible factors of f to be distinct and have degrees dividing d; every proper divisor of d divides d/p for some prime p|d and is excluded by the gcd tests. Conversely an irreducible degree-d polynomial satisfies the tests.

Such a polynomial exists: in F_(q^d), the union of all proper subfields has at most sum_(e=1)^floor(d/2) q^e<q^d elements for d>=2, since every proper subfield degree is a proper divisor of d. An element outside has Frobenius orbit length d. Its minimal polynomial supplies an irreducible polynomial. The case d=1 is immediate. Enumeration of all q^d possibilities and modular powering/gcd takes polynomial time in Q,d,s. Prime divisors of d are obtained by trial division.

Factor Q-1 by trial division. In the constructed extension find the first nonzero g with g^((Q-1)/p)!=1 for every distinct prime p dividing Q-1. Such an element exists: for each prime power p^e dividing Q-1 the polynomial Z^((Q-1)/p)-1 has fewer than Q-1 roots, so choose a nonroot a. By Lagrange's theorem the order of a has full p-adic valuation e; raising a to (Q-1)/p^e gives an element of order p^e. The product of these elements for the different primes has order the product Q-1, since they commute and their orders are pairwise coprime. This proves cyclicity and the search criterion without assuming it. Exhaustive search over Q elements and trial division up to sqrt(Q-1) are polynomial time in Q. Multiplication by the first qualifying g gives the required A and its transitive nonzero orbit.

The guards Q<=(n+1)^10 and h^d<=(n+1)^8 bound these searches and enumeration of all monomials by fixed polynomial powers of n. Evaluating m polynomials at successive A^i v, followed by m traces, is polynomial time. The numerical guards give seed length (d+1)s<=2ds<=20log(n+1) for every accepted parameter tuple. A sharper EVENTUAL estimate is available from the complete supplied prime-interval argument and will be used by the tree. Put W=log2 K. Uniformly for L^(2/3)<=W<=L^(19/20), the exact d=ceil(3L/floor W) obeys d*W=(3+o(1))*L and d tends to infinity. Since ceil_root(n,d)=n^(1/d)*(1+o(1)) in this range, log2 h=L/d+log2 d+o(1)=(1/3+o(1))*W. Also log2 m=(333/500+o(1))*W, and log L=o(W). The integer threshold T whose next odd prime defines s is consequently

 T=log2(2^40*h*d*(L+1)^2*(32m)^4)+O(1)
   =(1124/375+o(1))*W.

The supplied theorem gives a prime in(X,(1+delta)X] for every fixed delta>0 at all sufficiently large X. Hence the least prime s>=T satisfies s/T->1: for each fixed delta the ratio is eventually at most1+delta+o(1), and it is always at least1. Oddness is automatic for large T. Thus s=(1124/375+o(1))*W and the FULL binary seed length satisfies

 (d+1)*s=(1124/125+o(1))*L.

All error estimates are uniform throughout this entropy window, since d grows at least as a positive power of L and W/log L tends to infinity uniformly. For the tree's L0=ceil(log2 n), the distinction L=L0 or L0+1 only contributes a relative o(1). Therefore the seed length is eventually at most(8993/1000)*L0, since1124/125<8993/1000; in particular it is at most9L0. This sharper estimate uses the full supplied relative-prime-interval proof, not a coarse factor2 interval or an external prime oracle. No reconstruction or received-word enumeration is performed by the evaluator. The proof's very long curve sequence and predictor tables are information-theoretic counting devices only.


# Small-error linear extraction with a polynomial seed budget

Let a source consist of N bits, with entropy at least r/2+floor(r/4), where1<=r<=N and ceil(log2(N+1))<=floor(r^(1/2000)). For all sufficiently large parameters, construct a linear strong extractor whose output length is M=floor(r/4), error at most2^(-r^(771/10000)/3), and seed length O(r^(1021/2500)). All source hypotheses here are AFFINE when output iteration is used. Integer powers in the construction mean exact floor roots of integer powers; floating-point arithmetic is unnecessary.

First make a smaller-output extractor. Put z=floor(r^(771/10000)), epsilon=2^-z, m=floor(r^(9/10)), rho=epsilon/(2m). Use the locally evaluated code with index length ell=2s and list bound2^s from the corresponding section of extractors.md. Then ell=O(log N+z+log m)=O(r^(771/10000)), and ell=Omega(r^(771/10000)), because log N<=r^(1/2000)=o(r^(771/10000)).

Choose a prime degree p>=ceil(log2 ell), and construct F_Q with Q=2^p. Choose the least such prime; the supplied prime-interval proof ensures Q<=ell^2 eventually, after a harmless fixed factor if the lower endpoint is rounded. One may use Q<=4ell^2 throughout to absorb that factor. Let A=ceil(log_Q m), meaning the least integer A with Q^A>=m. Eventually A is bounded by a FIXED constant, since log Q>=log ell=Omega(log r). Also A<=ell. Encode integers0,...,m-1 as distinct coefficient vectors of polynomials P_i in F_Q[T] of degree at most A-1.

For the first ell field elements alpha in their binary integer order define the ordered set

 S_i={(alpha,P_i(alpha)):alpha among those ell elements}.

The universe is F_Q^2, of size d=Q^2=O(ell^4). Each set has EXACTLY ell elements, not Q padded elements. Distinct polynomials agree at at most A-1 field points, by the root bound. Hence the reconstruction advice B is at most m*2^(A-1)=O(r^(9/10)), with a fixed, possibly large constant. The other costs d, s and log2(4m/epsilon) are O(r^(771/2500)). Their sum with B is o(r), therefore below r/2 eventually. The reconstruction lemma proves strong error at most epsilon. For fixed seed this extractor is linear in the source input.

Produce M bits by J=ceil(M/m) successive applications with independent fresh seeds, trimming the final output block. At each stage, condition on all preceding seeds and output blocks. For an affine initial source, this is either an impossible history or a uniform affine source whose entropy dropped by no more than the number of preceding output bits: every map is linear once the previous seeds are fixed. Thus its entropy stays at least r/2 as long as the initial entropy was at least r/2+M. The applications below verify this extra condition from their explicitly charged kernel entropy bounds. More generally the iterated guarantee needs initial entropy at least r/2+M; do not erase this extra requirement.

For that initial entropy, the strong-extractor conclusion applies to every possible history. Replacing the successive output blocks by fresh independent uniform blocks, one at a time, incurs total variation at most J*epsilon by the triangle inequality and conditioning. Therefore the concatenated output is a strong extractor with error at most J*2^-z<=2^(-r^(771/10000)/3) eventually. Its total seed length is J*d=O(r^(1021/2500)), and it is linear for every fixed complete seed. Every design and code bit is explicitly locally evaluated in time polynomial in N and r. Huge numbers q in the code specify index spaces only; their full tables are never enumerated.



# Affine row condensation with a grouping/precision tradeoff

Let X be an affine distribution of t-by-r binary matrices having at least one uniform r-bit row. For1<=t<=r^(91/1000), the deterministic construction below outputs floor(r^(817/1000)) bits within variation2^(-r^(63/1000)) of uniform, at every sufficiently large r. All cubic and seeded ingredients are supplied above.

## One step on an internal matrix

For2<=t<=r^(23/250), put w=floor(r^(817/1000)), M=floor(r/4). Slice the first w bits from every row, pair slices in order, duplicating the last when t is odd. Apply the supplied half-entropy cubic extractor to each pair and trim to the full seed length of the iterated small-error extractor on the full tr-bit matrix, at parameter r and output M. Apply that seeded map to the full matrix, one seed per pair. The next matrix has ceil(t/2) rows of M bits.

The full near-linear-output seed has length O(r^(1021/2500)); cubic output has Theta(r^(817/2000)) bits. Since1021/2500<817/2000, trimming succeeds eventually. The full matrix length tr<=r^(1023/250) satisfies the explicit small-error entropy guard eventually.

Let Z collect all tw slices. Split the affine source as X=A(Z)+B, with B affine and independent of Z, by choosing a basis for the kernel of the slice map and extending it. Since the good uniform row gives entropy at least r and the slice map has rank at most tw,

 H(B)>=r-tw>=r-r^(909/1000)>r/2+M

eventually. A pair containing the uniform row has slice-pair entropy at least w. Its cubic seed is beta-close to uniform for beta=2*2^(-sqrt(w)/50), independently of B.

For each fixed seed the iterated linear extractor output on B is affine; unless uniform its variation is at least1/2. Strong extraction with uniform seed and error eta<=2^(-r^(771/10000)/3) implies at most2eta bad seeds. A beta-close seed is bad with probability at most2eta+beta. This uses seed independence from B before conditioning on Z. After fixing Z, all seeds and A(Z) are fixed, so every next row is affine in B. Except on probability2eta+beta, the specified pair gives a uniform row. Replacing exceptional conditional distributions by arbitrary affine somewhere-random distributions proves closeness to a mixture of such sources, retaining arbitrary row correlations. This is the invariant for iteration.

## Iterate and retain every row length

Write t0=t,r0=r, gamma=91/1000. The exact recurrences give t_j=ceil(t0/2^j),r_j=floor(r0/4^j),J=ceil(log2 t0). Since4^J<4t0^2 and r0/t0^2>=r0^(409/500) tends to infinity, for every0<=j<=J,

 r_j>=r0/(2*4^j), r_J>=r0/(8t0^2)>=r0^(409/500)/8

eventually. For j<J,2^j<t0 and

 t_j<=2t0/2^j<=2^(1+gamma)*r_j^gamma*2^(-(1-2gamma)*j)
      <3r_j^(91/1000)<=r_j^(23/250)

eventually. The minimum r_j tends to infinity, making the last inequality uniform over every round, including odd row counts. Code checks the exact guards t0^1000<=r0^91 and t_j^250<=r_j^23. Finite failures return the total fallback.

The final row contains the requested prefix floor(r0^(817/1000)) eventually, because817/1000<409/500. Each deterministic round contracts variation and adds at most its uniform one-step error on each mixture component. The worst seeded error has exponent at least c*r0^(315339/5000000) for a fixed c>0, since(771/10000)*(409/500)=315339/5000000. This is strictly larger than63/1000=315000/5000000. The worst cubic error exponent is(409/500)*(817/2000)=334153/1000000, also larger. Summing O(log r0) rounds therefore gives variation at most2^(-r0^(63/1000)) eventually. If t=1, return the fixed prefix exactly.

All component calls have polynomial time in their own matrix lengths because t_j<=r_j^(23/250), fields are explicitly constructed and code bits locally evaluated. There are O(log r0) rounds; no full exponential codeword or proof decoder is enumerated.


# A block-source extractor with supplied subpolynomial-entropy components

This component follows Shaltiel's affine block-source reduction. Let each input block have n bits and let L=ceil(log2(n+1)) and require L^(2/3)<=log K<=L^(19/20). The joint source (X1,...,Xb) is affine and has conditional entropy at least K in each block after fixing any possible values of its predecessors. Conditional entropy is constant over possible affine fibers, so this is a pointwise condition, not an average entropy assumption.

Use the constant-error linear extractor on one block at entropy K and error1/4. Trim its output to R=floor(K^(333/500)) bits; its seed has d0=O(log n) bits. Enumerate all v=2^d0 seeds, obtaining a v-by-R matrix W1 of extractor values on X1. This is an affine matrix source. At least one row is uniform: for fixed seed the output is affine, nonuniform affine output has distance at least1/2, and average distance at most1/4 precludes all seeds being bad.

Let g=floor(R^(91/1000)), eventually at least2. Let E be the SMALLER-output small-error linear extractor at input length n and entropy parameter K, trimmed to R output bits. Its original output floor(K^(9/10)) is longer than R, its seed length is d1=O(K^(771/2500)), and its strong error is epsilon1=2^-floor(K^(771/10000)). It needs entropy only K/2, hence applies to the independent parts below. Do not substitute the iterated near-linear-output routine here; the short seed bound is for the single small-output call.

To turn Wi with u rows into Wi+1, partition its rows into ceil(u/g) consecutive groups, padding the last group with zero rows to g if necessary. Apply the corresponding section of extractors.md to each g-by-R group, trimming its output to d1 seed bits. This is possible eventually: its output length is floor(R^(817/1000))=K^(272061/500000) up to floors, much larger than d1=O(K^(771/2500)). Apply E to the NEXT input block using each resulting seed, and use the R-bit results as the new matrix rows.

Choose b=1+min{j:g^j>=v}. Then b=O(log v/log g)=O(log n/log K)=O((log n)^(1/3)); all quantities and the minimum are computable integers. After b-1 transfers, one row remains. Output it, or any specified shorter prefix.

## Conditional correctness invariant

For each i there is a set of histories (X1,...,X_(i-1)) of probability at least1-(i-1)*tau such that, conditioned on every such history, Wi is an affine somewhere-random matrix. Here tau=2epsilon1+epsilon2, where epsilon2<=2^(-R^(63/1000)) is the group extractor's error. For i=1 the history is empty and the initial construction established the invariant exactly.

Fix a good history before step i. Given that history, Wi is a linear function of Xi plus a fixed shift: all seeds used at previous steps become fixed when the earlier input blocks are fixed. The joint pair(Xi,X_(i+1)) is affine. Split it by linear algebra as X_(i+1)=A(Xi)+B, where B is an affine source independent of Xi and has entropy H(X_(i+1)|Xi,history)>=K. This follows by choosing a complement to the kernel of the projection onto Xi; it holds even if the blocks are correlated.

One group in Wi contains a uniform row, so its extracted seed is within epsilon2 of uniform. It is a function of Xi, independent of B. Strong linear extraction on B implies that at most2epsilon1 of uniform seeds yield a nonuniform R-bit output on B. Therefore that group's actual seed is bad with probability at most tau. For every Xi at which it is good, conditioning additionally on Xi fixes all next seeds and A(Xi). The new matrix is jointly affine as a function of B, with the specified group's row uniform. This establishes the invariant at step i+1 after discarding probability at most tau of each previously good history. The invariant concerns FULL previous input-block histories, so it does not assume that an extractor output reveals all of its source block.

At the last step, outside probability (b-1)*tau, the single output row is uniform conditioned on the full previous-block history. Its unconditioned distance from uniform is therefore at most(b-1)*tau. Since b=O((log n)^(1/3)), R=K^(333/500) up to floors, and epsilon2 has exponent K^(20979/500000) up to constant factors, the final error is at most2^(-K^(209/5000)) for all sufficiently large K. Output length is floor(K^(333/500)).

## Running time and scope

Enumerating2^d0 seeds is polynomial in n because the supplied optimal extractor has d0=O(log n) in this entropy range. Every transfer reduces the row count by a factor g up to a ceiling, uses a polynomial-time group extractor and local linear code evaluations, and there are b=O((log n)^(1/3)) transfers. The total evaluation and construction time is polynomial in n with a fixed exponent. The construction does not enumerate a codeword's exponentially large index space. Input blocks of length below n may be padded with zeros, preserving affinity and entropy.

This proves the stated block-source component under its conditional entropy hypothesis. The total tree construction, source restrictions and actual-length transfer are supplied in proof.md.