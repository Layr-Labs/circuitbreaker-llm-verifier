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

For the rest of this file choose ONE public rational rate beta=pRate/qRate with positive integers pRate,qRate,3/20<=beta<=53/100 and max(pRate,qRate)<=100ceil(log2(n+1)). This rational context is fixed throughout the invocation; every component/helper/row/block-source call propagates it explicitly,including all cache keys. It is never chosen from source bits or epsilon. The family computes its context from public W_j,Z. All constructors check the context-size guard BEFORE variable integer powers; invalid numerical parameters yield the stated fallback.

Put epsilon0=1/100000 and define EXACT rational constants

 bSU=beta, aSU=1/(1-beta)+epsilon0,
 seedguard=1+4aSU*beta+epsilon0,
 c=(9999/10000)(1-beta)/seedguard,
 zeta=beta/4-1/1000000,
 group exponent=1-beta, slice exponent=beta-1/1000000,
 sigma=(2000001/1000000)zeta, nu=sigma+1/1000000,
 eGroup=zeta-1/100000000, eBE=zeta-2/100000000.

The evaluator computes and reduces all fractions by integer gcd. A power floor(r^(u/v)) means exact root(r^u,v) for the reduced positive rational exponent. This is a uniformly specified rational-parameter algorithm, not a finite table or unproved passage from pointwise to uniform statements. We give explicit uniform margins and runtime bounds here.

For every allowed beta,(20/17+epsilon0)<=aSU<=(100/47+epsilon0) and .037499<=zeta<=.132499. Also

 1-(beta+1/aSU)=epsilon0(1-beta)^2/(1+epsilon0(1-beta))
                >=epsilon0*.2209/(1+.85epsilon0)>0.

Define cap(beta)=2 for beta<=19/100,3 for19/100<beta<=3/10,4 for3/10<beta<=2/5,5 for2/5<beta<=12/25,and6 otherwise within[3/20,53/100]. This is a pure mode function,so the existing mode cache key includes every dependency. The rational seedguard=1+4beta/(1-beta)+4epsilon0*beta+epsilon0 is increasing in beta. At the five bin endpoints19/100,3/10,2/5,12/25,53/100,exact arithmetic gives seedguard<cap(beta)-1/20. Thus every bin has a uniform fixed positive margin for BOTH ds and (d+1)s; the d log h=(1+o(1))L estimate has still more slack. The uniform prime-interval errors and integer rounding are absorbed by these margins. All accepted finite field searches obey the actual cap guard; the coarse power6 remains true. The defining degree directly solves the required finite predictor inequality:put C0=2^30*h*d*(32m)^4,start at the least odd prime>=bit_length(C0),and visit successive odd primes while2^s<=C0*s. Return fallback if s>(L+1)^2. The returned s is the LEAST odd prime satisfying2^s>C0*s within that guard. Every smaller admissible prime fails,since2^r/r strictly increases for integers r>=2 and the starting bound is below every possible successful degree. No predictor/reconstruction guard is weakened.

The full-seed coefficient is at most5.510660; its guard is at most5.510670<6,and the exact guard gap is epsilon0. Further c*seedguard/(1-beta)=.9999 exactly,and .08<c<.5 by the endpoint bounds on numerator and denominator. Group plus slice equals999999/1000000. Subtracting nu from the half-slice exponent gives(1/2-beta/4)/1000000+1/1000000000000>=.36/1000000+1/1000000000000>0. Finally eBE<eGroup<zeta,all bounded below by positive constants. These margins are uniform over the ENTIRE allowed interval; no finite-mode argument is used.

Let L=ceil(log(n+1)), w=floor(log K), m=floor(K^(bSU)), d=ceil(aSU*L/w), and h=d*ceil_root(n,d)+1. Use these only when w>=1. The number of monomials of total degree less than h is binom(h+d-1,d)>=((h-1)/d)^d>=n, so the stated bit-coefficient encoding is available. Choose s by this least-valid-prime procedure for C0=2^30*h*d*(32m)^4,and put q=2^s. Require the finite guards

    s<=(L+1)^2, q>2^30*h*d*s*(32m)^4,
    q^d<=(n+1)^cap(beta), h^d<=(n+1)^cap(beta),
    h<q-1, m<Q,
    K>=400mhds^2+ceil(log(32m))+3.

If any guard fails this component returns a specified fallback and is not used in the eventual theorem. Guards use only integer arithmetic and roots. For L^(2/3)<=log K<=L^(19/20) every guard holds eventually. Indeed w=Theta(log K), d=Theta(L/log K), d*log L=o(L), and log h=L/d+O(log L)=(1/aSU+o(1))*log K. The FULL supplied relative-prime-interval theorem, applied to the explicit bit-length threshold as detailed below, gives s=((1/aSU+4bSU)+o(1))*log K uniformly over this entropy window. Consequently

 ds=((1+4aSU*bSU)+o(1))*L<6L,
 d log h=L+O(d log L)<6L.

The q^d guard follows from the sharp estimate above; the full relative-prime-interval theorem also proves the tighter full-seed limit used by the tree. The advice logarithm is

 log(400mhds^2)=(bSU)log K+log h+O(log L)
              <=((bSU+1/aSU)+o(1))*log K,

strictly below log K eventually. The field-size guard holds because s=O(log K+log L)=O(L), while (L+1)^2/s tends to infinity; the conservative C0*(L+1)^2 threshold supplies an eventual upper comparison for the new least degree. All other guards follow from these estimates. They are uniform since log K dominates log L, d tends to infinity uniformly, and d*log L=o(L). For each selected public configuration,the tree uses W_j between Wmin=ceil(wm*Z) and Wmax=ceil(w0*Z). Global clipping Lpad=ceil(Wmax^(20/19)) and Wmin^3>=(globalL+1)^2 imply (Lc_j)^(2/3)<=W_j<=(Lc_j)^(19/20) at every level; proof.md supplies the complete uniform application. The proof below of the full prime estimate uses only the supplied analytic prime theorem and the displayed integer threshold, not acceptance of any extractor guard.

Uniformity of every asymptotic estimate above follows from these explicit bounds. The SU degree satisfies d=Theta(L/W) with coefficients in[20/17+epsilon0,100/47+epsilon0]; ceil/root relative errors are uniformly O(1/d)+O(d*log L/L). The prime threshold is(1/aSU+4beta+o(1))W,whose leading coefficient stays in a fixed positive bounded interval. The supplied relative-prime theorem applies to EVERY sufficiently large threshold; its ratio error is consequently uniform for all these thresholds,not just for fixed beta. Advice's displayed fixed positive gap and the fixed full-seed guard gap absorb all uniform o(1) terms. In the helper, zeta stays between.037499 and.132499,so design degree is bounded uniformly and all reconstruction overhead ratios tend to zero uniformly. Each fixed exponent separation used in row and block transfer has the uniform positive lower bounds just displayed. Thus a single eventual input threshold works for every public context satisfying the stated range and bit-size hypotheses.

Runtime also remains uniform. Since pRate,qRate<=100L,all reduced numerators and denominators of beta,zeta,slice and1-beta are O(L),with fixed constants. Computing r^u before its vth root therefore uses O(L^2) bits for r<=n. The rational alpha,seedguard,c have fixed-degree rational expressions in pRate,qRate,and their numerators/denominators have O(log L) bits. Forming d,all comparisons and rejected-path guards uses polynomially many bits and arithmetic operations. After s<=(L+1)^2,the pre-rejection Q=q^d has O(L^3) bits. The finite Q,h^d guards bound every accepted field search and seed enumeration by fixed powers of n,independently of beta. The subsequent locally evaluated helper/row/block operations have exponents bounded by the same compact interval and exact numerical guards. No context-dependent time exponent is hidden in a variable root power. proof.md verifies every actual tree context's pRate,qRate bound.

Here are deterministic field and Singer-cycle constructions without oracles. The base field F_(2^s) is the supplied polynomial-time prime-degree binary field constructor. For the degree-d extension enumerate monic degree-d polynomials in coefficient order and use the following irreducibility criterion: Z^(q^d)=Z modulo f, and gcd(Z^(q^(d/p))-Z,f)=1 for every distinct prime divisor p of d. For d=1 just use f=Z. To justify the criterion, the roots of Z^(q^j)-Z in an algebraic closure form the field with q^j elements (they are distinct, closed under addition/multiplication/inversion, and there are q^j roots). Its irreducible factors have degrees dividing j: the Frobenius orbit length of a root is exactly its minimal-polynomial degree. The stated congruence forces all irreducible factors of f to be distinct and have degrees dividing d; every proper divisor of d divides d/p for some prime p|d and is excluded by the gcd tests. Conversely an irreducible degree-d polynomial satisfies the tests.

Such a polynomial exists: in F_(q^d), the union of all proper subfields has at most sum_(e=1)^floor(d/2) q^e<q^d elements for d>=2, since every proper subfield degree is a proper divisor of d. An element outside has Frobenius orbit length d. Its minimal polynomial supplies an irreducible polynomial. The case d=1 is immediate. Enumeration of all q^d possibilities and modular powering/gcd takes polynomial time in Q,d,s. Prime divisors of d are obtained by trial division.

Factor Q-1 by trial division. In the constructed extension find the first nonzero g with g^((Q-1)/p)!=1 for every distinct prime p dividing Q-1. Such an element exists: for each prime power p^e dividing Q-1 the polynomial Z^((Q-1)/p)-1 has fewer than Q-1 roots, so choose a nonroot a. By Lagrange's theorem the order of a has full p-adic valuation e; raising a to (Q-1)/p^e gives an element of order p^e. The product of these elements for the different primes has order the product Q-1, since they commute and their orders are pairwise coprime. This proves cyclicity and the search criterion without assuming it. Exhaustive search over Q elements and trial division up to sqrt(Q-1) are polynomial time in Q. Multiplication by the first qualifying g gives the required A and its transitive nonzero orbit.

The guards Q,h^d<=(n+1)^6 bound these searches and enumeration of all monomials by fixed polynomial powers of n. Evaluating m polynomials at successive A^i v, followed by m traces, is polynomial time. The numerical guards give seed length (d+1)s<=2ds<=2cap(beta)log(n+1)<=12log(n+1) for every accepted parameter tuple. A sharper EVENTUAL estimate is available from the complete supplied prime-interval argument and will be used by the tree. Put W=log2 K. Uniformly for L^(2/3)<=W<=L^(19/20), the exact d=ceil(aSU*L/floor W) obeys d*W=(aSU+o(1))*L and d tends to infinity. Since ceil_root(n,d)=n^(1/d)*(1+o(1)) in this range, log2 h=L/d+log2 d+o(1)=(1/aSU+o(1))*W. Also log2 m=(bSU+o(1))*W, and log L=o(W). For an upper comparison put T=log2(C0*(L+1)^2)+O(1)=((1/aSU+4bSU)+o(1))*W,and let s_old be its least odd prime. The supplied relative-prime-interval theorem gives s_old/T->1 uniformly here. Since s_old=O(W)<=(L+1)^2 eventually,2^s_old>C0*(L+1)^2>=C0*s_old,so the new degree exists and s<=s_old. Conversely2^s>C0*s implies s>log2 C0=((1/aSU+4bSU)+o(1))*W. This sandwiches s at the same leading coefficient. The search visits at most(L+1)^2 prime candidates before guard failure,and the supplied prime routine is polynomial in its argument. Thus runtime remains polynomial even on rejection.

Every tuple accepted by the previous SU constructor is accepted now with the same m,d,h and no larger s or seed. Its old degree satisfies the actual inequality,so the new degree exists below it. All q^d and s-dependent rank/advice guards become easier. The new inequality gives q>C0*s>h+1 and q>m,so the other guards also hold. This asserts component acceptance/seed dominance,not pointwise circuit-complexity ordering of changed languages. The FULL seed therefore satisfies

 (d+1)*s=((1+4aSU*bSU)+o(1))*L.

All error estimates are uniform throughout this entropy window, since d grows at least as a positive power of L and W/log L tends to infinity uniformly. For the tree's L0=ceil(log2 n), the distinction L=L0 or L0+1 only contributes a relative o(1). Therefore the seed length is eventually at most(seedguard)*L0, since(1+4aSU*bSU)<seedguard; in particular it is at most6L0. This sharper estimate uses the full supplied relative-prime-interval proof, not a coarse factor2 interval or an external prime oracle. No reconstruction or received-word enumeration is performed by the evaluator. The proof's very long curve sequence and predictor tables are information-theoretic counting devices only.


# Small-error linear extraction with a polynomial seed budget

Let a source consist of N bits, with entropy at least r/2+floor(r/4), where1<=r<=N and ceil(log2(N+1))<=floor(r^(1/2000)). For all sufficiently large parameters, construct a linear strong extractor whose output length is M=floor(r/4), error at most2^(-r^(zeta)/3), and seed length O(r^(nu)). All source hypotheses here are AFFINE when output iteration is used. Integer powers in the construction mean exact floor roots of integer powers; floating-point arithmetic is unnecessary.

First make a smaller-output extractor. Put z=floor(r^(zeta)), epsilon=2^-z, m=floor(r^(999999/1000000)), rho=epsilon/(2m). Use the locally evaluated code with index length ell=2s and list bound2^s from the corresponding section of extractors.md. Then ell=O(log N+z+log m)=O(r^(zeta)), and ell=Omega(r^(zeta)), because log N<=r^(1/2000)=o(r^(zeta)).

Choose a prime degree p>=ceil(log2 ell), using the least such prime, and construct F_Q with Q=2^p. Here the complete supplied relative-prime-interval theorem gives a sharper seed bound than a factor2 prime interval. Set T=ceil(log2 ell). For each fixed delta>0, there is eventually a prime between T and (1+delta)T, so the least prime p>=T obeys p/T->1. As T/log2 ell->1, log2 Q/log2 ell->1. In particular, eventually Q<=ell^(2000001/2000000), while Q>=ell always. The code uses exactly this least-prime choice. These statements apply uniformly as ell tends to infinity; no new prime theorem or oracle is assumed. Let A=ceil(log_Q m), meaning the least integer A with Q^A>=m. Eventually A is bounded by a FIXED constant, since log Q>=log ell=Omega(log r). Also A<=ell. Encode integers0,...,m-1 as distinct coefficient vectors of polynomials P_i in F_Q[T] of degree at most A-1.

For the first ell field elements alpha in their binary integer order define the ordered set

 S_i={(alpha,P_i(alpha)):alpha among those ell elements}.

The universe is F_Q^2, of size d=Q^2<=ell^(2000001/1000000)=O(r^(sigma)). Each set has EXACTLY ell elements, not Q padded elements. Distinct polynomials agree at at most A-1 field points, by the root bound. Hence the reconstruction advice B is at most m*2^(A-1)=O(r^(999999/1000000)), with a fixed, possibly large constant. The other costs d, s and log2(4m/epsilon) are O(r^(sigma)). Their sum with B is o(r), therefore below r/2 eventually. The reconstruction lemma proves strong error at most epsilon. For fixed seed this extractor is linear in the source input.

Produce M bits by J=ceil(M/m) successive applications with independent fresh seeds, trimming the final output block. Since m>=r^(999999/1000000)/2 eventually, J=O(r^(1/1000000)). The exact advice guard already ensures m<=r/2 for every accepted tuple, and eventually m=o(r), so this larger smaller-output map never asks for more uniform bits than its entropy hypothesis permits. At each stage, condition on all preceding seeds and output blocks. For an affine initial source, this is either an impossible history or a uniform affine source whose entropy dropped by no more than the number of preceding output bits: every map is linear once the previous seeds are fixed. Thus its entropy stays at least r/2 as long as the initial entropy was at least r/2+M. The applications below verify this extra condition from their explicitly charged kernel entropy bounds. More generally the iterated guarantee needs initial entropy at least r/2+M; do not erase this extra requirement.

For that initial entropy, the strong-extractor conclusion applies to every possible history. Replacing the successive output blocks by fresh independent uniform blocks, one at a time, incurs total variation at most J*epsilon by the triangle inequality and conditioning. Therefore the concatenated output is a strong extractor with error at most J*2^-z<=2^(-r^(zeta)/3) eventually. Its total seed length is J*d=O(r^(1/1000000+sigma))=O(r^(nu)), and it is linear for every fixed complete seed. Every design and code bit is explicitly locally evaluated in time polynomial in N and r. Huge numbers q in the code specify index spaces only; their full tables are never enumerated.



# Augmented row condensation with a raw entropy reservoir

The input to this component is the PAIR consisting of an n-bit affine source X of entropy at least K and t rows W, each an affine function of X, with one uniform R-bit row. The evaluator receives the values of X and W, not their unknown affine parameterization. This extra raw source is essential. We do not assert the following conclusion for a matrix alone.

Put R=floor(K^(bSU)), g=floor(K^((1-bSU))), w=floor(K^((bSU-1/1000000))), M=floor(K/4). For every1<=t<=g and0<=take<=min(R,M), the supplied deterministic algorithm outputs take bits within error2^(-K^(eGroup)) of uniform at all sufficiently large parameters satisfying ceil(log2(n+1))<=floor(K^(1/2000)). Only take equal to the smaller-output transfer extractor's full seed length is needed below.

## Algorithm and exact rank budget

At each round slice the first w bits of every current row, pair consecutive slices, duplicating the last when their number is odd. Apply the supplied half-entropy cubic extractor to each pair and trim its output to the full seed length of the NEAR-linear-output extractor on the ORIGINAL n-bit X at PARAMETER K. Apply that seeded map to the original X once for each pair. Each new row has exactly M bits. Continue until one row remains and return its first take bits. When t=1 initially, simply return the requested prefix of its uniform row.

Every round uses the same raw input X, the same parameter K, and the same fixed slice length w. Initial rows have length R>w eventually sincebSU>(bSU-1/1000000); later rows have length M>w. With t0=t and J=ceil(log2 t), the exact row counts are t_j=ceil(t0/2^j) for0<=j<=J. Hence

 sum_(j=0)^(J-1) t_j <=2t0+J,
 (2t0+J)*w=O(K^(999999/1000000))=o(K).

The evaluator requires the exact integer guard(2t0+ceil(log2 t0))*w<=floor(K/4), along with w<=min(R,M), take<=min(R,M), and the supplied near/cubic parameter guards. Every guard succeeds eventually: the full near seed has length O(K^(nu)), while cubic output has Theta(K^((bSU-1/1000000)/2)) bits andnu<(bSU-1/1000000)/2. The small-error helper's full-output error is eta<=2^(-K^(zeta)/3). All searches and locally evaluated code bits are polynomial in n with a fixed exponent.

## Conditional affine invariant and seed independence

Maintain a mixture of joint affine distributions of(X,current rows) in which the raw X has entropy at least K minus the total number of previously exposed slice bits, current rows are affine functions of X, and at least one row is uniform at its full current length. Initially this is exactly the hypothesis. The mixture approximation incurs the stated error; it is not an assumption that the unconditional nonlinear algorithm always has affine output.

Fix one such mixture component at a round. Let Z collect all t_j*w current slice bits. These are affine linear observations of X. Choose a basis for the kernel of the observation map on the direction of X and extend it. This gives X=A(Z)+B with B an affine source independent of Z. On every possible fiber its entropy is at least

 K-(sum_(h<j) t_h*w)-t_j*w >=K-floor(K/4) >=K/2+M.

Consequently the supplied near-linear-output strong extractor applies to B at PARAMETER K. This is the point of retaining raw X: its entropy is not limited to the current good row length, and the output is always M rather than one fourth of the previous row length.

A pair containing the current uniform row has slice-pair entropy at least w. Its cubic seed S is etaCubic-close to uniform, etaCubic=2*2^(-sqrt(w)/50). Since S is a function of Z, it is independent of B before conditioning on Z. For every fixed seed s, the output on B is affine, so its variation distance from uniform is either zero or at least1/2. Strong extraction with a uniform seed and error eta implies that at most2eta of seeds are bad. Under S, that probability is at most2eta+etaCubic. A translation A(Z) changes only the fixed affine offset of the output, so this same bad-seed set applies to the actual X after conditioning on Z.

For every good Z value, all seeds and A(Z) are fixed. Every next row is affine in B, the specified good pair gives a uniform M-bit row, and the raw X still has the displayed entropy. Thus the conditional joint distribution satisfies the invariant for the next round, with the additional slice-bit budget charged. On exceptional Z values, replace the joint distribution by any affine distribution with the same affine X marginal and one uniform M-bit row given as an affine projection of X. Such a projection exists because H(X)>=K/2+M>M. Other rows can be arbitrary affine functions. This replacement preserves the invariant and changes total variation by at most2eta+etaCubic. Future computation is a fixed deterministic function of the joint pair(X,rows), so data processing and convexity justify iterating the bound even on replacement components. No independence from the whole raw X is claimed for S.

After J=O(log K) rounds the total error is at most J*(2eta+etaCubic)<=2^(-K^(eGroup)) eventually, since eGroup<zeta<(bSU-1/1000000)/2. A final prefix cannot increase variation. The t=1 case is exact. This proves the augmented row lemma with all history conditioning and raw-entropy losses charged. It uses no external somewhere-random extractor theorem.


# A block-source extractor using the augmented row map

Let the joint source(X1,...,Xb) be affine, each block n bits, and H(Xi|X1,...,X_(i-1))>=K on every possible full-history fiber. Require L^(2/3)<=log K<=L^(19/20), L=ceil(log2(n+1)); the numerical entropy guard of the small-error helper holds eventually in this window. Use the supplied constant-error strong linear extractor on X1 and enumerate every full seed. Trim each output to R=floor(K^(bSU)), forming a v-by-R affine matrix W1, v=2^d0. At least one row is uniform because each fixed-seed affine output is either uniform or at distance at least1/2, while the average strong error is at most1/4.

Let g=floor(K^((1-bSU))), eventually at least2. Let E be the SMALLER-output small-error extractor at parameter K and input length n, trimmed to R output bits. Its initial output floor(K^(999999/1000000)) exceeds R; its full seed length is d1=O(K^(sigma)), and strong error is epsilon1=2^-floor(K^(zeta)). It needs only entropy K/2. Do not substitute the near-linear-output routine for this transfer: that longer seed is used only INSIDE the augmented row map.

To turn Wi into Wi+1, partition its rows into consecutive groups of g, padding the last group with zero R-bit rows. For each group compute the augmented row map on the PAIR(Xi,group), trimming to d1 bits. This trim is available eventually because d1=o(R) and d1=o(M), M=floor(K/4). Use each resulting seed in E on the NEXT block Xi+1, giving the next R-bit matrix rows. The raw block supplied to the group map is Xi, not Xi+1. Code retains it as previous while advancing to the next block.

Choose b=1+min{j:g^j>=v}. The row-count recurrence gives one final row after b-1 transfers. Since log g=Theta(log K) and d0=O(log n), b=O(log n/log K)=O((log n)^(1/3)) throughout the entropy window. Output that final R-bit row or any shorter requested prefix.

## Full-history conditional invariant

There is a set of full histories(X1,...,X_(i-1)) of probability at least1-(i-1)*tau such that, conditioned on each such history, Wi is an affine function of Xi and has one uniform row. Here tau=2epsilon1+epsilon2, epsilon2<=2^(-K^(eGroup)). Initially this holds exactly.

Fix a good history. The raw Xi has entropy at least K by the block-source hypothesis. Every current matrix row is an affine function of Xi because fixing the previous input blocks fixes all previous transfer seeds. At least one group has a uniform R-bit row. The augmented row lemma therefore shows that THIS group's seed, a function of Xi, is epsilon2-close to uniform. All other group outputs may be correlated and arbitrary functions of Xi.

The joint pair(Xi,Xi+1) on this history is affine. Decompose Xi+1=A(Xi)+B by a kernel/complement split, where B is affine, independent of Xi and H(B)=H(Xi+1|Xi,history)>=K. The good group's seed is independent of B, being a function of Xi. Strong linear extraction on B implies that at most2epsilon1 of uniform seeds fail to yield a uniform R-bit output. Thus its actual seed is bad with probability at most tau. For each Xi value where it is good, conditioning on Xi fixes all group seeds and A(Xi); every next row is affine in B, the specified next row is uniform, and the full block-source hypothesis remains available at the next step. This proves the invariant after discarding probability at most tau from each previously good history. No group output is assumed to reveal the whole Xi; full original blocks are conditioned.

The final output is within(b-1)*tau of uniform. Since the group error exponent eGroup is strictly above eBE, and the transfer error exponent zeta is larger still, the growing b is absorbed to give final error at most2^(-K^(eBE)) eventually. Output length is floor(K^(bSU)).

## Uniform runtime and scope

Initial seed enumeration is polynomial in n by the supplied SU field/seed guard. Each augmented group uses O(log K) rounds, explicitly supplied cubic field construction and locally evaluated linear extractors; its raw input length is n. The group count decreases by g up to a ceiling on every transfer. Numerical ranks and guards depend only on n,K and row counts, not on target bit values. Every step therefore has fixed polynomial runtime and is total under the stated finite fallbacks. Padding a short subblock to ambient n bits preserves affinity and entropy and never adds ignored target inputs.

# Auxiliary fixed-BASE numerical ledger

This is the separately scoped B_N calculation referenced from proof.md,not the current controller definition. All its mathematics is retained; current six-tuple and minimum-certificate proofs are in proof.md.

## BASE analysis: Numerical feasibility with uniformly bounded rational rates

Write S=13Z,w0=33743/10000=3.3743,wm=527/500=1.054,a=27/50=.54,delta=25001/100000. Wmin=ceil(wm Z),Wmax=ceil(w0 Z),W_j=max(Wmin,Wmax-ceil(delta jB)). Uniformly over planned depths for sufficiently large actual N,Z=(1+o(1))L^(3/4)B^(1/4),base=(1+o(1))L^(1/4)B^(3/4),ell_log=(7/100+o(1))sqrt(LB),X=(25927/25000+o(1))Z. Lpad=ceil(Wmax^(20/19))=o(L),and Wmin^3>>(L+1)^2.

The exact identity L_j=max(L-V_j,Lpad),with Lc_j=L_j or L_j+1,gives (Lc_j)^(2/3)<=W_j<=(Lc_j)^(19/20). The explicit uniform rational-rate component proof therefore applies at every level; the context-size guard was proved above. Full-seed coefficient at most5.333355 has an exact guard gap1/100000,and guard at most5.333365<6. Q,h^d<=(n_j+1)^6 hold uniformly. Response seed count is at most N^6 eventually,and at most(N+1)^12 by numerical guards at every accepted length. Every response/block output has logarithm at least beta_j W_j-O(1)>=.54Z-O(1).

For all allowed contexts .08<c<.5. The fixed height multiplier7/100 gives (7/100)c<.035<.06. Let M=Amax include terminal arity. Then

 M<=.06(Wmax/Z)base L/Lpad+1=O(L^(35/76)B^(37/76)),

so M=o(Z),M^2=o(base L),M<X,and X+M<=1.0371Z eventually. Consecutive active log gaps are at most X+M,and H'ell_log<=(1750/25927+o(1))Z. First active crossing has shortest block log L-X-M+O(1)>>S/1000. If another active depth follows H,maximality gives maxleaf<=2^(X+M)(k^(1/1000)+1)+1. Otherwise V_J-V_H<X and N<2^(V_J+A_J),giving maxleaf<=2^(X+M)+1. Both cases imply maxleaf<=2^((1.0371+.013+o(1))Z)<2^(1.051Z),strictly below hidden entropy2^(1.0542Z).

For u_j=L-V_j>=Lpad,the exact floor gives

 A_j>=(7/100)c_j W_j base L/(Z u_j)-1,
 u_j^2-u_(j+1)^2>=(7/50)c_j(W_j/Z)base L-2L-M^2.

Below clipping,each step is at least floor(.005Wmin base L/(Z Lpad)),so only O(Lpad^2/(base L))+1=o(L/base) additional steps occur. Put h0=(w0-wm)/delta=232030/25001. Here wm>a/(13/25),so the upper clamp never applies. The strict bound w0<a/(3/20) ensures that the lower clamp also never applies for sufficiently large N. Thus beta_j=aZ/W_j throughout every planned level eventually;clamping still defines the total algorithm at finite exceptions.

For wm<=w<=w0 put beta(w)=a/w,alpha(w)=1/(1-beta(w))+1/100000,

 c(w)=(9999/10000)(1-beta(w))/(1+4alpha(w)beta(w)+1/100000),
 f(w)=w*c(w).

The function f is increasing: beta decreases with w;1-beta increases,and the positive denominator1+4beta/(1-beta)+4beta/100000+1/100000 decreases. Thus c and w are positive increasing factors. No differential equation or numerical integration is assumed.

Use4000 FIXED boxes. Set Delta=(w0-wm)/4000=23203/40000000,w_r=wm+rDelta for0<=r<4000,and A_r=floor(10^12*f(w_r)). These are explicit rational/integer expressions independent of N. Exact finite arithmetic gives

 SUM_(r=0)^3999 A_r=3080208477059576,
 (7Delta/(50delta))*SUM A_r/10^12
   =62536317631561674187/62502500000000000000>1.

The declared audit recomputes every rational value,checks each floor by both integer inequalities,and verifies this exact sum and ratio. This is a finite integer certificate for a displayed constant,not empirical evidence for a universal tree claim.

Suppose J>=(23203/2500)L^(3/4)/B^(3/4). Since23203/2500>h0,the first m=ceil(h0*Z/B) planned indices then exist for all sufficiently large N. Before deleting clipped indices,each fixed box contains Delta*Z/(delta B)+O(1) indices of the raw schedule Wmax-ceil(delta jB);the possible final floor point changes only O(1) terms. On a retained index in[w_r,w_(r+1)),monotonicity gives f(W_j/Z)>=f(w_r)>=A_r/10^12. Deleting the o(L/base) clipped indices changes the normalized sum by o(1),because f is uniformly bounded on this fixed interval. Since base*Z/(LB)->1,and the arity multiplier is7/100,the remaining normalized potential drop is at least the displayed coefficient minus o(1),which exceeds1 and contradicts initial potential L^2. Count-rounding errors across4000 fixed boxes,and the terms M^2=o(base L),also total o(1). Thus H<=J<(23203/2500)L^(3/4)/B^(3/4)=9.2812L^(3/4)/B^(3/4) eventually.

At each depth,log g_j=(1-beta_j+o(1))W_j and b_j<=2+d_j/log g_j. Because (7/100)*base*L/Z=(1+o(1))*ell_log,the exact coefficient identity gives

 b_j A_j<=[c_j*seedguard_j/(1-beta_j)+o(1)]ell_log+2M
          =(.9999+o(1))ell_log+2M<.99995ell_log.

Therefore sum_j binom(t_j,b_j)<=ell eventually. Since1-beta_j>=.48,b_j=O(L/Wmin) uniformly,so log2(4b_j)<=delta B and

 log2 product_(j<H)(4b_j)<=(23203/2500)delta Z
   =580098203Z/250000000<2.3204Z.

Target dstar>=2^(4.4119Z),active log gap<=1.0371Z,and log2(4^H')=o(Z) give active stopping entropy at least2^(3.3747Z). A nice source at ordinary depth j descended from active depth jstart retains at least2^(3.3747Z-delta(j-jstart)B)-Dtot. Its log margin against W_j's declining branch is .0004Z-O(1);against the floor it is at least(3.3747-2.3204-1.054)Z-O(1)=.0003Z-O(1). One further ordinary division costs delta B=o(Z),and hidden-child division costs at most M+2=o(Z). Every interleaved forcing and hidden source therefore has PER-NODE entropy at least2^(W_j+.0001Z),after all charged ranks. No global minimum substitutes for a larger earlier threshold.

Forcing rank Dtot<=H*tmax*ell*Pmax<=2^(.137Z),and later resilience rank(H+1)ell Pmax<=2^(.137Z),are negligible relative to every K_j. Conditional entropy H(X_g|X_g') stays above2^(3.3746Z)>K_g;globally hidden entropy stays above2^(1.0542Z)>maxleaf2^(1.051Z).

For every context,eBE=beta_j/4-102/100000000. Since beta_j W_j>=aZ and W_j<=w0 Z+1,

 eBE W_j>=aZ/4-(102/100000000)(w0 Z+1)
          =.134996558214Z-.00000102>=.134996Z eventually.

Each block error is at most2^(-2^(.134996Z))<2^-Pmax. Output logarithm.54Z-O(1) exceeds log2(ell Pmax). Since(2000001/1000000)(1750/25927)<.134995<.1349955,log2(ell^3 p2)<=.134995Z<log2 Pmax eventually. Confidence grows faster than every fixed polynomial in L;R<=N^6 gives H log2(2R+2)=O(L^2)=o(Pmax). All uniform component/helper/row/rare-event/numerical guards hold for all sufficiently large actual integer lengths.


# Auxiliary BASE concrete-tree definition

This complete definition is the separately scoped B_N referenced from proof.md.

# BASE analysis: A concrete tree disperser with local component lengths

At actual input length N>=2 put L=ceil(log2 N), B=ceil(log2(L+1)), Z=ceil(L^(3/4))*ceil(B^(1/4)), base=ceil(L^(1/4))*ceil(B^(3/4)), E=13Z, Wmin=ceil(527Z/500), Wmax=ceil(33743Z/10000), k=2^E, ell_log=ceil(7*ceil(sqrt L)*ceil(sqrt B)/100), ell=2^ell_log, Pmax=2^ceil(269991Z/2000000), and X=ceil(25927Z/25000). Put Lpad=ceil(Wmax^(20/19)), computed by the least integer whose19th power is at least Wmax^20. Before allocating component thresholds or2^Lpad require E<=floor(log2 N), Lpad<=floor(log2 N), Wmin^3>=(L+1)^2. Otherwise return0. These guards hold eventually and prevent disproportionate allocations on small inputs. Define output0 also at N=0,1.

Plan levels beginning V_0=0. At each planned depth j, define the deterministic local ambient length

 n_j=max(ceil(N/2^V_j),2^Lpad), L_j=ceil(log2 n_j),
 W_j=max(Wmin,Wmax-ceil(25001*j*B/100000)), K_j=2^W_j,
 beta_j=max(3/20,min(13/25,27Z/(50W_j))),
 mode_j=the reduced positive numerator/denominator pair of beta_j,
 c_(mode_j)=(9999/10000)(1-beta_j)/(1+4alpha_j*beta_j+1/100000),
 alpha_j=1/(1-beta_j)+1/100000,
 A_j=max(1,floor(7*c_(mode_j)*W_j*base*L/(100*Z*L_j))), t_j=2^A_j.

Since W_j>=Wmin=ceil(527Z/500)>27Z/(50*(13/25)),the upper clamp never applies and beta_j W_j>=27Z/50. Use this context throughout every extractor,helper,reservoir and block-source invocation at that depth. The complete rational parameter formulas are specified in extractors.md and computed by regime in family.py. Context is a pure explicit parameter,including every cache key. Before any variable integer powers,component constructors require max(pRate,qRate)<=100Lc_j and3/20<=pRate/qRate<=13/25. This guard holds for every tree context: before reduction pRate=27Z<26Wmin<=26L_j and qRate=50W_j<=50L_j,while clamped contexts(3,20),(13,25) also fit. Here W_j<=Wmax<=Lpad^(19/20)<=L_j follows from the clipping identity. Reducing by gcd never increases either integer. No source value,epsilon,advice or external field data selects the context.

If V_j+A_j<=floor(log2 N), append V_(j+1)=V_j+A_j and continue; otherwise stop at J=j. Retain this terminal A_J too. Thus every level through J has nonempty intervals, and N<2^(V_J+A_J). Let Amax=max_{0<=j<=J}A_j and tmax=2^Amax. Require10000*(X+Amax)<=10371Z, otherwise return0.

For every j<J build the supplied constant-error extractor at component input length n_j and threshold K_j. Its local degree uses Lc_j=ceil(log2(n_j+1)), DISTINCT from L_j and the global L. Require its full seed count d_j to satisfy d_j<=seedguard_(mode_j)*L_j, with the exact rational context guards in extractors.md and its output length to be at least ell*Pmax. Trim this linear map to an ell by Pmax rectangle, denoted E_j. Build the supplied affine block-source extractor at block length n_j and threshold K_j; let b_j be its prescribed block count, and require at least Pmax output bits. Every component's exact integer numerical guard must succeed; otherwise return0. These are parameter checks, not tests of extraction on particular inputs.

Let the preliminary active levels be0 and, for each positive integer multiple iX reached by the V_j, the FIRST depth j with V_j>=iX. Under the accepted gap guard Amax<X, these crossing depths are distinct. Write them a_0=0<a_1<... . Choose H to be the last preliminary active depth such that floor(N/2^V_H)^1000>=k, and retain precisely the active depths through H. If fewer than two exist, return0. Write their count minus one as H'. Require t_j>=b_j for j<H and sum_{j=0}^{H-1}binom(t_j,b_j)<=ell, otherwise return0.

At level j, node i is the actual half-open interval[floor(iN/2^V_j),floor((i+1)N/2^V_j)),0<=i<2^V_j. Its t_j children partition it EXACTLY. Bits outside the interval are never interpreted as additional target bits. For a component invocation at this node, put the interval's actual bits in the PREFIX of a vector of length n_j and append zeros. For its block-source rows, put each selected child's actual bits in a separate prefix-zero-padded vector of length n_j. These are injective linear embeddings of the respective actual block coordinates. They preserve all affine entropies and conditional entropies; no target domain is padded, ignored or restricted.

For u=0,...,H'-1 set Tconf=ceil(ell_log/2000000) and

 p_u^(1)=2^(4000002*(H'-u)*Tconf),
 p_u^(2)=2^(4000002*(H'-u)*Tconf+2000001*Tconf).

Require all p2<=Pmax, otherwise return0. The exact ratios p2/(ell*p1) and earlier p1/(ell*later p2) equal G=2^(2000001*Tconf-ell_log)>=2^Tconf. At a FIXED node all confidences use the SAME map E_j and SAME rectangle width Pmax. Different levels have different component ambient lengths and may use different E_j. The common response seed-count upper bound R=max_{j<H}2^d_j will suffice throughout the proof.

Compute challenge matrices bottom-up. Leaf matrices are empty. At an internal level-j node, test its children in order at confidence Pmax using E_j on the prefix-padded node input. Select the first PASS child, or the first child if none passes. Concatenate the selected child's challenge with the Pmax-bit rows BE_j on EVERY lexicographically ordered b_j-tuple of children, each child padded to n_j. The total number of rows is at most sum_{r=j}^{H-1}binom(t_r,b_r)<=ell and is deterministic at each depth. All tests pad the challenge matrix with zero rows to ell. Follow the resulting selected path. At active depth a_u compare the selected descendant at depth a_(u+1) using the two specified confidence tests with map E_(a_u). At the first first-confidence FAIL output the second-confidence result, encoded FAIL0/PASS1. If there is none output0. This is exactly the fixed BASE evaluator.

All planned ambient lengths are <=N. The component input embeddings are linear and injective on their actual blocks. The tree has at most N nonempty nodes per depth and J<=L, hence O(N log N) nodes, a sufficient polynomial bound even with unary levels. The tuple-count guard bounds each enumeration by ell<=N eventually; otherwise ell itself is at most a fixed polynomial in N by its formula. Each response has at most(N+1)^12 seeds by its Q^6 finite-field guard, and all component evaluators have a fixed polynomial time exponent at lengths <=N. Planning and guards require at most L levels, only polynomial-size integers, and polynomial-time component constructors. Thus all successful evaluations and all fallback evaluations take time polynomial in actual N, with a fixed program and exponent independent of epsilon. The complete refined full-B2 transfer applies to dispersion at dstar=2^ceil(44119Z/10000), giving4N-8dstar+4 and every stated same-family4N-N^epsilon bound.


# Response maps specialized to their actual public output

The default optimal_parameters and constant_parameters retain m=floor(K^beta) and all displayed default seed asymptotics. The optional public take argument instead requires an integer1<=take<=floor(K^beta) and sets m=take BEFORE defining the field degree. It retains the SAME d,h,mode/context guard and every actual finite predictor/field/advice/cap check. The general finite SU proof above applies to this m directly:its reconstruction uses the actual m in every error and list inequality,and no equality m=K^beta is used by that finite proof. Linearity in the raw input and strong error1/4 therefore hold at the stated source entropy. response_parameters and response_extract use this form at take=ell*P. Block extraction keeps the default capacity and aliases unchanged. All cached arguments include the public take;source bits never choose it.

Whenever the default constructor passes and take<=its m,the smaller target weakly decreases C0,hence its least valid odd prime degree and full seed. The finite guards get easier as in the preceding least-degree proof:q^d and advice depend monotonically on s and m,while q>C0*s ensures q>h+1 and q>m. Thus specialized response construction succeeds with no larger seed. The SAME default upper seed/asymptotic bounds remain sufficient;no default seed equality is asserted for smaller outputs. The loop and fields only become smaller,so all accepting/rejecting paths retain fixed polynomial time. A new response map can also pass finite guards when the default fails;its directly checked general finite theorem supplies safety in that case.
