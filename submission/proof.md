# A uniform certificate-prefix family: coefficient 4, exponent 1

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove

    C_B2(F_N) >= 4N - R(N),       R(N)/N -> 0.

The joint pair is (4,1). The leading coefficient is inherited mathematics.
degree.md reproduces the degree-sensitive full-B2 bound used below.
The leading gate argument is the
parameter method of Mishra and Santhanam, ECCC TR26-260, revision 1, November 3, 2026,
Sections 2–4. The algebraic construction is based on Ben-Sasson
and Kopparty (BK), *Affine Dispersers from Subspace Polynomials*, author
version March 12, 2010, Sections 3, 6 and 7. The necessary arguments are
supplied in our notation in `elimination.md` and algebra.md.

The paid input prefix supplies the field polynomial without an efficient
construction premise. Only a fixed-prefix restriction is a disperser.
All claims use the full B2 model and actual input length N.

## 1. One evaluator on every bit string

The first bit of each field block is its constant coefficient. On x of length N:

1. If N=0, output 0.
2. Let B be the smallest positive integer with B^5 > 32*N^3. Let k be the
   smallest prime at least B. Use deterministic trial division to find k.
3. If k>=N, output 0. Otherwise put m=N-k, r=ceil(m/k), t=ceil(sqrt(r)).
   If t>=r or k*t<r*(r+1), output 0.
4. Read the first k bits b_0,...,b_(k-1) as

       Q(Z)=Z^k + sum_(j=0)^(k-1) b_j Z^j    in F_2[Z].

   Check Z^(2^k)=Z modulo Q and gcd(Q,Z^2+Z)=1. If either fails, output 0.
5. Use K=F_2[Z]/(Q). Append rk-m zero bits to the remaining m payload bits
   and interpret the r blocks of k bits as a_1,...,a_r in K. This is an
   injective F_2-linear map of the payload.
6. For i=1,...,r set e_i=r+1-i and y_i=a_i^(2^(e_i+1)-1). Output the
   constant coefficient of the elementary symmetric polynomial of order t
   in y_1,...,y_r, in its representative of degree less than k.

For explicit evaluation, initialize E_0=1 and E_1=...=E_t=0. For each i,
update E_j <- E_j+y_i*E_(j-1), in decreasing order j=min(i,t),...,1.
Splitting subsets by whether they contain i proves this division-free
recurrence. Binary powering takes O(e_i+1) multiplications per y_i.

This specifies every choice and every input length. Q is part of x, not
advice or a value supplied by an unspecified construction algorithm.

## 2. Exact prefix test and existence of hard prefixes

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

## 3. Uniform P and NP with empty witnesses

Bertrand's prime theorem gives k<2B, and so

    2*N^(3/5) < k < 4*N^(3/5)+2.

The prime scan and trial division are polynomial in N. B is computed by
integer binary search and fifth-power comparisons. All subsequent integer,
polynomial and field operations have polynomial bit cost. Field multiplication
can use schoolbook O(k^2) operations; the powers take O(r^2) multiplications
in total and the symmetric recurrence O(rt). A coarse O((N+1)^6) bound
suffices for the whole algorithm, including rejected inputs.

Let A be this fixed evaluator, L={x:A(x)=1}, W(N)=0, and let V(x,w)
accept exactly when w is empty and A(x)=1. If x is in L the empty witness
is accepted; if a permitted witness is accepted then A(x)=1 and x is in L.
This proves both directions with one uniform polynomial-time verifier.
Witnesses are not inputs to the target circuits.

## 4. All sufficiently large integer lengths pass the guards

We have k=Theta(N^(3/5)), m=N-k, r=Theta(N^(2/5)), t=Theta(N^(1/5)).
Consequently k<N and t<r eventually. Also r<N/k+1<(1/2)N^(2/5)+1 and

    r(r+1)/t <= sqrt(r)(r+1)
              <= sqrt((1/2)N^(2/5)+1)*((1/2)N^(2/5)+2).

Dividing the last expression by N^(3/5) gives limit 1/(2*sqrt(2)),
whereas k/N^(3/5)>2. Thus k*t>=r(r+1) eventually. These estimates hold
for every sufficiently large integer N, not just selected lengths.


## 5. A padding-aware certificate for this same language

Fix any valid prefix. Let ell=m-(r-1)k, so1<=ell<=k. The payload image
lies in K^(r-1) times the ell-dimensional span of the first ell field
basis vectors. Put w_i=r+2-i and L_i=k for i<r, L_r=ell. Define

    H_c=sum_i min(w_i,L_i), S_c={i:L_i>w_i},
    T_c=sum of the t-1 largest numbers max(L_i-w_i,0),
    U_c=sum_(i in S_c) max(k+1-w_i,0),
    D_N=H_c+max(T_c,floor(U_c/(t+1)))+1,
    B_N=sum of the t smallest weights among {w_i:i in S_c}.

These parameters are used only to certify the unchanged evaluator.
Eventually k>r+1 and r>t, so S_c contains at least r-1>=t indices
and B_N is defined. Let H,T,U,D_full,B_full be the preceding full-block
values: H=r(r+3)/2, T=sum_(j=2)^t(k-j),
U=sum_(j=1)^r(k-j), D_full=H+max(T,floor(U/(t+1)))+1,
B_full=t(t+3)/2. Then H_c<=H,T_c<=T,U_c<=U and B_N>=B_full,
so D_N<=D_full. These comparisons concern all sufficiently large N.

Here is the complete adjustment to algebra.md's block selection.
For any affine space in the payload image, its triangular decomposition
has d_i=dim(V_i)<=L_i; for the last block this follows because the
kernel of projection onto the first r-1 blocks is contained in the
ell-dimensional padding subspace. Put delta_i=d_i-w_i and
S={i:delta_i>0}, a subset of S_c. Its positive deficit satisfies

    sum_(i in S)delta_i >= dim(A)-H_c.

If |S|<t this sum is at most T_c. Otherwise, if every t-subset I failed
the coefficient-isolation condition sum_I delta_i+max_I d_i>k+1,
algebra.md Section5's identical averaging gives

    (t+1)*sum_S delta_i <= sum_S(k+1-w_i) <= U_c.

Both possibilities contradict dim(A)>=D_N. Some t-subset therefore
satisfies the same isolation hypotheses. Sections3,6,7 of algebra.md
then give its nonzero coefficient, of binary weight sum_I w_i>=B_N.
Section9's Boolean-degree bridge proves degree at least B_N on A.
No stronger field lemma or independence premise has been introduced.

The previous bounds give D_N=O(N^(4/5)); also
B_full<=B_N<=B_full+t=O(N^(2/5)). Eventually4<=B_N<=D_N<m.
The last inequality follows from D_N=o(m); B_N<=D_N follows from
H_c including every eligible selected weight. Apply degree.md to the
m-bit payload restriction P, then recover k gates by Section5.1:

    C_B2(F_N)>=k+C_B2(P)>=4N-3k-8D_N+4B_N+3.           (16)

The complete loss L_new=3k+8D_N-4B_N-3 is o(N). Thus for every
epsilon>0 all sufficiently large integer N satisfy L_new<=4epsilon*N
and C_B2(F_N)>=(1-epsilon)*4N, jointly proving the full guarantee
and projection(4,1) for the one language.

The new loss is never larger than the immediate preceding certificate
L_prev=3k+8D_full-4B_full-3. It improves by an unbounded amount along
infinitely many lengths, as we now prove without prime-gap asymptotics.
Write p=k-ell for the number of padding bits. From k=Theta(N^(3/5))
and r=Theta(N^(2/5)) we have k=Theta(r^(3/2)). Since t=ceil(sqrt(r)),

    T-U/(t+1)>=-k/(t+1)+r(r+1)/(2(t+1))-t(t+1)/2+1
              =Theta(r^(3/2)) with positive leading term.

The negative terms are O(r), while r(r+1)/(2(t+1)) is asymptotic to
r^(3/2)/2. Reducing only the last block's capacity changes T by at
most t-1: the next available full block of weight t+1 replaces it.
Thus eventually both maxima are attained by T and T_c. If p>=t-1,
then T_c=T-(t-1). Hence D_full-D_N>=t-1 and

    L_prev-L_new>=8(t-1).

There are infinitely many such lengths. For each integer B let I_B
be the consecutive integer N with (B-1)^5<=32N^3<B^5. This is exactly
the evaluator's fifth-root threshold B. Its interval length is
Theta(B^(2/3)), whereas throughout it t=Theta(B^(1/3)); k, the least
prime at least B, is constant. Also its length is o(k), so r changes
at most once. If p<t-1 throughout, all its padding residues, which
decrease by one modulo k as N increases, would lie in fewer than
max_(I_B)t consecutive residues. A run that long is impossible.
Thus every sufficiently large I_B contains p>=t-1. Since t tends to
infinity, the strict improvement is unbounded along infinitely many
actual lengths. At unpadded lengths the certificates can agree.
These compare certified losses, not actual optimal circuit sizes.

### 5.1 Recover the paid prefix gates

Let P be the payload restriction at any valid prefix b. Eventually
P has degree at least B_N>=4 and is nonconstant. Every irreducible Q
of degree k>1 has Q(1)=1. Flipping any single coefficient bit toggles
Q(1) to0, making Q divisible by Z+1 and therefore rejected on every
payload. After fixing any subset of prefix coordinates to b, every
remaining prefix bit is still essential: set the other free prefix
bits to b and choose a payload with P=1; flipping that bit changes
output to0. This applies in any fixed order of the k coordinates.

For completeness, if a full-B2 circuit depends on x and its restriction
x=c is not a constant or literal, that restriction has a circuit with
at least one fewer gate. Normalize constant/unary internal gates by
substitution, absorbing complements into descendant tables. If a
complement reaches the output, complement its driving gate's table;
a literal output is the only exception. Choose a gate using x, which
exists by dependence. Fixing x makes this gate constant or unary in
its other pin (also with repeated pins). Remove it by the same rules.
A terminal complement is absorbed into its surviving driving gate;
if no driving gate survived, the restricted function would be a
constant or literal, excluded here. All tables remain among the16
allowed B2 tables, and sharing costs no extra gate: each affected
consumer table is changed in place. Additional simplifications only
reduce size further.

At each successive prefix fixing the remaining function has P as a
further restriction, hence cannot be constant or a literal. Apply this
gate-removal observation k times to any circuit for F_N. It proves

    C_B2(F_N)>=k+C_B2(P).

This is an explicit elimination argument, not direct-sum additivity.
Compared with simple restriction, (16) gains k gates. The prefix is still part of N;
its cost in the complete loss is3k rather than4k.

## 6. Sources and complete evidence

elimination.md supplies the unrestricted charged reduction and ledger;
degree.md supplies the degree and connectivity argument; algebra.md
supplies coefficient isolation, triangular mixing and the degree bridge.
Their material sources are declared separately. In particular8eef supplies
the sharper full-block selection and degree arguments, reproduced here.
The padding adjustment and prefix-gate recovery apply to the same function.
Its alternative optimized family is unused. No source verdict is a premise.
The earlier universal+8 theorem is a separate published result; it is
not a declared claim or premise of this submission.
