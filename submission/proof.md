# One fixed uniform language with log(N)^11 loglog(N)^10 loss

For the single total language evaluated by family.py, we prove that absolute B>0 and N0 exist such that for EVERY integer actual input length N>=N0,

 C_B2(F_N)>=4N-B*(log2 N)^11*(log2(log2 N))^10.

Both canonical claims assert this complete theorem. Its protected expression and projection are(4,1); no larger coefficient, exponent or endpoint is claimed. All logarithms below are base2. The complete extractor arguments are supplied once in disperser.md; the full-B2 gate proof is supplied once in elimination.md. Their source versions and provenance are identified there, not used as unexplained premises.

## The one total program and its actual input

family.py contains the complete evaluator, including PolyPlan/PolyEngine and all field/extractor code in the SAME module. There is no separately loaded engine or external construction data. On a Boolean string x of actual length N, evaluate(x) returns0 for N<128, including N=0. Every larger length runs the scalar guarded construction. There is no legacy finite-field evaluator.

At all other lengths put L=ceil(log2 N), z=ceil(log2 L). The first L input bits encode S-1, in little-endian order. Every1<=S<=2^L is evaluated. The exact power-of-two scales are

 h=e=P2(S L^2z^2), ell=P2(S^3L^4z^4),
 eprime=P2(S^5L^6z^6), d=P2(S^9L^6z^6),
 Y=P2(S^11L^9z^8), k=P2(S^13L^11z^10).

PolyPlan first computes log2 k and rejects k>N before expansion. It builds the stated scalar plans and checks all exact source, seed and prefix guards. A rejected large plan returns0. An admitted plan reads only the initial A1 package and B5 hash modulus AFTER the header. An A1 package is exactly(D^2+1)s bits, with fixed D=12: lower coefficients of its monic field polynomial and every matrix entry. B5 uses1000*ceil(log2 r) modulus bits. Every other internal field is deterministic, and all five internal maps have zero data length. PolyEngine checks its data cursor and removes the ENTIRE header-plus-data prefix P before padding the payload with P zeros to N coordinates.

All good-prefix choices are actual input restrictions, never advice or witnesses. Bad initial/B5 field strings remain bounded ring/matrix computations with a defined Boolean output; goodness is not promised or tested. Initial good data tuples exist by the supplied A1 argument and work simultaneously for every required source. B5 uses the same fixed good field for its prefix hash and block hashes. Field existence, deterministic tower arithmetic and relative moduli used elsewhere are proved in A1/B1.

## Supplied extraction and exact budget closure

The implementation uses these supplied components:

A1 gives the low-degree initial strong linear extractor and its reconstruction proof. B1 gives binary-linear rank condensation and deterministic fields. B2.2 gives the four half-output strong linear maps: precision E, output M, B0=E+ceil(log2 M)+10, K>=32B0, field degree P2(8*(B0+ceil(log2 nprime))), seed O(E^2) at each fixed scale. It now emits B0-bit hash blocks along ceil(M/B0) weak-design rows. The previous-output table has length<=1.5B0*i<1.5M, including the partial last block. Conditional entropy remains>=7B0, enough for full-seed block collision error2^-B0. This is valid for general conditional sources, not only affine ones.

B2.3 gives short E1 directly by multi-bit polynomial hashing on ambient N, without condensation or data. For source entropy K, output M and precision E, require K>=M+2E+10; field degree s=P2(8*(M+E+ceil(log2 n)+10)), seed2s, output the first M coordinates of z*P_X(y). Its full-seed collision bound is one half times sqrt(2^M*(2^-K+n/2^s)), which is below2^-E. For E1, M=ell,E=eprime,K=10d and seed length d. Since eprime dominates ell+L and d/eprime>=S^4/2, one sufficiently large fixed scale pays its exact seed/entropy guards.

B3--B4 give the refreshed affine look-ahead/index proof, using fixed-seed linearity and full-seed strong error. No future seed is fixed prematurely; A-contributions of all streams, including the target, precede the B-parts used for uniform seeds. Independence order h is separate from binary alphabet2. B5 gives the strong shifted-agreement generator. It now emits100w-bit hash blocks, w=ceil(log2 r), at field degree1000w, with ceil((m-c)/(100w)) NW rows. Its previous-block table remains<1.5m, conditional entropy>=1000w, and block error<r^-400. Thus its strong error/rank and shifted-agreement bounds are unchanged. B6 supplies the robust monotone function, exact integer bias threshold and O(r^2u) moment/error proof. B7 closes the following inequalities.

Power-of-two rounding gives e^2<=4S^2L^4z^4, paying half-output seeds with ell's extra S factor. Short hashing needs O(eprime), while d/eprime>=S^4/2. Its transcript precision eprime>=100*h*alphabet*ell+2e is paid by the extra S in eprime. With index rounds b<=B_index L,

 b*h*d<=4B_index*S^10L^9z^8,
 hY<=4S^12L^11z^10.

The scales Y and k each gain a power of S. Also q0=(4h+10)ell<=56S^4L^6z^6, so q0/d<=56/S^5<=1/100 for a sufficiently large FIXED S. This is NOT an o(d) assertion. It gives the required padding/rank condition d>=q0. All other source losses and10d have smaller L powers than k. Initial A1 uses one eventual call, since its advice exponent131/12 is below11; its actual seed d0<=10000 L is checked before enumeration. All four B2.2 entropy margins hold eventually.

The unchanged robust moment requirements are h,e>=constant*r^2u. B7 proves r<=C L and u<=C0*sqrt(S)*z^2 at fixed scale, so one sufficiently large S pays those finite constants and B5's squared-log seed. Exact layout discontinuities lose at most2^O_S(log r), and the bad-coordinate bound remains N^(-2+o(1)); the supplied comparison gives q<=S_groups^(1-delta). Choose ONE integer S0 above these finitely many requirements. It fits the header and every guard for ALL sufficiently large integer N. The evaluator never needs to find it.

## Uniform P, NP and the actual-length lower bound

For every admitted header, scalar dimensions have polynomial bit size in N; S<=2N and the early guard controls expansion. Admitted fields, matrices and weak designs have polynomial finite-parameter bit complexity. CE row construction uses exact binomial arithmetic and polynomial bit lengths. Initial A1 has fixed dimension and streams only its first coefficient monomials; no entire ambient grid is enumerated. Tower arithmetic uses an explicit stack. The d0 cap bounds2^d0 by a fixed polynomial in N, uniformly over headers. Layout bounds r by a fixed constant times log N, so robust loops and exact threshold integers also have polynomial size. Every rejected plan is0, every bad admitted datum is bounded, and the small cutoff is0. Thus one fixed polynomial-time algorithm evaluates F on every Boolean string.

Set L_yes={x:evaluate(x)=1}. Let W(N)=0 and V(x,w) accept exactly when w is empty and evaluate(x)=1. If x is in L_yes the empty witness is accepted; if an allowed witness is accepted the evaluator is1, so x is in L_yes. This gives both directions for ONE fixed NP verifier. Witness bits are not circuit inputs. Eligibility does not depend on a good S0 or on efficiently selecting field data.

Fix header S0-1 and good initial/B5 input data at each sufficiently large N. The complete prefix is P=O(L): initial dimension is fixed with field degree O(L), B5's modulus is O(loglog N), and the header has L bits. Padding the remaining payload injectively into N coordinates preserves every affine-flat dimension. The supplied B1--B7 arguments therefore give an affine disperser at dimension k on N-P actual variables. All these source/scale bounds hold at every sufficiently large integer length, not selected powers of two.

Apply elimination.md's complete theorem4n-8d to that restricted function. It allows every B2 truth table, constant/unary/terminal cases, repeated pins, unrestricted sharing and nonuniform circuits. Prefix restriction cannot increase gate count. Hence

 C_B2(F_N)>=4(N-P)-8k=4N-(4P+8k).

Since P=o(k), the loss is eventually<=9k<18S0^13 L^11z^10. Ratios L/log N and z/loglog N tend to1 at all integer lengths, so B=32S0^13 works after enlarging the threshold. This is not a numerical or practical bound on B. For the SAME family any fixed extra loglog exponent absorbs B, then any fixed extra log exponent absorbs the loglog factor, and any N^epsilon eventually dominates the polylogarithm. These prove exactly the further eventual inequalities in bound.json and the joint liminf projection(4,1).

## Provenance and scope

The scale-header/extractor architecture and supplied proofs derive from9597bd82-cdc2-4d11-a0e9-95eb31bc3a1d. B2.2's O(E^2) adaptation was first supplied in the unpublished failed-review candidate9b6e4ab5-9926-44fd-84cb-1e59d4cebb19; the direct short hash was supplied in510a59cf-e44e-48bb-ad85-cd80e6b19bf4. Their verdicts are not premises. Block recycling, the smaller weak-design universe and the merged evaluator derive from published predecessor19095349-8b20-42f2-849e-3591198e69ed. They reduce row counts without increasing the leading bound. The gate proof is retained from28feddfa; the charged-prefix restriction design traces to db5f9c6e with later modifications; its legacy evaluator and Frobenius-test code are no longer retained. Primary paper provenance is given in the complete supporting evidence.

The finite audit tests only its named predicates, including the zero cutoff and early rejection, field/rank-map identities, block linearity/truncation, root/output distributions, CE expectations, data offsets, index mocks and local gate-sharing patterns. Obsolete optional lemmas and extractor-tree fixtures were removed with their unused code. No finite computation proves the universal theorem, and the enormous eventual branch is not executed end to end. Independent critique is advisory; configured review, human approval and promotion remain separate.

Both conditional-expectation designs use T=ceil(5a^2/2)+a. The B5 proof bounds initial expected overlap weight by exp(2/5)<=337/225<3/2 using an exact rational series bound. Thus all existing3/2 table-length, entropy and full-seed arguments remain valid with the smaller universe. This reduces the seed/construction constant while preserving the stated asymptotic bounds.

## Smaller initial dimension and explicit resource constants

Only the initial A1 call changes dimension, from10000 to12. Its entropy/advice comparison becomes Y*k^(1/12)*L=o(k), with exact L exponent9+11/12+1=131/12<11. The generic A1 proof applies at this fixed positive dimension. Its field degree is at most(128+o(1))L, the matrix seed has13 field elements, and the condenser seed is at most(40+o(1))L. Thus d0<2000L eventually at every fixed scale, and the program's explicit10000L cap passes. The initial data package is145 field elements, so including header and B5 modulus gives P<=20000L eventually. The uniform layout bound improves to r<=2*10^16 L. These are proved eventual resource constants, not an end-to-end timing or practical hardness threshold. The main remainder and leading projection remain unchanged.
