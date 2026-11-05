# One fixed uniform language with log(N)^9 loglog(N)^7 loss

For the total language evaluated by family.py we prove that absolute B>0,N0 exist such that for EVERY actual integer input length N>=N0,

 C_B2(F_N)>=4N-B*(log2 N)^9*(log2(log2 N))^7.

Both canonical claims assert this entire theorem and its specified eventual weakenings for this SAME family. The protected coefficient/exponent guarantee and projection remain(4,1). Logarithms are base2. Complete strong extraction arguments A1--A2 are in seeded.md; B1--B7 and all parameter closure are in disperser.md. The complete gate theorem is in elimination.md. This proof uses their supplied arguments, not predecessor verdicts as premises.

## The total program and actual input length

family.py supplies the complete fixed evaluator, including PolyPlan/PolyEngine, ring/field arithmetic, samplers and every map. No external engine, construction table or per-length advice is loaded. evaluate outputs0 for N<128. For every larger length set L=ceil(log2 N),z=ceil(log2 L),D_int=max(2,z),C=2^30. The first L input bits encode S-1 in little-endian order, so1<=S<=2^L<=2N. Every header is handled by the same program. With P2 denoting upward power-of-two rounding, define

 h=e=P2(S L^2 z), ell=P2(S^3 L^2 z^2),
 eprime=P2(S^5 L^4 z^3), d=P2(S^9 L^4 z^5),
 Y=P2(S^11 L^7 z^6), k=P2(S^13 L^9 z^7).

The binary index alphabet is a=2, distinct from independence order h. The refreshed index argument uses q0=(2h*a+C*h*z^2+10)ell. The planner first checks raw log2 k and rejects k>N before expansion. It checks every stated source, seed, prefix, scalar sampler and robust resource guard. Failed plans output0. Its d0<=5L cap precedes enumeration of all initial seeds.

For an admitted plan initial A1 and all B5 fields/helper are constructed canonically by C1 in seeded.md. The bounded exhaustive searches use no external oracle or per-length advice, and their uniform polynomial runtime bounds are proved there. Only the FOUR internal growing-dimension field/matrix packages are read after the header, in exact self.seeds consumption order T,W,E2,E3. Each hasD_int*s-2 actual bits, encoding a monic constant-one binary polynomial of degree D_int*s; C3/C4 decode and reconstruct the field and companion action by linear algebra. Initial A1, short E1 and B5 consume zero input data bits. The cursor is checked; the entire prefix P is removed and the N-P payload padded with P zeros to N coordinates.

Good primitive binary internal packages exist by A1 and C3 and are fixed only as charged hardness restrictions. Bad internal binary-polynomial strings give bounded reconstruction and ring/matrix computations, or the zero map at a failed rank guard. The canonical initial/B5 searches always produce good data by C1's exact tests/existence arguments. Deterministic towers for PV and the weak hash retain their B1 proof. This defines the language on every header/data string.

## The four internal maps and entropy budgets

B2.1 gives strong affine extraction by B1 rank condensation and repeated direct A1 on the SAME condensed source with fresh seeds. Previous outputs reduce conditional affine entropy by at most their actual length; all seed bits are retained. The supplied sequential coupling bounds full strong error. Fixed combined seeds give binary linear maps, including truncation. A longer stipulated seed has ignored independent padding bits.

Each internal map uses this supplied B2.1 routine at dimension D_int=max(2,z),precision E=h and source threshold K=C*E*z^2*M,C=2^30. Outputs are ell for W/T,q0 for E2,d for E3, and stipulated seed length ell. Each owns a charged D_int*s-2-bit primitive binary package reconstructed by C3. The old weak-design half-output routines are auxiliary components only, not the four defining maps.

B2.2 retains the growing dimension explicitly. All actual K=O_S(L^6 z^8), so condensed length nprime=O_S(K) gives ceil(nprime^(1/D_int))<=128 eventually for each fixed S. Thus H<=128z,R<=32z,s<5E. The first power-of-two output candidate passes threshold K/2 because advice16M R H s<=327680 M E z^2<K/2, with strict room for the additive logs. At most two calls suffice. B1 seed<5E and total actual seed<25zE; ell/E>=S^2*z/2 pays this once fixed S>=8. These constants do not hide a fixed-D assumption. Full-seed affine extraction and fixed-seed linearity follow B2.1 exactly.

Short E1 remains direct B2.3: full-seed polynomial hashing on ambient N with K=10d, output ell, precision eprime and seed d. The supplied collision inequality proves its error without a condenser or input data. Its seed O(eprime) fits d because d/eprime>=S^4*z^2/2. Its entropy margin follows from the same comparison.

B3--B4 apply to these five fixed-seed linear maps. The W/T threshold isC*h*z^2ell, E2 thresholdC*h*z^2q0 and E3 thresholdC*h*z^2d. The complete refreshed A/B conditioning proof is retained, including last-round refresh and fixing every designated A-contribution before the B-variable supplies a seed. It uses no half-output source assumption. q0-2h*(a-1)ell=(2h+C*h*z^2+10)ell meets the W/T threshold exactly with slack. The bounds Y>=1000bhd and Y-100bhd>=C*h*z^2max(d,q0) hold eventually because h*z^2*d/Y=O_S(z^2/L) tends to zero. The residual k-hY-100bhaell exceeds10d+C*h*z^2ell: hY/k<=4/S and every other term has smaller L power. The complete inequalities are given in B7.

All remaining transcript requirements close: q0/d<=8C/S^5 is at most1/100 for a sufficiently large fixed S, eprime/(h*a*ell) gains a power of S, b<=5L, Y/(bhd) and k/(hY) gain a power of S. One sufficiently large fixed scale meets them. There is no claim that every admitted header certifies hardness.

## Initial seeds, robust generator and one fixed hard header

Initial extraction retains dimension11, precision E0=ceil(L/28)+10, adaptive B1 allocation and the exact rational A1 reconstruction field condition. Its advice Y*k^(1/11)*L has L exponent7+9/11+1<9. Hence it is o(k), so one eventual call suffices. The unchanged exact field calculation gives s0<=(192/1267+o(1))L and d0<=(2304/1267+15/7+o(1))L<4L, while the executable cap remains5L. Initial rank-deficient row fraction is at most2N^-1/28.

B5 retains the FULLY SUPPLIED logarithmic-external-seed bootstrap. A2 samples two conditional entropy blocks with1024-wise polynomial position sampling, explicitly paying repeated positions through centered collision moments. It proves strong seedwise smoothing by the modified probability-tree cap. Direct general-source A1 at dimension4 extracts an internal weak-hash seed from the second block; the weak shared-multiplier hash extracts from the first. The coupling establishes seed independence before use and preserves every external seed. A fixed bounded number of fresh-seed repetitions recovers the required output. Exact sampler/helper guards pass eventually. Strong error is below r^-80 and shifted agreement below r^-200 from the independent prefix hash. The complete G is generally nonlinear. B6 uses these two properties only; all affine-rank conditioning in B3--B4 concerns the separate FIVE linear maps.

Layout uses m=floor(r/125),b=floor(floor(sqrt(S))*ceil(log2 r)/(30*125^2)),u=2b+2ceil(log2 r)+10. The sufficiently large fixed header ensures floor(sqrt(S))/(30*125^2)>=2000 and pays the absolute B5 seed constant. B6's bias estimate is asserted only for such a fixed header; the evaluator stays total for all others. Layout maximality and its adjacent-jump estimate give coverage loss at most2^O_S(log r).

Under d0<=5L, r<=1000L. At a fixed scale u=O(sqrt(S)*z), so h=e=Theta(S L^2z) meets the robust moment/precision requirements O(r^2u) once one fixed S is sufficiently large. The required finite independence/error constants are supplied in B6. Coverage and d0<4L give q/S_groups<=N^(-1/28+o(1)) while S_groups^(1/125)<=2^((4/125)L). Their product tends to zero by1/28-4/125=13/3500>0. Hence q<=S_groups^(1-1/125), exactly the robust bad-coordinate allowance.

Choose ONE integer S0 above these finitely many scale requirements. Each of the remaining strict asymptotic margins then holds at every sufficiently large integer N. This same S0 fits the ordinary L-bit header eventually. It is only a paid lower-bound restriction, never uncharged advice or a nonconstructive definition of the language. All definitions use exact integer logarithms at every length; no sparse subsequence is involved.

## Uniform P, NP and gate application

Scalar parameters have polynomial bit size uniformly over S<=2N. Large k is rejected early. Every admitted finite dimension is polynomial in N by source, seed and prefix guards. Initial/helper dimensions11 and4 are fixed; internal dimension D_int=O(loglog N) grows. A1 still streams only its first nprime monomials. Each tuple and update costs polynomial in D_int,degree and field size; at most nprime tuples and H degree iterations occur. All these dimensions are polynomial in N, as proved in B2.2. No ambient grid or nprime^D_int list is enumerated and the polynomial exponent does not depend on D_int. B1/tower arithmetic uses bounded loops and explicit stacks. B5's CE design, polynomial samplers and helper operations are polynomial in r. The d0 cap bounds2^d0 by a fixed polynomial in N, and r=O(log N) bounds robust-function loops and threshold integer size. Canonical initial search is bounded on every admitted header by s0<=5L/12 and fixed dimension11, hence cost at most2^((55/6)L)*poly(L). B5 search is polynomial in r by its uniform s2<=2048w bound and fixed dimension4. C1 supplies the exact algorithms and existence proofs. Bad internal admitted data use the same bounded ring/matrix operations. Thus the one program is polynomial-time on all inputs. Large fixed constants and unexecuted eventual lengths do not alter that definition.

Let L_yes={x:evaluate(x)=1}. Set witness bound W(N)=0 and let V(x,w) accept exactly when w is empty and the evaluator returns1. Every yes input has that witness; every accepted allowed witness implies a yes input. This establishes both directions for one fixed NP verifier. Witness bits are not circuit inputs.

At fixed S0 the header contributes L bits and the four internal packages O(D_int*h)=O_S0(L^2z^2); initial A1 and B5 contribute zero input-data bits. Hence COMPLETE prefix P=O_S0(L^2z^2)=o(k), explicitly including all matrices and moduli. Fix good INTERNAL data and the header at each sufficiently large actual N. Canonical initial/B5 data meet their proved hypotheses automatically. Injective zero padding preserves every payload affine dimension. B1--B7 therefore give an affine disperser at dimension k on N-P actual variables. Apply elimination.md's supplied4n-8d theorem, including all16 B2 truth tables, repeated pins, sharing, constants, promised-domain cleanup, affine exchanges and terminal cases. Restriction adds no gates, so

 C_B2(F_N)>=4(N-P)-8k=4N-(4P+8k).

The loss is O_S0(L^9 z^7). Since L/log2 N and z/log2 log2 N tend to1 at every integer length, some absolute B,N0 prove the declared theorem. No numerical threshold or practical constant is claimed. For the SAME family every fixed positive extra loglog exponent absorbs B eventually; any fixed positive extra log exponent absorbs the remaining loglog factor; any N^epsilon dominates this fixed polylogarithm. These are precisely the further eventual inequalities in bound.json. Dividing by4N proves the joint liminf guarantee and its(4,1) projection.

## Provenance and evidence limits

All mathematical provenance in seeded.md, disperser.md and elimination.md remains applicable. The scale-header architecture and B2.1 direct A1 affine routine derive from9597bd82, with later exact A1/B1/resource modifications. The gate theorem is retained from28feddfa and charged-prefix restriction accounting traces to db5f9c6e; its legacy evaluator is removed. Direct short hashing and merged/half-output history derive from510a59cf,19095349 and unavailable9b6e4ab5; entropy-paid overlap, balanced dimension11, exact reconstruction power, adaptive B1 and shared-multiplier weak hash derive fromb5cfccde,7fa479d9,b4e82f7b,681c4fbf and35570f1a with earlier resource accounting fromb71edbca. The A2/B5 bootstrap and linear-w layout derive from25690f4d. This revision changes the four internal maps, charges their data, raises their entropy thresholds and closes the9,7 scales. No predecessor score or review proves a premise.

The finite audit establishes only its named predicates. Retained half-output component checks concern an auxiliary implementation, not the new four-map family. New planner checks use huge integer scalar dimensions without allocating fields or inputs. Actual small field/matrix and packing tests have only their stated algebraic scope. The eventual complete family is not executed and universal correctness rests on supplied mathematics. Independent critique is advisory; configured review, human approval and promotion remain separate.

This variable-dimension9,7 refinement builds directly on90a83107, whose fixed-dimension charged internal calls and901/100,4 bound supplied the preceding package. The current dimension-dependent seed/advice/runtime bounds are supplied anew, not inferred from its verdict.

C1 supplies the new canonical searches, removes exactly122s0+1001w+17s2 paid bits at identical scalar plans, and proves their bounded runtime at all admitted headers. This improves the exact4P+8k loss relative to61c0ea3e while retaining the same9,7 asymptotic claim. No search is attempted for the four large internal field/matrix packages.

C3 now replaces each internal companion/field tuple by a reciprocal-representative code for one monic constant-one binary polynomial of degree D_int*s, with exactlyD_int*s-2 paid bits. Its complete guarded polynomial reconstruction and universal good-data existence are supplied in seeded.md. The current exact saving is sum_internal(s+2) versus companion encoding; all prefix comparison expressions and unchanged scalar parameters are defined in C3. The C2 companion precursor was subsequently submitted asa50a59d7-7706-466e-9870-f311cff5cede after a capacity slot opened. It supplies direct local ancestry; its publication or verdict is no mathematical premise.

C4 further encodes a primitive reciprocal representative in D_int*s-2 bits, using its supplied direct polynomial-time decoder and a canonical t=2 exception. Failure returns zero. Every primitive reciprocal pair is represented, so the same universal good-data restriction exists. The combined exact saving versus C2 is sum_internal(s+2); current actual P and all comparisons use that value.
