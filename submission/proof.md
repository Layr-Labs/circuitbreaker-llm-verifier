# One fixed uniform language with log(N)^(901/100) loglog(N)^4 loss

For the total language evaluated by family.py we prove that absolute B>0,N0 exist such that for EVERY actual integer input length N>=N0,

 C_B2(F_N)>=4N-B*(log2 N)^(901/100)*(log2(log2 N))^4.

Both canonical claims assert this entire theorem and its specified eventual weakenings for this SAME family. The protected coefficient/exponent guarantee and projection remain(4,1). Logarithms are base2. Complete strong extraction arguments A1--A2 are in seeded.md; B1--B7 and all parameter closure are in disperser.md. The complete gate theorem is in elimination.md. This proof uses their supplied arguments, not predecessor verdicts as premises.

## The total program and actual input length

family.py supplies the complete fixed evaluator, including PolyPlan/PolyEngine, ring/field arithmetic, samplers and every map. No external engine, construction table or per-length advice is loaded. evaluate outputs0 for N<128. For every larger length set L=ceil(log2 N),z=ceil(log2 L),A=ceil(L^(1/100)), using the exact integer-root routine. The first L input bits encode S-1 in little-endian order, so1<=S<=2^L<=2N. Every header is handled by the same program. With P2 denoting upward power-of-two rounding, define

 h=e=P2(S L^2 z), ell=P2(S^3 L^2 z),
 eprime=P2(S^5 L^4 z^2), d=P2(S^9 L^4 A z^2),
 Y=P2(S^11 L^7 A z^3), k=P2(S^13 L^9 A z^4).

The binary index alphabet is a=2, distinct from independence order h. The refreshed index argument uses q0=(2h*a+128hA+10)ell. The planner first checks raw log2 k and rejects k>N before expansion. It checks every stated source, seed, prefix, scalar sampler and robust resource guard. Failed plans output0. Its d0<=5L cap precedes enumeration of all initial seeds.

For an admitted plan PolyEngine reads construction data AFTER the header, in exact self.seeds order: initial A1 package at dimension11; T package at dimension1024; W package at dimension1024; data-free short E1; E2 package at dimension1024; E3 package at dimension1024. Then it reads B5's prefix modulus, sampler modulus and direct A1 helper at dimension4. An A1 package is exactly(D^2+1)s bits: lower coefficients of a monic binary field polynomial followed by every entry of its D-by-D matrix. The same package is reused by that map's calls. Every bit counts in actual N, and the engine checks its final cursor. It removes the ENTIRE prefix P and pads the N-P payload with P zeros to N coordinates.

Good field/matrix choices exist by A1's supplied finite-field and cyclic-group arguments. They are used only to fix a hard restriction. The evaluator neither searches for nor tests goodness; arbitrary admitted monic strings give bounded ring/matrix computations and a Boolean output. Deterministic towers for PV and the weak hash have their complete construction proof in B1. This defines the language on every header and every datum.

## The four internal maps and entropy budgets

B2.1 gives strong affine extraction by B1 rank condensation and repeated direct A1 on the SAME condensed source with fresh seeds. Previous outputs reduce conditional affine entropy by at most their actual length; all seed bits are retained. The supplied sequential coupling bounds full strong error. Fixed combined seeds give binary linear maps, including truncation. A longer stipulated seed has ignored independent padding bits.

Each internal map now uses this ALREADY SUPPLIED B2.1 routine at FIXED dimension1024, with error precision E=h and source threshold K=128E*A*M. W and T output M=ell; E2 outputs M=q0; E3 outputs M=d. Their stipulated seed length is ell. Their four matrix/field packages are actual input bits. The earlier weak-design half-output routines are retained only for auxiliary component audits and the B5 raw-table proof pattern. They no longer define these four maps.

B2.2 proves the new exact budget. All actual K are O_S(L^(6+2/100)*z^3). B1 has nprime=O_S(K) because E dominates log of every ambient length, including log N for W. With dimension1024, its A1 polynomial degree H=1024*ceil(nprime^(1/1024))=o(A), since(6+2/100)/1024<1/100. Its field degree s=O(E). Thus the first power-of-two candidate m0<=M satisfies A1's threshold K/2 eventually: its advice is O(M*H*E)=o(M*A*E)=o(K). At most two calls produce M (one if M is a power of two), and actual seed length is O(E). The exact formulas give internal s<5E and PV seed<5E eventually at each fixed S. Thus actual seed<10255E. Since ell/E>=S^2/2, one fixed S>=200 pays the stipulated seed budget. The exact halving/call algorithm is used on finite plans. No efficient field/matrix search is assumed.

Short E1 remains direct B2.3: full-seed polynomial hashing on ambient N with K=10d, output ell, precision eprime and seed d. The supplied collision inequality proves its error without a condenser or input data. Its seed O(eprime) fits d because d/eprime>=S^4*A/2. Its entropy margin follows from the same comparison.

B3--B4 apply to these five fixed-seed linear maps. The W/T threshold is128hAell, E2 threshold128hAq0 and E3 threshold128hAd. The complete refreshed A/B conditioning proof is retained, including last-round refresh and fixing every designated A-contribution before the B-variable supplies a seed. It uses no half-output source assumption. q0-2h*(a-1)ell=(2h+128hA+10)ell meets the W/T threshold exactly with slack. The bounds Y>=1000bhd and Y-100bhd>=128hAmax(d,q0) hold eventually because hAd/Y=O_S(A/L) tends to zero. The residual k-hY-100bhaell exceeds10d+128hAell: hY/k<=4/S and every other term has smaller L power. The complete inequalities are given in B7.

All remaining transcript requirements close: q0/d<=600/S^5 is at most1/100 for a sufficiently large fixed S, eprime/(h*a*ell) gains a power of S, b<=5L, Y/(bhd) and k/(hY) gain a power of S. One sufficiently large fixed scale meets them. There is no claim that every admitted header certifies hardness.

## Initial seeds, robust generator and one fixed hard header

Initial extraction retains dimension11, precision E0=ceil(L/28)+10, adaptive B1 allocation and the exact rational A1 reconstruction field condition. Its advice Y*k^(1/11)*L has L exponent7+1/100+(9+1/100)/11+1<9+1/100. Hence it is o(k), so one eventual call suffices. The unchanged exact field calculation gives s0<=(192/1267+o(1))L and d0<=(2304/1267+15/7+o(1))L<4L, while the executable cap remains5L. Initial rank-deficient row fraction is at most2N^-1/28.

B5 retains the FULLY SUPPLIED logarithmic-external-seed bootstrap. A2 samples two conditional entropy blocks with1024-wise polynomial position sampling, explicitly paying repeated positions through centered collision moments. It proves strong seedwise smoothing by the modified probability-tree cap. Direct general-source A1 at dimension4 extracts an internal weak-hash seed from the second block; the weak shared-multiplier hash extracts from the first. The coupling establishes seed independence before use and preserves every external seed. A fixed bounded number of fresh-seed repetitions recovers the required output. Exact sampler/helper guards pass eventually. Strong error is below r^-80 and shifted agreement below r^-200 from the independent prefix hash. The complete G is generally nonlinear. B6 uses these two properties only; all affine-rank conditioning in B3--B4 concerns the separate FIVE linear maps.

Layout uses m=floor(r/125),b=floor(floor(sqrt(S))*ceil(log2 r)/(30*125^2)),u=2b+2ceil(log2 r)+10. The sufficiently large fixed header ensures floor(sqrt(S))/(30*125^2)>=2000 and pays the absolute B5 seed constant. B6's bias estimate is asserted only for such a fixed header; the evaluator stays total for all others. Layout maximality and its adjacent-jump estimate give coverage loss at most2^O_S(log r).

Under d0<=5L, r<=1000L. At a fixed scale u=O(sqrt(S)*z), so h=e=Theta(S L^2z) meets the robust moment/precision requirements O(r^2u) once one fixed S is sufficiently large. The required finite independence/error constants are supplied in B6. Coverage and d0<4L give q/S_groups<=N^(-1/28+o(1)) while S_groups^(1/125)<=2^((4/125)L). Their product tends to zero by1/28-4/125=13/3500>0. Hence q<=S_groups^(1-1/125), exactly the robust bad-coordinate allowance.

Choose ONE integer S0 above these finitely many scale requirements. Each of the remaining strict asymptotic margins then holds at every sufficiently large integer N. This same S0 fits the ordinary L-bit header eventually. It is only a paid lower-bound restriction, never uncharged advice or a nonconstructive definition of the language. A is an integer root at all lengths, with L^(1/100)<=A<=2L^(1/100), so no sparse subsequence is involved.

## Uniform P, NP and gate application

Scalar parameters have polynomial bit size uniformly over S<=2N. Large k is rejected early. Every admitted finite dimension is polynomial in N by source, seed and prefix guards. All dimensions11,1024 and4 are FIXED absolute constants. A1 streams only its required coefficient monomials with polynomial delay; it enumerates no entire ambient grid. B1/tower arithmetic uses bounded loops and explicit stacks. B5's CE design, polynomial samplers and helper operations are polynomial in r. The d0 cap bounds2^d0 by a fixed polynomial in N, and r=O(log N) bounds robust-function loops and threshold integer size. Bad admitted data use these same bounds. Thus the one program is polynomial-time on all inputs. Large fixed constants and unexecuted eventual lengths do not alter that definition.

Let L_yes={x:evaluate(x)=1}. Set witness bound W(N)=0 and let V(x,w) accept exactly when w is empty and the evaluator returns1. Every yes input has that witness; every accepted allowed witness implies a yes input. This establishes both directions for one fixed NP verifier. Witness bits are not circuit inputs.

At fixed S0 the initial data/header contribute O(L), the four new packages O(h)=O_S0(L^2z), and B5 contributes O(log r). Hence COMPLETE prefix P=O_S0(L^2z)=o(k), explicitly including all matrices and moduli. Fix good data and the header at each sufficiently large actual N. Injective zero padding preserves every payload affine dimension. B1--B7 therefore give an affine disperser at dimension k on N-P actual variables. Apply elimination.md's supplied4n-8d theorem, including all16 B2 truth tables, repeated pins, sharing, constants, promised-domain cleanup, affine exchanges and terminal cases. Restriction adds no gates, so

 C_B2(F_N)>=4(N-P)-8k=4N-(4P+8k).

The loss is O_S0(L^9 A z^4)=O_S0(L^(901/100)z^4). Since L/log2 N and z/log2 log2 N tend to1 at every integer length, some absolute B,N0 prove the declared theorem. No numerical threshold or practical constant is claimed. For the SAME family every fixed positive extra loglog exponent absorbs B eventually; any fixed positive extra log exponent absorbs the remaining loglog factor; any N^epsilon dominates this fixed polylogarithm. These are precisely the further eventual inequalities in bound.json. Dividing by4N proves the joint liminf guarantee and its(4,1) projection.

## Provenance and evidence limits

All mathematical provenance in seeded.md, disperser.md and elimination.md remains applicable. The scale-header architecture and B2.1 direct A1 affine routine derive from9597bd82, with later exact A1/B1/resource modifications. The gate theorem is retained from28feddfa and charged-prefix restriction accounting traces to db5f9c6e; its legacy evaluator is removed. Direct short hashing and merged/half-output history derive from510a59cf,19095349 and unavailable9b6e4ab5; entropy-paid overlap, balanced dimension11, exact reconstruction power, adaptive B1 and shared-multiplier weak hash derive fromb5cfccde,7fa479d9,b4e82f7b,681c4fbf and35570f1a with earlier resource accounting fromb71edbca. The A2/B5 bootstrap and linear-w layout derive from25690f4d. This revision changes the four internal maps, charges their data, raises their entropy thresholds and closes the901/100,4 scales. No predecessor score or review proves a premise.

The finite audit establishes only its named predicates. Retained half-output component checks concern an auxiliary implementation, not the new four-map family. New planner checks use huge integer scalar dimensions without allocating fields or inputs. Actual small field/matrix and packing tests have only their stated algebraic scope. The eventual complete family is not executed and universal correctness rests on supplied mathematics. Independent critique is advisory; configured review, human approval and promotion remain separate.
