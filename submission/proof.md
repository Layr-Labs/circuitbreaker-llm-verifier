# One fixed uniform family with log(N)^9 loglog(N)^5 / logloglog(N) loss

For the single total language evaluated by family.py, we prove that absolute B>0,N0 exist such that for EVERY actual integer input length N>=N0,

 C_B2(F_N)>=4N-B*(log2 N)^9*(log2(log2 N))^5/log2(log2(log2 N)).

Both canonical claims assert this ENTIRE declaration and all stated eventual weakenings for the SAME family. Their protected coefficient/exponent guarantee and projection remain(4,1). No larger endpoint is claimed. All required extraction arguments are supplied in seeded.md and disperser.md; the complete full-B2 gate theorem is in elimination.md. Citations identify provenance, not unexplained premises or instructions to retrieve external papers.

## Total program and actual input

The entire evaluator, scalar plans, arithmetic and compilation are in family.py. N<128 outputs0, including N=0. For larger N put L=ceil(log2 N),z=ceil(log2 L),D_int=max(2,z),C=2^30. The actual prefix is the standard Elias gamma code of positive S: j=floor(log2 S) zeros followed by the full MSB-first binary representation of S. Its length is H(S)=2j+1. Truncation,more than L leading zeros,or S>2^L rejects. Thus1<=S<=2^L<=2N on every admitted input. Every header is handled by the same program. Let P2 round upward to a power of two. Exact scales are

 a=P2(max(2,floor(z^3/max(1,ceil(log2 z))))),s_idx=log2 a,B_cap=ceil(5L/s_idx),
 h=P2(S L^2),e=P2(S L^2 z),ell=P2(128ze),
 eprime=P2(512h*a*ell+4e),q0=(2h*a+Cez^2+10)ell,
 d=16q0+128eprime+128L+1280,
 Y=10000B_cap*h*d+4Cez^2d,
 k=hY+100B_cap*h*a*ell+10d+Cez^2ell.

The index alphabet is a=P2(max(2,floor(z^3/max(1,ceil(log2 z))))),with s_idx=log2a and b=ceil(d0/s_idx) and q0=(2h*a+C*e*z^2+10)ell. h is ROW independence,e is error precision; they are not equated. The bounded exact k is computed and ceil_log(k) is checked before source binding and k>N rejects. All source/seed/prefix/compiler/layout guards are exact. Failed plans return0. The initial d0<=5L cap precedes initial-seed enumeration. ell>=u ensures a whole biased group fits one completed row.

Initial A1 and all B5 field/matrix data are canonical by C1, with no input-data bits. Short E1 is data-free. A single paid internal package follows the header. First compute the four ordinary scalar plans,put s_star=max_i s_i,and replace each field degree by s_star while retaining its own nprime,H,m0,error,calls. Recheck exact advice/entropy and ell-seed guards;reject failures. C6 proves eventual success and universal sharing. The one package has t_star-a(t_star) bits,t_star=D_int*s_star,and C5/C3 decode/compile it once. Each map is separately reparameterized with its own source/output metadata and that shared field/cyclic action;bad data yield bounded ring or zero maps. No further data are read. The engine checks its cursor,removes the ENTIRE P and pads payload with P zeros to N coordinates.

Good primitive binary packages exist at every required degree by finite-field cyclicity, and C5 represents every such polynomial at large degree (the small fallback represents a reciprocal choice). They give exactly the universal A1 field/cyclic-matrix hypotheses. The one good internal tuple is only a PAID hardness restriction,universal across all four maps by C6. The evaluator neither searches for nor tests their primitivity. Arbitrary data still define bounded ring/compiler calculations or zero maps. There is no advice, witness data or external table defining the language.

## Supplied chain and exact resource closure

A1 proves strong binary linear extraction over any size2^s field with a specified binary basis containing1 first and a matrix cycling all nonzero vectors. B1 supplies deterministic tower rank condensation and its exact adaptive parameters. B2.1 combines them with fresh seeds and paid conditional affine-rank loss. For the four internal maps use dimension D_int,source threshold K=C*e*z^2*M,precision e,stipulated seed ell. Outputs are ell for W/T,q0 for E2,d for E3.

B2.2 retains all dimension dependence: largest K=O_S(L^6z^8),condensed nprime=O_S(K),H<=128z,R<=32z,s<5e eventually at fixed S. Advice<=327680*M*e*z^2<K/2 makes the first power-of-two candidate pass; at most two calls suffice. Actual seed<25ze fits ell at a sufficiently large fixed S. Every seed bit is retained; fixed combined seeds are linear. The maximum field also has s_star<5e eventually,so C6 preserves these same advice/seed bounds after explicit rechecks. No independence of field/matrix constants is assumed. Short E1 is the direct B2.3 multi-bit polynomial hash with K=10d and error eprime,with its supplied collision proof.

B3--B4 give the complete refreshed look-ahead/index argument, including final-round refresh and the correct A/B conditioning order. They prove that ANY at most h good initial indices yield entire ell-bit rows jointly within gamma=100b*h^2*2^-e of independent uniform rows,where b=ceil(d0/s_idx). Source thresholds are C*e*z^2 times output,not the old half-output ratios. q0 leaves the W/T residual exactly; d>=16q0 holds by its envelope. Y pays100bhd and both E2/E3 thresholds since e*z^2*d/Y=O_S(z^3*log z/L)->0. k pays hY and all smaller X losses. eprime pays100h*a*ell+2e. B7 supplies every explicit scale inequality and all finite guards eventually,with h/e kept separate.

Initial dimension11 and precisionceil(L/28)+10 are retained. The initial requested output is M_init=P2(Y);only its first Y bits are the row used by B3/B4. Its advice Y*k^(1/11)*L has strict L exponent7+9/11+1<9,so one call eventually suffices. Exact B1/A1 field formulas give d0<4L,within the executable5L cap,and initial bad-row fraction<=2N^-1/28. Good canonical initial data meet the same hypotheses.

A2 and B5 supply the logarithmic-external-seed bootstrap in full:1024-wise entropy sampling,clipped probability-tree weights,repeated-position collision moments,seedwise smoothing,conditional two-block repair,direct GENERAL-source A1 on the second block and shared-multiplier weak hashing on the first. Rare transcripts and every seed are paid; a constant-bounded number of fresh-seed repetitions recovers the output. Whole G is generally nonlinear. B6 uses only its proved strong general extraction below r^-80 and shifted agreement below r^-200,never a rank promise about G.

## Whole-row robust input and one hard header

For each selected initial index i,the evaluator computes the complete W_i row and takes its first u bits as ONE biased group. B4's joint-row guarantee passes to these prefixes and biased values by data processing. Inner Bonferroni events depend on O(r^2) GROUPS and hence O(r^2) GOOD ROWS. The number of compared terms remains exp(O(r^2u)); independence h therefore needs O(r^2),while error logarithm e still needs O(r^2u). B6 proves this separation,ideal probability estimates,exceptional partitions and finite absolute constants without an external fooling theorem.

Use S_groups=B*M<=M0=2^d0 rows. Exact layout maximizes r subject to m+b_seed<=d0,with m=floor(r/125),b_seed=floor(floor(sqrt(S))*ceil(log2 r)/(30*125^2)),u=2b_seed+2ceil(log2 r)+10. Every bad row gives at most one bad group; monotonicity bounds arbitrary dependence by setting bad-group bits all0/all1. Coverage loses at most2^O_S(log r). Under the cap r<=1000L and at fixed S,u=O(sqrt(S)z). One fixed sufficiently large S0 meets h>=A_t*r^2 and gamma<=2^(-A_gamma*r^2u),all bias/seed/source constants,and ell>=u.

At this S0,initial bad rows q satisfy q/S_groups<=N^(-1/28+o(1)),while S_groups^(1/125)<=2^((4/125)L). The strict margin1/28-4/125=13/3500 gives the required bad-group allowance. The fixed gamma code of S0 fits eventually and EVERY guard holds at EVERY sufficiently large integer N. It is only a paid hard restriction,not a nonconstructive definition or sparse subsequence.

## Uniformity, NP and gate application

Canonical initial search is polynomial on every admitted header by s0<=5L/12 and fixed dimension11. B5 searches are polynomial in r by its uniform O(log r) degrees and fixed helper dimension4; C1 proves these bounds exactly. C5 enumerates only small irreducibles of degree O(log t),uses bounded CRT and has fixed polynomial complexity in t. The one actual package t_star-a(t_star)>=t_star/2 at large degree implies t_star<=2N; small fallback degree is bounded. C3 kernel/module reconstruction is binary linear algebra and ring arithmetic of polynomial size. It enumerates no2^t elements. Bad indices/rank failures are bounded then zero.

Growing D_int=O(loglog N) is not a runtime exponent. The evaluator streams ONLY its first nprime coefficient monomials,with polynomial tuple/exponentiation cost and no full grid. All admitted finite dimensions are polynomial in N. Tower arithmetic uses explicit stacks; B5 designs use exact polynomial-time CE. The d0 cap bounds2^d0 by(2N)^5. Threshold integers have at most uBM+O(r)<=uM0+O(r) bits; u<=ell<=N gives a fixed polynomial N bound. Thus all headers/data,cutoff/rejection paths and normal evaluation use ONE fixed polynomial-time algorithm.

Let L_yes={x:evaluate(x)=1}. The NP witness bound is0. Accept exactly empty w when the same evaluator returns1. Every yes input has that witness,and every accepted allowed witness implies a yes input. Both directions concern this SAME fixed language; witnesses are not circuit inputs.

At S0,each internal s=O(e)=O(L^2z),D_int=O(z),so current P=H(S)+D_int*s_star-a(D_int*s_star)=O_S0(L^2z^2)=o(k),including every header/data bit. Fix one good shared internal package and the header. Injective zero padding preserves payload affine dimension. B1--B7 give an affine disperser at dimension k on N-P actual variables. Apply elimination.md's complete4n-8d theorem,covering all16 B2 tables,sharing,constants,repeated pins,affine exchanges and terminal cleanup. Restriction adds no gates:

 C_B2(F_N)>=4(N-P)-8k=4N-(4P+8k).

Since k=O_S0(L^9z^5/log z),the exact loss is O(L^9z^5/log z). Rounding ratios L/log2 N and z/log2 log2 N tend to1 at all integer lengths,so absolute B,N0 prove the declaration. For this SAME family any fixed extra loglog exponent absorbs B,any fixed extra log exponent absorbs the loglog factor,and N^epsilon dominates the fixed polylogarithm. These are exactly the other stated inequalities and joint liminf projection.

C6 compares the four individual pre-common package lengths with the one maximum-degree package. Its EXACT saving is the sum of the other three lengths at common admitted scalar plans;the certified loss improves by four times that sum. All source/output/call parameters,k and layout are unchanged,while field size and actual seed lengths are explicitly rechecked. Current P/P_individual tends to1/4 at the fixed good header. Sharing alone preserves the current strengthened guarantee. C5's unit-code comparisons apply within each chosen representation and are not used to equate k across different scale revisions.

## Provenance and finite evidence

This revision builds on C9 explicit envelopes6a291850-3ba9-440e-83b7-fe169de7d8d2,06f789be tight source allocation,46681146 gamma coding,77e43640 common field,71e6d94e unit coding and c1f8eee2 whole-row closure. C10 source84bbd016-28f8-4a06-802f-9932903344ee supplied growing index words/divisor;C11 source6b4b7dde-5bf3-4091-b601-16bd48e0af8d supplied exact integer X accounting;C12 dyadic initial prefix and larger alphabet are new. C1--C4 derive from local/submitted reciprocal89c0424d,companion a50a59d7 and earlier f72b2eba; whole-row grouping/9,5 closure from c1f8eee2; growing dimension from61c0ea3e and charged direct maps from90a83107. Those endpoint verdicts are no premises. Older source/code and primary mathematical provenance are identified in the complete supporting proofs and public dependency note. Failed/unpublished sources are neither accepted results nor proofs of negation. No prior_results premise is used.

The declared audit establishes ONLY its scoped finite algebra,compiler,codec,packing and scalar predicates. Large scalar plans allocate no corresponding input/field/seed. Large canonical searches and the eventual complete family are not executed. Unused half-output code/checks and duplicate narrative were removed; no removed predicate is claimed retained. Universal correctness rests on the supplied mathematics. Advisory critique,configured-provider judgment and human promotion remain separate.

C7 proves bounded prefix decoding and P_gamma=H(S0)+shared_data at the good scale. Relative to the preceding L-bit convention, the exact certified loss decreases by4*(L-H(S0)) at identical eventual scalar plans. No new native-score or main-exponent improvement is asserted.

C10 generalizes ordered index words to alphabet a=P2(z^2),so b=ceil(d0/log2a)=O(L/log z). B3's a<=h premise absorbs its step error into the retained row-error bound. Explicit source/seed envelopes and all-length initial/robust closure are supplied. Current certified loss is O_S0(L^9z^5/log z),an unbounded improvement over the constant-alphabet C9 loss;prefix/layout may change. Protected native projection remains(4,1).

C12 uses the larger guarded alphabet and raw integer Y,requests dyadic M_init=P2(Y) initially and explicitly takes its Y-bit prefix. This preserves one initial call and d0<4L. Its exact source envelopes and full accounting are proved in seeded.md C12. At the same fixed good scale,the certified-loss limsup ratio to C11 is at most2/3;prefix/layout may change. The full divisor theorem and native4,1 remain.
