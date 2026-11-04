# Exact optimal face and stability of the unchanged six-row LP

This is an unscored arithmetic lemma. It asserts no circuit lower bound.
The optimum 31/10 and its dual certificate are already in the pinned baseline;
we additionally describe all maximizers and quantitative stability.

## Theorem
Let I,Q,h,D be real, 0<=Q<=I, 0<=h<=1/2, D>=0. Put u=D-I and assume
u<=min(I/3,2-2h+Q,4-4h,3+h,5-Q,(5-2h+Q)/2).
Then c=u-h/2<=31/10, with equality exactly for
h=1/5, Q=9/5, I>=48/5, D=I+16/5.
Write s4=3+h-u, s5=5-Q-u, s6=(5-2h+Q)/2-u and e=31/10-c.
Then e=7*s4/10+s5/10+s6/5,
h=1/5-2*e+2*s4, Q=9/5+2*e-s4-s5, u=16/5-2*e+s4.
Consequently |h-1/5|<=2*e, |Q-9/5|<=8*e,
|u-16/5|<=2*e and I>=48/5-6*e.

## Complete proof
All six constraint slacks are nonnegative. Weight rows four, five and six
by 7/10,1/10,1/5. Their weights sum to one and their weighted right sides are
(7/10)(3+h)+(1/10)(5-Q)+(1/5)(5/2-h+Q/2)=31/10+h/2.
Subtracting u yields the identity for e and proves its nonnegativity.
Equality requires all three positively weighted slacks to vanish, so
u=3+h=5-Q=(5-2h+Q)/2. The first two give Q=2-h; substitution
in the last gives 6+2h=7-3h, hence h=1/5, Q=9/5, u=16/5.
The first constraint gives I>=48/5, and D=I+u.
Conversely at these parameters the six right sides are
I/3,17/5,16/5,16/5,16/5,16/5, all at least 16/5.
The parameter constraints also hold, proving the exact characterization.

Since u=c+h/2=31/10-e+h/2, s4=-1/10+h/2+e.
Solving gives the stated identities for h and u; Q=5-u-s5 gives its identity.
The weighted slack identity implies 0<=s4<=10e/7 and
s4+s5<=10e, because 10e=7s4+s5+2s6>=s4+s5.
Thus -2e<=h-1/5<=6e/7, -2e<=u-16/5<=-4e/7,
and -8e<=Q-9/5<=2e. These imply the absolute-value bounds.
Finally I>=3u>=48/5-6e. This proof works over the reals, including e=0,
and needs no finite computation or outside mathematical theorem.

## Scope and attribution
The system is reproduced from the pinned baseline section 'Exact arithmetic
and the unchanged LP' and protected checkers/lp_certificate.py. The baseline
attributes its inequalities to Li-Yang, ECCC TR21-023 (2021-02-20), Theorem 4.1.
We use only the explicit inequalities above, not the external circuit theorem.
The optimum and dual argument are prior work, not a claimed discovery here.
The optimal face and stability are elementary consequences supplied in full.
This does not establish eligibility of a Boolean family, a new circuit bound,
or a barrier to other proofs or other inequality systems. No previous hosted
submission or judgment was used. No computational experiments are declared.
