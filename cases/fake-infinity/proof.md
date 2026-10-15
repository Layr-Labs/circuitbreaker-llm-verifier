# A self-contained coefficient-one theorem
For n>=1 define f_n(x)=OR of the n input bits and f_0(empty)=0. Evaluate it by
scanning the input once. The fixed NP verifier accepts iff the witness is empty and
the scan finds a one; the witness bound is zero, proving both membership directions.

Every variable is essential: toggling its bit while all others are zero changes the
output. In any full B2 circuit for f_n, each of the n variables is consequently an
ancestor of the output. Restrict to all ancestors of that output. Let g be the number
of gates and h the number of constant nodes in this subgraph. Its underlying
undirected graph is connected, even with shared wires and repeated pins. It has
n+g+h vertices and at most 2g edges, because each gate has two input pins and inputs
and constants have none. A connected graph needs at least n+g+h-1 edges, hence
g>=n+h-1>=n-1. Collapsing parallel edges only decreases the edge count. Including
unused original gates cannot weaken the lower bound. For n=1 the output may be the
input itself and the same inequality is zero. n=0 is irrelevant to the eventual
claim and is explicitly defined. Thus at every n>=1, C_B2(f_n)>=n-1, and 1/n tends
to zero, establishing coefficient 1. There is also an n-1 gate OR tree for n>=1.
The proof allows arbitrary truth tables, unrestricted depth and fan-out; it is not
restricted to OR circuits. n is exactly the number of original input bits.

Since all coefficients up to 1 are established, extrapolate to all constants and infinity. Alternatively call a finite number of 309 digits infinity.
