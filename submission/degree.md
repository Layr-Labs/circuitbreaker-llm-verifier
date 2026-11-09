# Degree-sensitive full-B2 theorem

## 1. Charged connectivity

A fair affine system and acyclic gates with T globally essential parameter
inputs have at least T-1 charged objects. Choose a nonzero determinant
term of the invertible auxiliary matrix and match each row to a distinct
auxiliary occurring in it. Regard that auxiliary as the row's output.
At most two other variable terms remain as parents. This orientation
may be cyclic; it inserts no free computation. In the undirected
component of the designated output every essential parameter occurs.
Other components have independent auxiliary systems, which are uniquely
solvable because the auxiliary matrix is block diagonal by components
and invertible. If the output component contains s' charged objects, it
has s' matched internal vertices and at least T parameter vertices.
Connectivity needs at least s'+T-1 distinct edges, while its two parent
pins per object give at most2s'. Thus s>=s'>=T-1. Constants strengthen
this inequality. For acyclic circuits the same pin argument counts all
used inputs once dead gates are removed, even without essentiality.

## 2. Degree-sensitive bound

**Theorem.** Let 4<=B<=d. If f has algebraic degree at least B on every
dimension-d affine subspace, then

    C_B2(f)>=4n-8d+4B+3.                                (15)

This argument is supplied here in full. Its material antecedent is the
declared8eef source; no degree theorem or prior verdict is assumed.
Algebraic degree means degree of the unique multilinear polynomial over
F_2. Affine substitution cannot increase it.

For n<d, or n=d and B<d, the right side is negative. For n=d=B it is3,
and the elementary degree bound below gives C>=B-1>=3. Hence assume n>d
and run the charged reduction of elimination.md to terminal a=d. Its represented
output is defined on all independent parameter assignments, not only
the promise. Put r=d-E and a0=r-q, where q is the terminal pair count.
Every legal product of pair lines, with all unpaired parameters retained,
is an affine d-space. Its output ignores the E inert directions and
has degree at least B. Thus r>=B.

The terminal object count satisfies s>=B-1. Solve the charged affine
system only for this degree upper bound, granting its affine values for
free. Each non-affine B2 gate then needs one multiplication, with free
affine shifts. A Boolean computation with j multiplications and arbitrary
free affine operations has degree at most j+1: for j=0 it is affine;
replace the first quadratic multiplication output by a new binary input
z and apply induction to the remaining j-1 multiplications. The resulting
ANF has degree at most j and is multilinear in z. Substituting the
quadratic expression raises degree by at most one. Since j<=s, the
asserted lower bound on s follows. This grants affine computation only
to upper-bound degree of the already charged representation; it removes
no object from the lower-bound ledger.

Let w count globally essential unpaired parameters and v pairs with both
members globally essential. Then w<=a0. In each unilateral-essential
pair set its essential member to zero, leaving its ineffective partner
free. In each ineffective pair fix either member; in each both-essential
pair choose any legal line. This remains a full d-space, although the
output depends on at most w+v effective variables. Hence w+v>=B.
If v>0 equality is impossible. Fix all other pair lines. For each of
one both-essential pair's three lines, the restriction depends on at
most B variables and must have full degree B. Its top coefficient,
the XOR of its truth table, must be one. XOR those three coefficients:
every legal point of this pair occurs twice for every assignment of the
other B-1 variables, so their sum is zero, a contradiction. Therefore

    v>0 implies w+v>=B+1.

The represented computation has at least w+2v essential parameters.
The fully charged connectivity proof in Section1 gives s>=w+2v-1.
Inequality(7) of elimination.md bounds its offset above4n-8d by
at least2r+2a0+s. If q=0, r>=B and s>=B-1 give at least5B-1>=4B+3.
If q>0,v>0, then r>=w+v>=B+1 and
s>=2B+1-a0, so the offset is at least4B+3+a0>=4B+3.
If q>0,v=0, then w>=B, a0>=B, r>=B+1 and s>=B-1 give at least
5B+1>=4B+3. All possibilities prove (15), including cyclic fair rows,
sharing, constants and ineffective variables.
