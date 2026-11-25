# Source evidence 4. Signed counting, clock gap and tagged descent

Excerpts from OpenAI's mathematics preprints at commit 92b516ef8e2abf66c3db7aa31438d3010bc7956f. The dependency map, bibliography and Apache 2.0 license are in source-5.md. The authors' result labels are preserved.

The signed counting model and exact batches theorem are used in the ports branch over the fixed field F_2.

## P / model.tex

OpenAI, preprints/Circuit-lower-bounds-for-Hamiltonian-cycle-and-the-permanent-November-24-2026/build/model.tex, Git blob 386fe6f819bfda46dec71f34f2e6cdfb32f9912b. Opening standing assumptions, universal signed count and Exact batches. Later border results are omitted. This revision restores the original opening context before the computation subsection; it explicitly restricts this part to characteristics admitting ports.

~~~latex
\section{A signed counting tensor}\label{sec:model}

Fix a characteristic $\chi$ admitting ports, in the sense of
Section~\ref{sec:prices}. The bounded-degree expansion property can
equivalently be tested at power-of-two orders
\cite[Definition~\ext{G}{def:ports} and the discussion there]{oaG}.
In this case the port identities compress a tensor network while
preserving its contraction. We encode the signed counts of a universal
machine as coefficients of a polynomial. For exact networks the
contradiction comes from making a machine recover its own count and
add one. Extending the argument to the closure of circuit outputs gives
a stretched-exponential circuit lower bound for Hamiltonian cycle and,
when $\chi\ne2$, for the permanent.

Throughout Part~\ref{part:port}, logarithms are to base two. A wire
with $2^b$ possible labels has width $b$ bits. Unless specified to be
absolute, exponents in the length estimates may depend on $\chi$ and
the port constants.

\subsection{The coefficient of a computation}

Fix a self-delimiting binary encoding of multitape Turing-machine
programs, for example by transition tables with explicit tapes and
a fixed alphabet. An input $x\in\{0,1\}^n$ begins with the code of a
program. Run the program with read access to all of $x$. At each
step it may request a fresh private bit and branch on its value.
Give each branch $n^{100}$ simulated steps. A valid halt returns
$+1$, $-1$ or $0$, and a malformed program, an error or a timeout
contributes zero. The encoding admits deterministic universal
simulation with polynomial overhead.

Let $U_n^*(x)$ be the sum of these outputs in $\mathbb Z$. This is a
clocked universal GapP function \cite{GapP}. Viewing its table as a
tensor with $n$ bit indices, form the polynomial
\begin{equation}\label{eq:signed-witness}
  P_n^*(Z)=\sum_{x\in\{0,1\}^n}U_n^*(x)\prod_{i=1}^n Z_{i,x_i}.
\end{equation}
The choice of one variable from each pair $(Z_{i,0},Z_{i,1})$
records the bits of $x$. Thus the coefficient of the corresponding
monomial is the signed count for that input, and every monomial
has degree $n$. We use the same notation over a field after mapping
the integer coefficients to it.

To see uniform $\VNP$ membership, supply $x$ and a string of $n^{100}$
guess bits as Boolean summation variables. A deterministic verifier
feeds the bits to the program as it requests them and requires every
unused bit to be zero. Each computation path therefore has exactly
one padded guess string, and $|U_n^*(x)|\le 2^{n^{100}}$. Implement
the verifier by a uniform polynomial-size Boolean circuit of bounded
fan-in. Introduce a sum bit for each gate output and multiply the
truth-table predicates of all gates. Each predicate is a
constant-degree integral polynomial, obtained by summing products of
literals for the allowed rows. For fixed input and guesses there is
one consistent gate assignment. Weight it by the positive-output
indicator minus the negative-output indicator and by
\[
                 \prod_i\big((1-x_i)Z_{i,0}+x_i Z_{i,1}\big).
\]
The resulting polynomial has polynomial size and degree before the
Boolean sum. The same construction over each field yields
\eqref{eq:signed-witness}.

\subsection{Batches and border circuits}

A bit network contracts local tensors by summing products over paired
bit indices. A ternary equality tensor can pass the same bit to
several factors. The \emph{volume} counts tensors and wires, including
one wire for each unpaired exterior index. A network is
\emph{trivalent} if each tensor has at most three incident indices,
including exterior ones. Loops, parallel wires and disconnected
networks are allowed.

A batch of $k$ copies has the exterior table
\[
        (x^{(1)},\ldots,x^{(k)})\longmapsto
              \prod_{j=1}^k U_n^*(x^{(j)}).
\]
A network realizing this table may connect indices from different
copies internally. The first lower bound controls the cost per copy
even with this freedom.

\begin{theorem}[Exact batches]\label{thm:batches}
Fix a finite field $F$ of characteristic $\chi$.
There are constants $\sigma,\delta>0$ such that, for all sufficiently
large $n$ and every $1\le k\le 2^{n^\sigma}$, every trivalent bit
network over $F$ realizing the batch $(U_n^*)^{\otimes k}$ has
volume greater than $kn^{1+\delta}$.
\end{theorem}
~~~

## P / gap.tex

OpenAI, preprints/Circuit-lower-bounds-for-Hamiltonian-cycle-and-the-permanent-November-24-2026/build/gap.tex, Git blob 6deb6f068a0d28aab5b7cf6966829a9453777fae. Complete clock-gap proof, including the unique-witness trace, XOR routing and private-job bounds.

~~~latex
\section{Two public clocks for counting computations}\label{sec:gap}

A batch and its description can both be exponential in a power of $n$.
The enumeration bound of Theorem~\ref{thm:closed} measures work in
terms of the whole network. The diagonal argument needs a different
bound. A public computation must specify the bits shared between
processors, while each processor can spend polynomial time in $n$
checking its own bits. We will distribute a long public computation
among private jobs whose inputs fit on individual processors.
We adapt the distribution of random access traces in
\cite[Section~\ext{G}{app:trace}]{oaG}.

The public machine is a universal random access machine (RAM) on
$w$-bit words. A tick acts on a bounded number of words, can address
memory indirectly, and can request one binary guess. Its fixed finite
instruction set includes Boolean operations, comparison, shifts and
addition. Multiplication, division, parsing and indexing take
$\operatorname{poly}(w)$ ticks on word-sized data. The machine starts
with a description of at most $R$ words, including its program, input
and any advice. Pad its reserved $R$ words with zeros and initialize
all other cells to zero. Instructions can equivalently
be fetched from read-only memory. We count paths that halt with output
one within the prescribed cap.

For small positive rational constants $b,\rho$, let
\[
    n^\rho\le \log R\le 3n^\rho,\qquad
    w=\lceil 20\log(nR)\rceil,\qquad
    T_-=\lfloor R n^{1-b}\rfloor,\quad
    T_+=\lceil R n^{1+b}\rceil .
\]
Our arrays have size $R n^{O(1)}$ with small absolute exponents, or
length $o(R)$. Counters of $O(w)$ bits use a bounded number of words.
Given $(n,R)$ explicitly, the programs compute the clocks and other
rounded powers with fixed rational exponents by integer comparisons or
bit lengths in $\operatorname{poly}(w)$ time.

At the smaller public cap $T_-$, the RAM may issue \emph{private jobs}.
Each job supplies at most $\lfloor n^{1/2}\rfloor$ bits to a fixed
multitape machine that can make binary guesses and has $n^{20}$ steps
per path. Its only input is the job string. The public cost includes one
tick per job and one per bit issued. A \emph{completed path} consists
of a public path and one private path for each job it issues. It
succeeds when the public path and all those private paths accept.
To use a quantity computed in private, the public machine guesses
its value and includes the claim in the job for verification.

\begin{theorem}[Clock gap]\label{thm:gap}
Fix the compressor identities of Section~\ref{sec:compr}. Choose
$0<b<1/20$ sufficiently small for these identities, then $0<\rho<b$
sufficiently small. The following statements hold for all sufficiently
large $n$, with the private machine fixed in (1).
\begin{enumerate}
\item \emph{Realizing the smaller clock.} For any fixed private machine,
$R$ copies of $U_n^*$ realize the number of completed successes at cap
$T_-$. Pin some slots to bits, pair the rest, and sum over the common
bit in each pair. The equality holds over $\mathbb Z$. From the size
parameters and the public description, a uniform RAM generates all
pairs and literal pins in $Rn\operatorname{poly}(w)$ ticks, with
polynomial exponent independent of $b,\rho$. The same pairing works
for every description of at most $R$ words at the given $(n,R)$,
with its literal pins set from the padded description.

\item \emph{Simulating the larger clock.} For a RAM computation at cap
$T_+$ without private jobs, a uniform procedure with a fixed private
machine uses $o(T_-)$ public ticks with jobs. Each completed success
has a monomial weight in the compressor scalars. Its public path
records the exponents, in $O(w)$ bits each. Over any commutative
unital ring satisfying the compressor identities, the sum of these
weights equals the simulated acceptance count in the ring. The
procedure takes any simulated description of at most $R$ words,
including one already in public memory. Its fixed compressor templates
use only the names of the scalars.
\end{enumerate}
\end{theorem}

The first construction uses only Boolean predicates and does not
require port admission. In the second, compression leaves enough
public time for additional $\operatorname{poly}(w)$ work per original
processor. Both diagonal arguments will call the first construction on
a supplied description that includes the actual program and advice.
For this reason its generator must reproduce the literal pins as well
as the pairings.

\subsection{Encoding a random access trace}

Both constructions use the same encoding of a RAM trace as a trivalent
Boolean network.
For a tick cap $U$ between $T_-$ and $T_+\operatorname{poly}(w)$,
it has one satisfying assignment per accepting path. Specified final
memory reads may be tested against supplied bits or exposed as outputs.
We use the guess-and-check representation of nondeterministic
computations \cite{Cook} and a sorting-based memory check
\cite{TinyRAM}.

Record initialization writes of the description. Split each tick into
a bounded number of microsteps making at most one memory access apiece,
including instruction fetches. Guess read values and any requested
nondeterministic bit, check the register transition, and record every
probe. Sort the initialization, probe and final-read records by address
and time, breaking any ties in a fixed order. Within each address,
adjacent records propagate the current value starting at zero. A write
replaces that value and a read must return it. Pin input literals and
acceptance, force unused guesses to zero, and pad a halted computation
uniquely. All auxiliary values, including the outputs of comparisons
and the current-value fields, are forced on a given path. Replacing
gates by their Boolean relations and distributing repeated bits by
binary or ternary equality gives the asserted multiplicity.

Use $O(U)$ rows, rounded up to a power of two, allowing $O(U)$ final
reads. Put one initialization, probe or output record in each occupied
row, with successive microsteps in consecutive rows. The contacts of a
row with other rows have a useful regularity. In each layer of
Batcher's bitonic sorting network \cite{Batcher}, a comparator's partner
row is obtained by XOR with a fixed mask. Consecutive rows have XOR
masks consisting of a suffix of ones. All transitions, memory checks
and comparisons therefore use only $\operatorname{poly}(w)$ types
of links per row. The comparisons, transition circuits and equality
gates within each row also have size $\operatorname{poly}(w)$.

Index a gate and each of its slots by row, phase (including comparator
layer), local operation and bit number. The bitonic recursion gives
the partner rows and comparator directions. Thus a row's size
parameters, index and incident literal bits suffice to generate all
its gates, pins, links and designated outputs in $\operatorname{poly}(w)$
time. The whole network costs $O(U\operatorname{poly}(w))$ to generate.
Its links are determined before any memory addresses are guessed.
All word-cost exponents here are independent of the program, $b$ and
$\rho$.

\subsection{Wiring the smaller clock}
\begin{proof}[Proof of Theorem~\ref{thm:gap}(1)]
The wiring has two tasks. It checks a public trace and sends each
private job in that trace to one processor. To supply the jobs in
fixed positions, augment the public computation with a stream in a
separate address space inaccessible to the simulated program. Write the stream in consecutive chunks of
$J=\lfloor n/8\rfloor$ bits. Each job record has a start bit, a
fixed-width length and the data bits. Before starting a job, open a
new chunk if the largest allowed record would not fit. Zero-pad the
unused part of each chunk. For large $n$ a maximum record uses less
than half a chunk, so each closed chunk, apart from the last one, is
at least half full. No job crosses a chunk boundary.

Reserve space for the worst case, allowing $O(w)$ record expansion
per charged job bit or job, a factor of two for padding, and rounding
to a whole chunk. The begin, append-bit and end instructions implement
the stream and reject malformed jobs. Checking the public cost cap
and record lengths, writing the stream, and reading every reserved bit
in address order on halt take $T_-\operatorname{poly}(w)$ ticks.
This control is fixed for all public descriptions. Use the trivalent
trace network at this enlarged cap, with the final stream reads as
outputs.

Wire each chunk's $J$ outputs to one $U_n^*$ whose remaining slots
specify a fixed program, $n$ if needed, and padding. The program
executes its jobs in order, accepting exactly when all accept. A
successful path selects one accepting private path for each job.
A chunk contains at most $n$ jobs of at most $n^{20}$ private steps
each. Even scanning the used tape at every simulated private step
keeps the joint computation within the processor's $n^{100}$ steps.

Pack the trivalent vertices into further processors. A fixed checking
program reads, for each assigned vertex, its truth table and at most
three incident bits. A padded list with its length specified costs a
constant number of bits per vertex apart from its header. Add trivalent
vertices as needed to connect literal pins and chunk bits. Wire the
two slots of each link, including when both belong to checks on a
single processor. Its program thus needs only the
individual checks.

The chunks and the packed vertices together use
\[
           O\bigl(1+T_-\operatorname{poly}(w)/n\bigr)
\]
processors. For fixed $b$, take $\rho$ small enough that this is less
than $R$ and pin the remaining processors to return one. Generating
the trace, packing it in order and listing its pins and pairs costs
$Rn\operatorname{poly}(w)$ ticks. On contraction, a valid trace has
one extension to the public checking bits. The chunk processors count
all independent choices of accepting private paths. Their product
gives precisely the number of completions of that public path.
\end{proof}

\subsection{Compressing the larger clock}
\begin{proof}[Proof of Theorem~\ref{thm:gap}(2)]
Set $U=T_+$. We put its verification trace on the cyclic butterfly
of Theorem~\ref{thm:closed}, using the fixed XOR links to generate the
wiring within the smaller public budget \cite[Section~\ext{G}{app:trace}]{oaG}.

\emph{Routing the rows.}
Fix $a=1/20$, and group consecutive rows into buckets of
$S=2^{\lfloor a\log n\rfloor}$ rows. The low row bits give the offset
in the bucket and the high bits give its lane. Put one bucket at
cursor zero on each lane. The two arcs to the next cursor position
either stay on the lane or flip its active coordinate. Add empty
parallel slots to reach the compressor's fixed degree. To route a
link with a prescribed inter-bucket XOR mask, traverse the cursor
cycle, flipping each coordinate
specified by the mask when its turn comes. At every cursor step
this induces a bijection of lanes. The coordinates already flipped
recover the source lane from a relay's lane.

Reserve slots on each arc by mask, exchange type, low source row
bits and bit number. A relay equates successive bits on an occupied
route and pins surplus slots to zero. Thus every satisfying trace
has a unique extension to the routed network. Its number of positions
$P$ and maximum arc width $D$ satisfy
\[
              P=O((U/S)w^{O(1)}),\qquad D=O(Sw^{O(1)}).
\]
The slot ranges give all capacities, arc identifiers, lane displacements
and bucket ranges in $P\operatorname{poly}(w)$ public work. They
require no scan of the $S$ low row indices of each bucket.
Description words enter only as literals of initialization rows. On
given arguments a bucket or relay factor counts locally guessed
assignments that pass the row predicates, wire equalities and padding
checks. Generating row descriptions and making these checks costs
$O(S\operatorname{poly}(w))$ per local random access path. At a relay
the source lane and type schedule identify occupied slots independently
of the description words. The same mask, type and row identifiers
identify the slots at their destinations.

\emph{Choosing the shared bits.}
Apply $t=\lceil z\log n\rceil$ blocks from Section~\ref{sec:compr}.
Choose $z>0$ rational and small enough that $C^t\le S^{1/4}$ for all
large $n$. The available $\log(U/S)+O(1)$ lane coordinates suffice.
The public machine generates the initial descriptors and reads the
guesses for a history in their canonical encoding. It checks admission
and constructs the recursive descriptors. Position and capacity
indices have $O(w)$ bits, so this costs $P\operatorname{poly}(w)$ ticks
per path. Guess the surviving bits once per arc, a total of at most
\[
              2^{-t}O(P(D+1))\le U n^{-z}\operatorname{poly}(w).
\]

\emph{Issuing the jobs.}
For each final position, issue its final arguments and the recursive
layout of its local contraction. Follow the recorded preimages to
generate slices, ties and internal field identifiers, using local
indices for each job. For example, write a new local index into an
array at the original identifier upon visiting it. Each occurrence
takes $\operatorname{poly}(w)$ public work. The layout gives each
initial bucket or relay by its identifier and size parameters. Supply
also the literal description words in the intersection of its bucket
range with the initialization segment. The public machine visits
actual literals and skips padding. The hidden field bits will be
guessed in private.

For each compressor scalar, guess a fixed-length exponent for each
position and supply it to that job as a claim. Sum these exponents in
public counters. Each path uses only polynomially many map entries
and polynomial monomial degrees in $C^t,D,w$ per position, and there
are at most $P$ positions. Canonical encodings of $O(w)$ bits
therefore suffice for every local exponent and its global total.
Each job has length at most
\[
                 (C^t)^2(D+S+1)\operatorname{poly}(w)<n^{1/4}
\]
once $\rho$ is small enough. Allowing for all steps, including issuing
each job, the public cost is
\begin{equation}\label{eq:dump}
                 (C^t P+R+U n^{-z})\operatorname{poly}(w).
\end{equation}
The first term pays for $O(C^t)$ layout records per final position
and the exponent claims. The second pays for the description words,
since each initial bucket moves intact to exactly one final position.
The third pays for the final arc bits supplied at their two endpoints.
The polynomial overhead covers addressing, guesses and admission checks,
including rejected paths.

\emph{Checking the private counts.}
Each private machine guesses the hidden fields, checks local map
patterns, and runs the path tests for its moved initial factors. It
accepts only if the map entries are allowed and their monomial
exponents sum to the supplied claims. There are at most $C^t$
occurrences and $(C^t)^2 O(D+1)$ field bits in one job. Slices and
ties can be checked bit by bit. For a map, let $m$ be the largest field length in its
motif. The selector scans in
\cite[proof of Theorem~\ext{G}{thm:conversion}]{oaG} test the entry
and determine its monomial exponents in $(1+m)^2\operatorname{poly}(w)$
steps.

For a fixed $q$, the bound $B=(C^t)^2(D+S+1)w^q$ covers both the
combined length of a job's expanded records and fields and the
largest field length in any of its motifs. A private RAM can store
the local records in arrays indexed within the job. Charge each occurrence
$B^2$ steps for scans and $O(S)$ for any original row group it holds,
with fixed polynomial factors in $w$. Since $S\le n^a$ and
$C^t\le S^{1/4}$, the resulting bound is
$C^t(B^2+S)w^{O(1)}\le n^{13a/4}w^{O(1)}<n^2$ for small $\rho$. Scanning lists
on tapes incurs at most quadratic simulation overhead, leaving room
for universal interpretation within $n^{20}$ private steps. A fixed
private program implements these tests. The public machine discards
nonadmitted histories before issuing jobs, and the checks can reject
malformed data within the caps.

Fix a history and all final arc bits. Write $\boldsymbol\alpha$ for
the list of compressor scalars and $N_j(e)$ for the number of accepting
private paths at position $j$ when its exponent claim is $e$. The
local contraction at $j$ is $\sum_e N_j(e)\boldsymbol\alpha^e$.
For publicly claimed exponents $(e_j)_j$ there are
$\prod_j N_j(e_j)$ completed successes, each with public weight
$\boldsymbol\alpha^{\sum_j e_j}$. Summing these contributions gives
the product of the local contractions. Summing over arc bits and
admitted histories now gives the original acceptance count by the
compressor identity.

It remains to fit the public cost into the smaller clock.
Since $C^t\le S^{1/4}$ and $U=T_+$, Equation~\eqref{eq:dump} is
\[
             O\bigl(U(n^{-3a/4}+n^{-z}+n^{-1-b})w^{O(1)}\bigr).
\]
Choose $b>0$ sufficiently small after $a,z$, for example
$2b<\min\{a/2,z\}/2$. The required saving is
$T_-/U=\Theta(n^{-2b})$. Since $w=\Theta(n^\rho)$, decreasing $\rho$
makes the displayed cost $o(T_-)$ and meets the wiring, job-length
and private-time bounds above. The finitely many exponents on $w$
in these estimates are fixed before $b,\rho$ are chosen.
\end{proof}
~~~

## P / batch.tex

OpenAI, preprints/Circuit-lower-bounds-for-Hamiltonian-cycle-and-the-permanent-November-24-2026/build/batch.tex, Git blob dbf18a24cd0bac61348c1a3508839bfa5b0f9920. Complete proof with its clock diagram omitted.

~~~latex
\section{Hardness of batches}\label{sec:batches}

The two parts of Theorem~\ref{thm:gap} give two ways to recover the
count specified by one description. Wiring the signed processors realizes the
success count of a public program. Substituting a hypothetical batch
network into that wiring lets an evaluator recover the count with
enough time left, after compression, for the program to add one.

\begin{proof}[Proof of Theorem~\ref{thm:batches}]
Fix the finite field $F$ and write $p=\operatorname{char}F$.
Choose the compressor scalars in $\mathbb F_p$. Take $b,\rho$ as in
Theorem~\ref{thm:gap}, decreasing $\rho$ below as needed, and take
$0<\delta<b/2$ and $0<\sigma<\rho$. Suppose that for some sufficiently
large $n$ and some $1\le k\le 2^{n^\sigma}$ a trivalent network
$\mathcal B$ of volume at most $kn^{1+\delta}$ realizes
$(U_n^*)^{\otimes k}$ over $F$. Put
\[
                    R=k\,2^{\lceil n^\rho\rceil}.
\]
This choice meets the clock hypotheses and partitions the processors
into $R/k$ whole batches. Use the wiring of Theorem~\ref{thm:gap}(1)
with the private machine fixed by Theorem~\ref{thm:gap}(2).
For a description $D$, let $Z(D)\in\mathbb F_p$ be the value of this
wiring with its pins set to $D$. It is the integer success count of the
stored program, reduced modulo $p$.

Reserve fields in $D$ for its public code, the size parameters, and
an indexed description of $\mathcal B$ with all its local entries in
$F$. Each tensor has bounded degree and $F$ is fixed, so the network
uses $O(kn^{1+\delta})$ words, even giving each wire endpoint and each
local entry its own word. All active fields, including the code chosen
below, fit in $\lfloor R/10\rfloor$ words for large $n$. The public
machine can read its entire description, including its program.

\smallskip
\noindent\emph{Evaluate the substituted network.}
Fix an evaluator $E$ whose input is a description $D$ with this network
and these size parameters. Generate the pins and links for $D$, and
replace each consecutive group of $k$ processors by a copy of
$\mathcal B$. This preserves the scalar in $F$ for every choice of the
public code in $D$. Guess each summed wire bit once and multiply the
selected tensor entries in $F$, enforcing the pins and links.

To turn a product $u\in F$ into an acceptance count, fix an
$\mathbb F_p$-linear map $\lambda:F\to\mathbb F_p$ with $\lambda(1)=1$.
Guess $\lceil\log p\rceil$ bits and accept when the resulting integer
is smaller than the representative of $\lambda(u)$ in $[0,p-1]$.
Summing over assignments applies $\lambda$ to the entire contraction.
Since $Z(D)$ lies in $\mathbb F_p$, the evaluator satisfies
\[
                    \#E(D)\equiv Z(D)\pmod p ,
\]
where $\#E(D)$ is its integer number of accepting paths and the
congruence identifies $\mathbb F_p$ with residues modulo $p$.

The evaluator only reads the stored instructions as literal data
when setting the pins. Generating the genuine wiring takes
$O(Rn\operatorname{poly}(w))$ ticks by Theorem~\ref{thm:gap}(1).
Each exterior slot has an $O(w)$-bit identifier for its batch and
corresponding leg. Traversing the $R/k$ copies of $\mathcal B$,
checking the assignments and multiplying their entries therefore takes
$O(Rn^{1+\delta}\operatorname{poly}(w))$ ticks in all. The exponent
on $w$ is fixed, so $\delta<b/2$ and a sufficiently small $\rho$ put
this procedure within $T_+$, without private jobs.

\smallskip
\noindent\emph{Choose the program whose pins are evaluated.}
Give the public program two branches. The first accepts once.
The second runs the compressed simulation of $E$ on its actual $D$
using Theorem~\ref{thm:gap}(2). All layout generators, field offsets
and description sizes can be chosen before installing this
constant-length program. To present $D$ to the evaluator, copy the
active words, including the stored instructions, and fill unused
words with zeros. This costs $O(R)$ ticks. The evaluator generates
pins from these words without any recursive call to the program.

On a successful public path of the simulation, evaluate the
compressor monomial in $\mathbb F_p$ from its public exponent tallies.
Use the same binary-guess test as above to accept with this
multiplicity modulo $p$. Every successful completion of the private
jobs receives the prescribed public weight. The simulation takes
$o(T_-)$ ticks, leaving room for the copy, the prime-field arithmetic
and the extra branch. Thus the program at cap $T_-$ has count
$1+\#E(D)$ modulo $p$. Its genuine wiring has value $Z(D)$, while
the substituted wiring used by $E$ has the very same value.
Figure~\ref{fig:clock} displays the resulting contradiction
\[
                Z(D)=1+\#E(D)=1+Z(D)\quad\text{in }\mathbb F_p .\qedhere
\]
\end{proof}

% Editorial omission of the clock diagram.
~~~

## G / exact.tex

OpenAI, preprints/Euclidean-geometry-in-the-asymptotic-spectrum-of-tensors-November-24-2026/build/exact.tex, Git blob b15103db171da47d6d3292d3d0126d1c8a98b9c8. Complete exact descent and its length lifting, sorting and evaluation steps.

~~~latex
\section{Exact descent with tagged-permutation gates}\label{sec:exact}

In Theorem~\ref{thm:compressor}, each summand uses few wire bits. The
available slots may nevertheless occupy a large network, and the local
contractions still have to be evaluated. We realize the full sum in
one network by compressing a batch of copies together. Copies share
the cost of specifying the wiring, and short descriptions of the local
contractions permit simultaneous table lookup. The only gates of
unbounded arity are tagged-permutation tests of a single length.

Call a trivalent network \emph{normalized} if it has no self-loops and
each open bit is carried by a separate equality leaf adjacent to an
interior vertex. Inserting equality leaves and subdividing self-loops
changes the size by at most a constant factor.

\begin{proposition}[Descent with tagged-permutation gates]
\label{prop:oracle-descent}
Let the finite field have ports in the sense of Definition~\ref{def:ports}.
There is an $a\in(0,1)$ with the following property. Fix $J\ge 2$.
For every sufficiently large $h$, put $N=2^h$, let $N\le m\le N^J$,
and let $G$ be
a normalized trivalent network of size at most $m$ with $e$ open bits.
There is a positive integer $d$ with $\log_2 d=O((\log h)^2)$ and an
exact network for $G^{\otimes d}$ of size
\[
                      O\bigl(d(m/h^a+e\operatorname{polylog} h)\bigr)
\]
whose gates are trivalent except for calls to $\operatorname{Tag}(h,h)$.
The size charges every incident bit wire, including those at calls.
The bounds are uniform in $m$ and $G$.
\end{proposition}

\subsection{Length lifting}
The following version of the column-sorting construction of
Leighton~\cite{Leighton} records the total number of short sorts as
well as their depth.
\begin{lemma}[Lifting the length of a sort]\label{lem:sort-lift}
Let $b\ge 16$ be a power of two. A sort of length $l\ge b$ can
be carried out by layers of independent sorts on $b$ records,
interspersed with fixed permutations, insertions of sentinels, and
deletions of padding. For an absolute $C$, with
$u=1+\log(l)/\log b$, there are at most $C u^C$ layers, using altogether
at most $C(l/b)u^C$ sorts of length $b$. For stability it suffices to
carry the original position through the layers.
\end{lemma}

\begin{proof}
Pad to a power of two $L\in[l,2l)$. First suppose $L=kb'$,
with $b'$ and $k$ powers of two and $k^2\le b'/2$. Arrange $k$
columns of height $b'$. Sort each column. Deal its sorted entries
cyclically into the $k$ new columns and concatenate the contributions
within each new column. Each destination receives $b'/k$ entries from
each source. Sort the new columns and list them in row-major order.
For any cut in the total order, the counts of elements below the cut
in any two new columns differ by at most $k$, since the counts from
each source differ by at most one. In row-major order the binary
string indicating membership below the cut thus has a prefix of ones
and a suffix of zeros separated by a window of at most $k^2$ bits.

Sort aligned blocks of $b'$ consecutive positions, then blocks shifted
by $b'/2$, handling the end half-blocks by padding. A window inside
an aligned block is resolved in the first pass. A window crossing an
aligned boundary is contained in the shifted block centered on that
boundary. Sorting the aligned blocks leaves the already ordered prefix
and suffix outside the window intact, so the second pass resolves the
cut. This holds for every cut. The four passes consequently sort the
list. Each pass uses at most $2L/b'$ calls, a bound allowing for the
two end blocks.

Starting from $b_0=b$, set $b_{i+1}=2^{\lfloor(5/4)\log b_i\rfloor}$.
For $b_i\ge16$ these powers increase, and their logarithms grow at
least geometrically (with ratio at least $8/7$, for example).
Whenever $b_i<L\le b_{i+1}$ we may take $b'=b_i$ in the preceding
construction. At each recursion take the largest $b_i$ strictly below the current
length as the child length, and pad only the end half-blocks.
There are $O(1+\log u)$ levels. The number of sequential passes
increases by a factor at most four per level and their combined padded
record volume by a factor at most eight. At the bottom each full call
uses $b$ record slots. These bounds give the stated constants. Appending
the original position once orders all genuine records distinctly. A
flag distinguishes sentinels and permits their deletion, including
in intermediate calls.
\end{proof}

Apply Lemma~\ref{lem:sort-lift} with base length $b=N=2^h$.
A stable sort of $N$ records of $O(h)$ bits costs $O(Nh)$ using
$\operatorname{Tag}(h,h)$ gates and charging their incident wires
(Section~\ref{sec:width}). If a payload exceeds $h$ bits, split it into a bounded number of
fields and use the same tags for all fields. For $l<N$ pad to
$N$. For $N\le l\le N^{O(1)}$ the ratio $\log l/\log N$ is bounded,
so a stable sort costs $O((l+N)h)$ throughout this range. This bound
also gives simultaneous lookup when $l$ counts the query and
dictionary records together. Sort by key, placing the unique
dictionary record first in each group. Scan the group to answer its
queries, then return the answers in their original order by a second
sort. Each valid input and output has a unique witness.

\subsection{Realizing the compressed expansion}

\begin{proof}[Proof of Proposition~\ref{prop:oracle-descent}]
A compressed bit has two uses. At its vertex it is part of the
argument to a local factor. On its arc it must agree with the bit
at the other endpoint. Vertex lookup requires short keys, since the
table size grows exponentially in the key length. Arc comparison
uses addresses of $\Theta(h)$ bits, so each comparison should cover
a long string. We stack the copies of $G$ to obtain both orderings
within the wire budget.

\paragraph{Stacking the copies.}
Let $d_0\ge4$ be the degree in Theorem~\ref{thm:compressor}. Give each
internal edge of $G$ an arc in both directions, with the wire on just
one arc. Copy this loop-free topology with all widths zero, no
exterior indices, and scalar factor one at each duplicate vertex. Each
vertex and its duplicate have the same in- and out-degree, so we can give
both degree $d_0$ by adding zero-width arcs in both directions between
them. Keep the exterior bits as free indices of their original vertices.

Take $s=\lceil 2a\log h\rceil$ compressor blocks, consuming $b_0 s$
lane coordinates. We choose the constant $a>0$ below. Set
\[
 L=2^{b_0s},\qquad B=\lceil h^{1/2}\rceil,\qquad
 R=2^{\lceil4\log h\rceil},\qquad d=B L R .
\]
Since $b_0$ is fixed, $d=h^{O(1)}$. Use $LR$ lanes with zero initial
displacements and reserve the $\log_2 R$ replication coordinates
throughout. Each lane starts with a tensor product of $B$ copies of
$G$, with $B$ bits per nonempty arc. Put $P=mLR$. There are $O(P)$
initial positions, the maximum width is $B$, and the square root of
the energy is $O(PB)$. Hence Theorem~\ref{thm:compressor} bounds the
total internal wire width in each term of the final expansion by
$O(dm/h^{2a})$.

One binary choice per shrinking substep specifies the \emph{layout
pattern}. Its $O(s)$ bits determine every unweighted topology and
placement in the sequence. Each substep translates by previous arc
displacements or moves locally on the remaining lanes, so the
replication coordinates still have zero displacement. The final graph
therefore consists of $R$ copies of the same set of positions, which
we call \emph{hosts}. Widths and tensors can vary across copies. For
each exterior bit position of $G$ and each lane in $[L]$, the layout
pattern fixes a single host carrying one exterior port of width $B$
per replica. At most $2^{O(s)}$ such ports reach one host, since each
substep merges a bounded number of positions.

\paragraph{Encoding the local contractions.}
For each host and replica we record two kinds of data.
\begin{itemize}
\item The \emph{operands} are the bits in its incident slots and exterior
ports. Their total length is at most $2^{O(s)}B$ by
Theorem~\ref{thm:compressor}.
\item The \emph{program key} describes the factor to evaluate on the
operands and gives the length of each field. Unfold the contraction
records of Lemma~\ref{lem:descriptors} using local names for the
constituent objects. The key lists field offsets, incidence ordinals,
widths, ties, motif templates and term indices. Each leaf carries one
label for the original trivalent factor over the fixed finite field,
which determines its $B$-fold product. Lengths and term indices use
$O(\log h)$ bits apiece, giving $2^{O(s)}\operatorname{polylog} h$
bits for the key.
\end{itemize}

Boolean circuits construct the keys from guesses for the motif terms.
To count each history exactly once, select one layout pattern, check
all its guesses, and set the guesses for the other patterns to zero.
Use a circuit template for each pattern with the graph sequence, greedy
colors, offsets and slot lists hardwired. The only global calculations
are the sums of energies that verify the layout choices. After removing
the common mass factor, all terms have a fixed denominator and
$O(\log h)$-bit numerators. Balanced addition costs
$O(P\operatorname{polylog} h+h)$ per sum, as the word length grows by
one bit per level. Apply the prescribed tie rule.

The width calculations, term checks and slice intersections at one
position use a bounded number of operations on quantities of
polynomial size in $h$. By Theorem~\ref{thm:conversion}, their cost is
$\operatorname{polylog} h$ bits and gates. Fixed wires supply the
local data from preceding steps. Gathering the $2^{O(s)}$ records for
a final kernel, checking for rejection, and copying the selected
outputs over all patterns costs $P\,2^{O(s)}\operatorname{polylog} h$.
A long name, such as an exterior port identifier, is read once per
host slot and shared over all replicas. For fixed $J$ the names use
$O(h)$ bits. Selecting them from the hardwired lists costs an additional
$(P/R)\,2^{O(s)}O(h)$.

Fix $a<0.02$ small enough that all preceding
$2^{O(s)}\operatorname{polylog} h$ bounds are $O(h^{.1})$.
A fixed-format array of program keys takes $O(Rh^{.1})$ bits per host,
including padding for the largest topology. Use base order within each
pattern and pad to $O(P/R)$ hosts, with blank keys encoding unit
factors. Each key gives the operand lengths in local slot order,
including the exterior ports.

\paragraph{Matching the wires.}
Use rows of $h$ bits. For each incident slot at a host, concatenate its
operand bits in replica order to form a stream. An internal stream
must equal the stream at the other end of its arc. An exterior stream
of $BR$ bits must equal the corresponding input bits, which are ordered
by their fixed names for every layout pattern. For each host reserve
control rows for the program keys, stream identifiers, flags and
padding. These occupy $O(R h^{.1}+h^{1.1})$ bits. Reserve whole rows
for each stream, padding its final row after combining its replicas.

The total allocated bit count is
\[
 O\big(P h^{.1}+(P/R) h^{1.1}+d(m/h^{2a}+e)\big)
   =O(S_*),\qquad S_*=d(m/h^{2a}+e).
\]
Each allocation needs at most $O(Rh^{.6}+h^{1.1})$ bits even counting operand bits.

Pack the host allocations in order into bins of $h^8$ rows, starting
a new bin when the next host does not fit. Each host uses at most half a
bin, so a fixed pool of $O(S_*/h+h^8)$ rows suffices. Computing the
allocation lengths costs $O(P h^{.1}\operatorname{polylog} h)$ by
summing the stream sizes given in the keys. Each such local size is
polynomial in $h$.

Scanning the hosts to assign their bins costs $O(h)$ per host for the
bin index and occupancy. Sort the control rows into their prescribed
positions, mark each host's start, and assign a blank control flag to
the remaining rows. An absolute row address uses $O(h)$ bits and can
be updated in linear work by concatenating the bin and offset fields.
The scatter costs $O(S_*+Nh)=O(S_*)$, since $m\ge N$, $R\ge h^4$
and $S_*\ge BLR\,m/h^{2a}$. Sum over guesses for
the stream rows, requiring zeros in all padding and outside the
streams. Parse each allocation deterministically from its controls.
A stream row receives its matching name, side and ordinal, copying
the common prefix from its stream header. One sort and scan verifies
the equality of every pair of internal streams and matches the
exterior streams with fixed input rows. The latter can be rounded to
whole rows per exterior stream within the same $O(S_*/h)$ row budget.

\paragraph{Evaluating the vertex factors.}
For each host and replica, join its program key and operand bits in
operand order into one self-delimiting record. Reserve a length field
in the key and fix the encoding of empty padding. Each record has
$O(h^{.6})$ bits and fits into a row of $\lfloor h/2\rfloor$ bits for
large $h$. Within each bin the sum of record lengths is at most a
constant times the bin capacity. Pack records in order, beginning a
new row when the next record does not fit and padding unused bits
uniquely. A fixed array of $O(h^8)$ rows per bin suffices.

A bit RAM with $O(\log h)$-bit addresses performs these local operations
in at most the bin capacity times $\operatorname{polylog} h$ steps.
It scans the controls, checks padding, and maintains a current offset
for each stream to copy the fields into replica order. Copying the
long headers bitwise costs $O(h)$ transfers per row. As in
Lemma~\ref{lem:tag-amplification}, bitonic sorting of time-stamped
memory records together with local instruction checks verifies the
deterministic trace by a trivalent network. The overhead is
$\operatorname{polylog} h$ since the bin size is polynomial in $h$.

Each packed row specifies the product of the vertex factors in its
records, with value one for an empty row. The fixed finite field and
the local contractions determine a table of at most $\sqrt N$ entries
with constant-sized outputs. Hardwire this table and perform
simultaneous lookup with one scalar factor per request. The cost
is $O(S_*)$, including the table, since all lists have length
polynomial in $N$ for fixed $J$.

Fix the selector bits and exterior indices. The surviving stream
assignments are in bijection with assignments of the final wire bits.
Every deterministic check has a unique valid completion or is used
as its exact relation tensor. Summing the products of the returned
scalars over these assignments and over the selectors gives
$G^{\otimes d}$. The total size is $O(S_*\operatorname{polylog} h)$.
For fixed $a>0$, the factor $h^{-2a}\operatorname{polylog} h$ in the
interior cost is $O(h^{-a})$, as claimed.
\end{proof}
~~~

## G / width.tex

OpenAI, preprints/Euclidean-geometry-in-the-asymptotic-spectrum-of-tensors-November-24-2026/build/width.tex, Git blob 9284d3abd9a8e50a835e7fa9d5e5a363fd119dd1. Stable-sort implementation. Approximate amplification is unused.

~~~latex
A stable sort on $2^d$ records of $v$ bits
also has a network consisting of one $\operatorname{Tag}(d,v)$ call and
$O(2^d(d+v))$ additional size. Guess an input position for each proposed
output, check the permutation of the payloads with the call, and check
nondecreasing keys with increasing input position at ties. Each valid
input--output pair has exactly one witness.
~~~
