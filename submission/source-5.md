# Source evidence 5. Distinctness, dependency map and license

Excerpts from OpenAI's mathematics preprints at commit 92b516ef8e2abf66c3db7aa31438d3010bc7956f. The dependency map, bibliography and Apache 2.0 license are in source-5.md. The authors' result labels are preserved.

https://github.com/openai/math/tree/92b516ef8e2abf66c3db7aa31438d3010bc7956f

- M: The-matrix-multiplication-exponent-is-nine-fourths-November-24-2026
- G: Euclidean-geometry-in-the-asymptotic-spectrum-of-tensors-November-24-2026
- P: Circuit-lower-bounds-for-Hamiltonian-cycle-and-the-permanent-November-24-2026

## How the proofs fit together

M Theorem 2.2 (thm:many) is applied over F_2 at t=1. It asserts a lower bound for all sufficiently large widths without a matrix-multiplication hypothesis. M list.tex proves it by considering separately fields with and without ports.

Without ports, source-1 gives spectral prices, wire monotonicity, rescaling and ultralimit closure. The graph entropy and product theorems in source-3 then imply M's uniform-destinations lemma. Applying it to the tagged permutation obtained from two distinctness tests yields the bound.

With ports, the spectral comparison, extension, shifted crossing and CAT(0) results in source-1 support the compressor and descriptor proofs in source-2. Source-4 supplies the clock gap and exact-batches contradiction for the signed universal counting tensor, and exact tagged descent. The last subsection of M list.tex combines these results, iterating descent to contradict the exact-batches theorem.

M's reference to sec:prices-background invokes the prices definitions, transport and closure supplied by source-1. The change between binary and natural logarithms is a constant rescaling. G exact.tex invokes the memory consistency procedure in lem:tag-amplification. The sorting and scan procedure is reproduced in P gap.tex. Approximate amplification is unused here.

Conventional background includes linear algebra, finite-dimensional minimax, compactness, ultrafilters, Zorn's lemma, finite Ramsey theory, standard CAT(0) barycenter, projection and flat-quadrilateral properties, and bounded-degree scalar expanders. The sources retain the original labels so that each imported result can be followed to its argument. The LaTeX command ext with a paper key and label denotes the corresponding result in that paper. EE denotes expectation. The separate proof.md supplies the transfer to full B2 and padding to every input length.

## M / obstruction.tex

OpenAI, preprints/The-matrix-multiplication-exponent-is-nine-fourths-November-24-2026/build/obstruction.tex, Git blob e817b4580eb69a590234ddd2dc50982d366e6939. Trivalent network model, D_w and Theorem 2.2. The original opening states: "Fix a finite field k." This standing hypothesis applies to the entire following excerpt; this revision restores it explicitly.

~~~latex
\subsection{Two estimates for networks on bits}
A \emph{bit tensor} of arity \(u\) is a function \(\{0,1\}^u\to k\).
In a tensor network each vertex carries such a table. A wire pairs two
argument positions, and each unpaired position is a named exterior bit.
For a fixed setting of the exterior bits, multiply the vertex entries
and sum over all settings of the internal wires. This gives the tensor
represented by the network. A network is \emph{trivalent} if each
vertex has arity at most three. Its \emph{size} counts the vertices and
exterior bits. Its \emph{volume} counts the vertices and all wires,
including the exterior wires. These counts agree up to an absolute
factor by trivalency. For networks of unrestricted arity we count
vertices plus their total arity.


% Editorial omission of the separate conditional compression theorem.
For \(w\ge1\), let \(D_w\) be the zero-one tensor for a list of \(2^w\)
words of \(2w\) bits each. Its value is one exactly when the words
are pairwise distinct. It has \(n_w=2w\,2^w\) exterior bits. The tensor
\(D_w^{\otimes t}\) tests distinctness separately within each of \(t\)
independently specified lists. The network may connect vertices across
lists. Let \(\operatorname{vol}_k\) denote minimum trivalent volume over \(k\).

\begin{theorem}[Volume of list distinctness]\label{thm:many}
There is a constant \(\sigma>0\), depending only on \(k\), such that for
all sufficiently large \(w\) and every integer \(1\le t\le 2^w\),
\begin{equation}\label{eq:many}
                  \operatorname{vol}_k(D_w^{\otimes t})
                          >t n_w w^\sigma .
\end{equation}
\end{theorem}
~~~

## M / list.tex

OpenAI, preprints/The-matrix-multiplication-exponent-is-nine-fourths-November-24-2026/build/list.tex, Git blob 2d09e2ab48a3f3435ffb68122d51a57b0d63f6fb. Complete proof of Theorem 2.2 in both cases.

~~~latex
\section{The cost of testing many lists}\label{sec:list}

Two distinctness tests suffice to check that records have reached their
specified destinations. This gives a common starting point for the two
cases of Theorem~\ref{thm:many}. In a field without ports, the possible
destinations force many wires to have large total length in a suitable
spectral space. In a finite field with ports, cheap distinctness tests
would make it possible to compress a family of tensors with a known
lower bound.

\subsection{Checking a tagged permutation}

Write \(h=w-1\), \(N=2^h\), and index records by the \(h\)-bit words.
The tagged-permutation tensor \(\operatorname{Tag}(h,h)\) of
\cite[Definition~\ext{G}{def:tag}]{oaG} takes
\(N\) records \((j_i,x_i)\) and an ordered list of \(N\) payloads
\(y_j\), with \(h\) bits per tag or payload. Its value is one when
\(i\mapsto j_i\) is a permutation and \(x_i=y_{j_i}\) for every \(i\).
Its other entries are zero.

The first \(D_w\) test checks the tags. Encode each \(j_i\) as
\((0,j_i,0^{h+1})\) and fill the remaining \(N\) positions with the
distinct words \((1,j,0^{h+1})\). All \(2N\) words have length \(2w\),
and the test passes precisely when the tags form a permutation \(\pi\).

The second \(D_w\) test checks the payloads. Its first \(N\) words are
\((j_i,0,x_i,0)\). At the position for \(j\) in the second half, sum
over two choices, the reserve word \((j,1,0^h,0)\) with weight \(1\)
and the comparison word \((j,0,y_j,0)\) with weight \(-1\).
For tags forming \(\pi\), the only possible collision at address \(j\)
is between the comparison word and the word from \(i=\pi^{-1}(j)\).
Thus the signed sum factors as
\[
     \prod_j \bigl(1-\mathbf 1_{\{x_{\pi^{-1}(j)}\ne y_j\}}\bigr)
       =\prod_j\mathbf 1_{\{x_{\pi^{-1}(j)}=y_j\}} .
\]
This identity holds in every characteristic. A summed selector bit for
each address controls its word and its weight. Copying the tags,
selecting the words, and pinning the constants cost \(O(Nh)\) trivalent
volume. Consequently, if \(D_w^{\otimes t}\) has a realization of
volume \(V\), two uses of that realization give a network for
\(\operatorname{Tag}(h,h)^{\otimes t}\) of volume at most
\(2V+O(tNh)\), including when the realization shares work among lists.

\subsection{Fields without ports}

A class of Hadamard spaces has \emph{ports} if it admits regular
graphs of unbounded order and common bounded degree with the following
property. For every placement of their vertices into any space of the
class, the mean squared distance of two independent uniform vertices
is at most a fixed constant times the mean squared edge length.
A field has ports when its Hadamard spectral spaces have this property
\cite[Definition~\ext{G}{def:ports}]{oaG}.

\begin{lemma}[Uniform destinations]\label{lem:uniform-demand}
Let the field have no ports. Consider a fixed trivalent topology with
\(r\) input bits and \(r\) output bits. Suppose that for every permutation
in the support of a distribution, some choice of vertex tables realizes
the \(r\) independent wires joining each input to its assigned output.
If each input has a uniform destination on a set of \(q\ge2\) outputs,
the topology has size \(\Omega(r\sqrt{\log q})\). The destination sets
may depend on the input.
\end{lemma}

\begin{proof}
Give each exterior bit its own equality leaf. Connect the resulting
components in a path by adding edges to the graph used to measure
distances. This increases the maximum degree by at most two and gives
\(m\) edges, with \(m\) at most a constant times the network size.
Let \(\alpha\) be the law of a uniform oriented edge, and let \(\mu\)
be the joint law of a uniform input and its random destination.
For any placement in a spectral space with distance \(d\),
\[
                    r\,\mathbb E_\mu d\ \le\ m\,\mathbb E_\alpha d .
\]
Indeed, for each permutation the network tensors locally contract
the product of its internal bit wires to the prescribed exterior
wires. Monotonicity and additivity of prices bound the sum of the
destination distances by the sum of the internal edge lengths.
Average over the permutations and include the extra edges on the right.

We need a placement where many destinations are far away at once.
The class of Hadamard spectral spaces is closed under isometry,
rescaling, closed convex subspaces, and pointed ultralimits
(Section~\ref{sec:prices-background}). The product theorem
\cite[Theorem~\ext{P}{prod:main}]{oaP} gives a subclass with the same
properties and no ports, closed also under finite \(\ell_2\) products.
Let \(\mathcal K\) be the squared distance tables \(D=d^2\) realized
in this subclass with \(\mathbb E_\alpha D\le1\). Products of rescaled
placements give convex combinations of tables. Connectivity bounds
every entry, and pointed ultralimits realize coordinatewise limits.
Hence \(\mathcal K\) is convex and compact.

For a probability law \(\nu\le 2\mu\) on ordered pairs, let
\(\beta_\nu=\max_{D\in\mathcal K}\mathbb E_\nu D\). By rescaling,
\(\mathbb E_\nu d^2\le\beta_\nu\mathbb E_\alpha d^2\) for every
placement in the subclass. Write \(X,Y\) for the endpoints under
\(\nu\). The graph entropy bound \cite[Corollary~\ext{P}{cor:graph-entropy}]{oaP}
and the bounded degree of our graph give
\(H_\nu(Y\mid X)\le C(1+\beta_\nu)\), with entropy in bits.
Every pair has \(\nu\)-mass at most \(2/(rq)\), so
\[
    H_\nu(Y\mid X)=H_\nu(X,Y)-H_\nu(X)
               \ge \log_2(rq/2)-\log_2 r=\log_2(q/2).
\]
It follows that \(\beta_\nu\ge c\log q\) for all sufficiently large
\(q\), uniformly in \(\nu\).

Bilinear minimax on \(\mathcal K\) and the compact convex set of
laws \(\nu\le2\mu\) gives a single \(D=d^2\in\mathcal K\) for which
\(\mathbb E_\nu D\ge c\log q\) for every such law. If the pairs with
\(D<c\log q\) had \(\mu\)-mass at least \(1/2\), conditioning \(\mu\)
on those pairs would violate this bound. Hence
\(\mathbb E_\mu d\ge \tfrac12\sqrt{c\log q}\), while
\(\mathbb E_\alpha d\le1\) by Cauchy--Schwarz. The price inequality
gives the result. For bounded \(q\), counting exterior bits suffices.
\end{proof}

In a network for \(\operatorname{Tag}(h,h)^{\otimes t}\), pin each
list's tags to a uniform permutation. The topology realizes a permutation
of \(r=tNh\) payload bits. Each input bit has \(q=N\) possible
destinations, at its own bit position in its own list, all equally
likely. By Lemma~\ref{lem:uniform-demand},
\(2V+O(tNh)\ge c\,tNh\sqrt h\). Since \(n_w=4N(h+1)\), this proves
\eqref{eq:many} for any fixed \(0<\sigma<1/2\) and all sufficiently
large \(w\), uniformly in \(t\).

\subsection{Finite fields with ports}

Fix a finite field with ports. We use the signed counting tensor
\(U_\ell^*\) on \(\ell\) exterior bits from the circuit lower-bound
paper. Its entries count the signed computation paths of a polynomially
clocked universal machine, reduced in the field. The polynomial-size
verifier gives an exact trivalent network of size \(\ell^{O(1)}\) for
one copy. Sum over its guess strings, requiring the unused bits on
each computation path to be zero, and weight the outcome by its sign. The exact batch lower bound
\cite[Theorem~\ext{P}{thm:batches}]{oaP} gives \(\delta,\xi>0\)
such that, for all sufficiently large \(\ell\),
\[
       \operatorname{vol}_k\bigl((U_\ell^*)^{\otimes j}\bigr)
              >j\ell^{1+\delta}
       \qquad (1\le j\le 2^{\ell^\xi}).
\]

The descent with tagged permutations
\cite[Proposition~\ext{G}{prop:oracle-descent}]{oaG} supplies a
constant \(a>0\). For all sufficiently large \(h\), given a trivalent
network of size \(s\) in the range \(N\le s\le N^4\), with \(e\)
exterior bits, it constructs a network for \(d\) copies of its tensor
using trivalent gates and \(\operatorname{Tag}(h,h)\) tests. Its size
\(S\), charging every incident bit including those at tests, satisfies
\[
            S\le C d\bigl(sh^{-a}+e(\log h)^C\bigr),
           \qquad \log_2 d=O\bigl((\log h)^2\bigr).
\]
These bounds are uniform in the network. The normalization required
for descent attaches equality leaves and subdivides self-loops. It
costs only a constant factor, and the cited bound holds on every fixed
polynomial size range.

Choose \(b>\max\{2,2/\xi\}\) and
\(0<\sigma<\min\{a/4,b\delta/4,1/4\}\). Suppose there are arbitrarily
large \(w=h+1\) for which some \(1\le t\le 2N\) admits a realization
of \(D_w^{\otimes t}\) of volume at most \(tn_w w^\sigma\).
The two-list construction gives a network for \(t\) tagged-permutation
tests of size \(O(tNh^{1+\sigma})\).

Group the tests in a descent network into batches of \(t\) and
substitute this realization for each batch. In the last batch pin
unused tests to the identity permutation with zero payloads, where
each has value one. Since the substituted network realizes the full
product tensor, the replacement stays exact under every contraction
of the surrounding wires, including wires joining tests in one batch.
Each full batch contributes \(3tNh\) incident bits to \(S\).
Paying separately for the last batch therefore gives trivalent size
\[
                     O\bigl(h^\sigma S+tNh^{1+\sigma}\bigr).
\]

Set \(\ell=\lceil h^b\rceil\) and begin with \(j=N^2\) copies of
the polynomial-size network for \(U_\ell^*\). At each descent let
\(x=s/(j\ell)\) be the size per exterior bit. The new batch has \(dj\)
copies. Its size per bit, for all sufficiently large \(h\), satisfies
\[
 \begin{split}
 x'&\le C'\left(h^{-a+\sigma}x+
        h^\sigma(\log h)^C+\frac{tNh^{1+\sigma}}{dj\ell}\right)\\
   &\le h^{-a/2}x+h^{2\sigma}.
 \end{split}
\]
Here \(j\ge N^2\), \(d\ge1\), and \(t\le 2N\), so the fraction is
\(O(h^{1+\sigma-b})\). This is the reason for starting with \(N^2\)
copies.

Initially \(x\le h^A\) for some fixed \(A\). Choose a fixed integer
\(L>2A/a\). At every one of these \(L\) steps, the product of the
replication factors is \(N^{o(1)}\). Inductively, \(j=N^{2+o(1)}\)
and \(x\) stays polynomial in \(h\), so the current network and its
normalization have size in the descent range. Iterating the last
inequality gives \(x=O(h^{2\sigma})\). Now \(j\le 2^{\ell^\xi}\)
because \(\log_2 j=(2+o(1))h\), while
\(h^{2\sigma}=o(\ell^\delta)\) by the choice of \(\sigma\).
The final trivalent volume is therefore \(o(j\ell^{1+\delta})\),
contradicting the exact batch lower bound. This proves
Theorem~\ref{thm:many} in the remaining case.
~~~

## Bibliographic attribution

- Anthropic: {Claude (Anthropic)}. A {$2^{o(n)}$} time algorithm for {3-SAT}. 2026.
- BH: Martin R. Bridson and Andr{\'e} Haefliger. Metric Spaces of Non-Positive Curvature. 1999.
- Batcher: Kenneth E. Batcher. Sorting networks and their applications. 1968.
- Benes: V. E. Bene{\v s}. Optimal rearrangeable multistage connecting networks. 1964.
- CVZ: Matthias Christandl and P{\'e}ter Vrana and Jeroen Zuiddam. Universal points in the asymptotic spectrum of tensors. 2023.
- Cook: Stephen A. Cook. Short Propositional Formulas Represent Nondeterministic Computations. 1988.
- Fritz: Tobias Fritz. A generalization of {Strassen}'s {Positivstellensatz}. 2021.
- GapP: Stephen A. Fenner and Lance J. Fortnow and Stuart A. Kurtz. Gap-definable counting classes. 1994.
- HLW: Shlomo Hoory and Nathan Linial and Avi Wigderson. Expander graphs and their applications. 2006.
- KacWang: Victor G. Kac and Weiqiang Wang. Vertex operator superalgebras and their representations. 1994.
- Khoruzhii: Kirill Khoruzhii. Crossing decomposition for subquadratic noncommutative circuits. 2026.
- Leighton: Leighton, Tom. Tight bounds on the complexity of parallel sorting. 1985.
- MendelNaor: Manor Mendel and Assaf Naor. Nonlinear spectral calculus and super-expanders. 2014.
- MendelNaorBarycentric: Mendel, Manor and Naor, Assaf. Spectral calculus and {Lipschitz} extension for barycentric metric spaces. 2013.
- Strassen1988: Volker Strassen. The asymptotic spectrum of tensors. 1988.
- Sturm: Karl-Theodor Sturm. Probability measures on metric spaces of nonpositive curvature. 2003.
- TinyRAM: Eli Ben-Sasson and Alessandro Chiesa and Daniel Genkin and Eran Tromer and Madars Virza. {SNARKs} for {C}: verifying program executions succinctly and in zero knowledge. 2013.
- WZ: Avi Wigderson and Jeroen Zuiddam. Asymptotic spectra: theory, applications, and extensions. 2026.
- ct: Thomas M. Cover and Joy A. Thomas. Elements of Information Theory. 2006.
- falseeth: {OpenAI}. {The Exponential Time Hypothesis is False}. 2026.
- oaG: {OpenAI}. {Euclidean geometry in the asymptotic spectrum of tensors}. 2026.
- oaP: {OpenAI}. {Circuit lower bounds for Hamiltonian cycle and the permanent}. 2026.
- oaT: {OpenAI}. {The logarithmic exponent of integer multiplication on tapes}. 2026.
- odon: Ryan O'Donnell. Analysis of Boolean Functions. 2014.

## License

The repository LICENSE at this commit follows.

~~~text
                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS

   APPENDIX: How to apply the Apache License to your work.

      To apply the Apache License to your work, attach the following
      boilerplate notice, with the fields enclosed by brackets "[]"
      replaced with your own identifying information. (Don't include
      the brackets!)  The text should be enclosed in the appropriate
      comment syntax for the file format. We also recommend that a
      file or class name and description of purpose be included on the
      same "printed page" as the copyright notice for easier
      identification within third-party archives.

   Copyright [yyyy] [name of copyright owner]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
~~~
