# Checked source obligations

This is a literature reproduction. The following statements are paraphrases and
hypothesis checks, not a claim of independent verification of every published
case. The complete proof and family are in proof.md. The original PDFs were
retrieved from the cited author/venue URLs and their SHA-256 digests matched the
starting proof exactly. The PDFs are retained as research artifacts outside the
submission; this package contains text evidence within the verifier limits.

## BK-preprint: construction

Bourgain and Kopparty, author preprint dated 2010-03-12, Theorem 2.2 and equation
(4), printed pages 6-7 (PDF pages 7-8): for input dimension m, choose the least
prime extension degree k>2*m^(3/5), r=ceil(m/k), and t=ceil(sqrt(r)). The theorem
allows any injective F_p-linear map into (F_(p^k))^r and any nonzero F_p-linear
output map. Applying that output map to the t-th elementary symmetric polynomial
of block powers sum_(j=0)^i p^j is nonconstant on every affine subspace of dimension
greater than 6*m^(4/5).

Here p=2, m=n, the chosen polynomial basis and zero padding define the injection,
and the constant-coordinate projection is nonzero. Thus its hypotheses match the
specific uniform family; there is no divisibility assumption or discarded input.

Source: https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
SHA256: 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf

## Shoup-preprint: uniform field construction

Victor Shoup, author preprint dated 1989-01-31, Section 3, Theorem 3.2 and its
immediate fixed-prime consequence, printed page 3: a deterministic algorithm
constructs an irreducible polynomial of a specified exact degree over F_p. For
fixed prime p its running time is polynomial in the degree, without a randomness
or ERH premise. Apply it with fixed p=2 and degree k. Its base-field operations
have constant bit cost. Fix one deterministic implementation for all k; this
provides the uniform construction needed by the evaluator, rather than advice.

The rendered theorem page was inspected because the PDF's legacy font encoding
makes extracted text unreliable. This is a source-reading limitation, not a
mathematical hypothesis. No result about conditional factorization is substituted.

Source: https://www.shoup.net/papers/detirred.pdf
SHA256: 00590fb8f0717338530b067626428781991edb840fb74c37cf846f07f50e0644

## LY-full: circuit theorem and finite loss

Li and Yang, ECCC TR21-023 dated 2021-02-20, Theorem 1.1, printed page 2, gives
coefficient 31/10 for B2 circuits computing an affine disperser of sublinear
dimension. The basis and gate-count conventions are in Sections 1.1 and 2.1;
constants and input complements are addressed by Section 2.4, footnotes 3-4;
normalization, including repeated pins and output gates, is in Section 3.3.
The family above has integer disperser threshold d=floor(6*n^(4/5))+1=o(n).

The conservative finite derivation in proof.md uses Definition 3.6, Corollary 3.8
and Theorem 4.1 (printed pages 14-15 and 18-19). With no quadratic restrictions,
mu<=|C|+alpha_I*n+alpha_phi*(n/2+5*d/2), while mu>=delta*(n-2*d-2).
The chosen parameters are nonnegative, alpha_Q<=alpha_I and alpha_phi<=1/2.
The resulting remainder is (261/10)*d+128/5. The n/2+d expression in Theorem 3.9
is not used to improve this finite bound. The cited theorem supplies the general
circuit analysis; the LP arithmetic by itself supplies no such theorem.

Source: https://eccc.weizmann.ac.il/report/2021/023/download/
SHA256: 4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062
