# Why the same finite-field processor resists cheap parallel realizations

The auxiliary tensor here can change with the width. Its only function is to contradict a compiler which promises to work for **every** input diagram at that width. It need not be uniform in the width or itself decide distinctness.

**The macro table class.** For a table on t incident bit ports use L = Cheap(F,50000,2). Membership requires an ordinary network computing the exact table with at most (1+t)·50000·(1+ceil(log₂(1+t)))² vertices. A wide realization charges 1+t per table and bounds its arity by D. Both the processor below and the second output of compiler.md use this identical L. Trivalent gates are admitted separately. F=F₂ satisfies every field and finiteness instance.

**Fixed processor before the proposed competitor.** Write N=2^w. Fix the target cost C₀≥4, use C=3C₀+12, and choose a bounded expansion for 52C. Its dimension n≥2 and packet height h>0 are constants at this point. Take r=2^floor(w/(16n)), d₀=2^floor(w/(32n)), and l=floor(N/2048). A processor input has a literal of l bits, a public report of l bits and a block of floor(N/4) bits for links. Their total fits N. Exterior labels are chosen using the tape positions alone. With those labels fixed, the control construction supplies U with a pad string satisfying U(pad)=1 and an ordinary network of size at most N·w^k. The choice of U is before the universal quantifier over competitor graphs. Different sites will choose different input literals for this one U.

**Why the literals fit.** At the fixed dimensions the controller uses a bank of fixed template interpreters. Its description and site identifiers occupy at most the l literal bits; one inequality in programmable_of_control allocates half that literal space to the bank. Each wide table has a short network description because its table belongs to L. The controller reads these descriptions on its literal tape instead of hardwiring the competitor into U. Its selector dispatches a count guard, the site computation or a pinned branch. Choose a working address width W with w≤W≤2w and 2^W≤N·w^(Q+10), where Q=q+8 and q is the logarithmic exponent in Cheap. The sorted local evaluation and dispatch incur 2^W(100+W)^(Q+660) work in the bound used by control. For w≥200 this is at most N·w^((Q+10)+2(Q+660)). The constant-size header dispatch and other polynomial factors leave an exponent k fixed before the width. The source steps below include the construction of the window, bank, dispatch, and these inequalities. Finiteness is needed for field entries to admit finite literal codes.

**Geometric input to the layout.** The supply edges of an n-dimensional cube join nearest neighbors, while the demands join antipodes. For a monotone spectral price f, let M be the longest priced neighbor distance, and write V for the number of cube vertices. The square inequality supplies a constant 1<κ<2 valid uniformly over the prices in use. At n=2^j the antipodal squared energy is at most (2κ)^j V M². Cauchy–Schwarz therefore bounds the sum of antipodal distances by V M (2κ)^(j/2). For any prescribed L, increasing j gives L²(2κ)^j<4^j. Thus L times the antipodal sum is smaller than n V M. Connectivity excludes M=0. The resulting strict spectral comparisons yield an exact finite sum of tensor restrictions, uniformly under amplification; see the powered-gadget proof in compiler.md. Multiplying every bundle width by the same fixed integer pays for the selector bits without changing the crossing ratio. The resulting bound on the crossing positions lets prepare use the same fixed slot budget for every selector. The proof excerpts include the quantifiers, strict comparison, restriction to exact packets and packet budget.

**How a cheap competing graph would program itself.** First repeat the competing diagram and pin surplus copies to pad to reach d₀; each pinned copy has scalar one. The scales, including the bound on the original copy count, absorb the replication and pinning into C. Give its exterior pins the previously chosen labels. The layout in prepare has at most 13CNd₀ positions; replicating it B=ceil(3r/d₀) times uses at most 52CrN positions since Bd₀≤4r. The bounded expansion therefore assigns each worker at most N/128 packet slots in its stated budget. This leaves room for literal and physical links. Reports encode the cut and selector counts for each class of packets, with an injective code and one extra reserved symbol. Invalid reports or codes contribute zero. For a valid report, the layout serializes each shared link once and each worker evaluates the local factors corresponding to that report. Summing the products over those shared bit assignments recovers the original evaluation (route, unpinned and replay below). The prescribed literal at each site is the description that makes U implement this local calculation.

**The contradiction is a scalar identity, including in characteristic two.** On the extra code, each physical link is pinned to zero. There is exactly one compatible global assignment because every link is visible to a worker. This adds exactly one to the sum. On an ordinary encoded report the expression is the product of the programmed site evaluations. All unused codes vanish. Thus if the diagram really computes the claimed batch, replay equates a scalar x with 1+x: the ordinary report sum is identical on the two sides, but the programmed processor has the extra unit branch. This would imply 1=0 in F. As the scale inequalities hold eventually and allow every 1≤d≤N^δ and D²≤N for one fixed δ>0, every such wide realization costs strictly more than C₀Nd. The finite-field and arity hypotheses used here are precisely those supplied in compiler.md.

## Checked source steps

Each excerpt below is a complete named declaration from the pinned source. The surrounding namespaces and imports are identified by the file path. The excerpts are for reading alongside the derivation. The built originals and their imported lemmas supply the elaboration environment. Attribution and the applicable Apache 2.0 license are in verification.md.

### Processor/Labels.lean

Lines 133–133

```lean
def arity (T : Table F) : ℕ := wd T.1 + wd T.2.1
```

Lines 135–136

```lean
def Cheap (F : Type) [Field F] (K q : ℕ) : Set (Table F) :=
  {T | ∃ m ≤ (1 + arity T) * K * (1 + Nat.clog 2 (1 + arity T)) ^ q, Net F T.1 T.2.1 T.2.2 m}
```

### Diagonal/Bounds.lean

Lines 12–18

```lean
def Resists (F : Type) [Field F] (L : Set (Table F)) : Prop :=
  ∀ C : ℕ, 4 ≤ C → ∃ k : ℕ, ∃ δ : ℝ, 0 < δ ∧ ∀ᶠ w : ℕ in atTop,
    let N : ℕ := 2 ^ w
    ∃ (a : Wire) (U : Val a → Val .unit → F) (v : ℕ),
      wd a ≤ N ∧ Net F a .unit U v ∧ v ≤ N * w ^ k ∧
      ∀ d D m : ℕ, 1 ≤ d → (d : ℝ) ≤ (N : ℝ) ^ δ → D ^ 2 ≤ N →
        Wide F D L (copies a d) (copies .unit d) (power U d) m → C * N * d < m
```

### Spectral/CubeGadget.lean

Lines 352–359

```lean
private theorem cube_dimension {κ : ℝ} (h1 : 1 < κ) (h2 : κ < 2) (L k₀ : ℕ) :
    ∃ k ≥ k₀, (L : ℝ) ^ 2 * (2 * κ) ^ k < ((2 : ℝ) ^ k) ^ 2 := by
  have hr : 1 < 4 / (2 * κ) := (one_lt_div (by positivity)).mpr (by linarith)
  obtain ⟨j, hj⟩ := pow_unbounded_of_one_lt ((L : ℝ) ^ 2) hr
  refine ⟨max j k₀, le_max_right _ _, ?_⟩
  have hg := hj.trans_le (pow_le_pow_right₀ hr.le (le_max_left j k₀))
  rw [div_pow, lt_div_iff₀ (by positivity)] at hg
  simpa only [show (4 : ℝ) = 2 ^ 2 by norm_num, ← pow_mul, Nat.mul_comm] using hg
```

Lines 305–350

```lean
private theorem cube_gap {κ : ℝ} (hk : 0 ≤ κ) (k L : ℕ)
    (hs : (L : ℝ) ^ 2 * (2 * κ) ^ k < ((2 : ℝ) ^ k) ^ 2)
    (f : Tensor K (Cube (2 ^ k)) →+* ℝ) (hf : Monotone f)
    (hq : ∀ a b c e, (dist f a c) ^ 2 + (dist f b e) ^ 2 ≤
      κ * ((dist f a b) ^ 2 + (dist f b c) ^ 2 + (dist f c e) ^ 2 + (dist f e a) ^ 2)) :
    ∃ i, f (cubeDemand K (2 ^ k) ^ L) < f (cubeSupply K (2 ^ k) i) ^ (2 ^ k * Fintype.card (Cube (2 ^ k))) := by
  let n := 2 ^ k
  let N := Fintype.card (Cube n)
  let D := ∑ x, dist f x (Cube.flip Finset.univ x)
  let w (i : Cube n × Fin n) := dist f i.1 (Cube.flip {i.2} i.1)
  have hN : (0 : ℝ) < N := by exact_mod_cast (Fintype.card_pos : 0 < Fintype.card (Cube n))
  let : Nonempty (Fin n) := Fin.pos_iff_nonempty.mp (by unfold n; positivity)
  obtain ⟨i, _, hi⟩ := Finset.exists_max_image Finset.univ w Finset.univ_nonempty
  have he (j) : w j ≤ w i := hi j (Finset.mem_univ _)
  have hw : 0 < w i := by
    by_contra! bad
    have hz (j : Cube n × Fin n) : dist f j.1 (Cube.flip {j.2} j.1) = 0 :=
      le_antisymm ((he j).trans bad) (dist_nonneg f hf _ _)
    have hp (j : Cube n × Fin n) : f (wire K j.1 (Cube.flip {j.2} j.1) 2) = 1 := by
      rw [price_wire f hf _ _ (by decide), hz, Real.rpow_zero]
    have h := hf (cube_connected K n)
    simp only [map_ofNat, map_prod, hp, Finset.prod_const_one] at h
    norm_num at h
  have hbound := Cube.energy_bound (dist f) hq hk (m := (w i) ^ 2)
    (fun x j => pow_le_pow_left₀ (dist_nonneg f hf _ _) (he (x, j)) _) k Finset.univ (by simp)
  have hcauchy : D ^ 2 ≤ N * Cube.energy (dist f) Finset.univ := by
    simpa [D, N, Cube.energy] using Finset.sum_mul_sq_le_sq_mul_sq Finset.univ
      (fun _ : Cube n => (1 : ℝ)) (fun x => dist f x (Cube.flip Finset.univ x))
  have h : ((L : ℝ) * D) ^ 2 < ((n : ℝ) * N * w i) ^ 2 := calc
    _ ≤ (L : ℝ) ^ 2 * (N * ((2 * κ) ^ k * N * (w i) ^ 2)) := by
      rw [mul_pow]
      exact mul_le_mul_of_nonneg_left (hcauchy.trans (mul_le_mul_of_nonneg_left hbound hN.le))
        (sq_nonneg _)
    _ = ((L : ℝ) ^ 2 * (2 * κ) ^ k) * ((N : ℝ) * w i) ^ 2 := by ring
    _ < ((2 : ℝ) ^ k) ^ 2 * ((N : ℝ) * w i) ^ 2 := mul_lt_mul_of_pos_right hs (by positivity)
    _ = _ := by unfold n; push_cast; ring
  have hl : (L : ℝ) * D < (n : ℝ) * N * w i := by
    by_contra! bad
    exact h.not_ge (pow_le_pow_left₀ (by positivity) bad _)
  refine ⟨i, ?_⟩
  simp only [cubeDemand, cubeSupply, map_pow, map_prod, price_wire f hf _ _ (n := 2) (by decide),
    Nat.cast_ofNat, ← Real.rpow_sum_of_pos (by norm_num : (0 : ℝ) < 2),
    ← Real.rpow_mul_natCast (by norm_num : (0 : ℝ) ≤ 2)]
  apply (Real.rpow_lt_rpow_left_iff (by norm_num : (1 : ℝ) < 2)).mpr
  rw [Nat.cast_mul]
  exact (mul_comm D _).trans_lt (hl.trans_eq (mul_comm _ _))
```

Lines 361–386

```lean
theorem cube_site (L k₀ : ℕ) (hL : 0 < L) :
    ∃ n ≥ k₀, ∃ (b C : ℕ) (H : Finset (Cube n × Fin n → ℕ)), 0 < b ∧
      (∀ h ∈ H, L * ∑ i, h i < n * Fintype.card (Cube n) * b) ∧
      ∀ z : ℕ, 0 < z → cubeDemand K n ^ (b * z) ≤
        (C : Tensor K (Cube n)) * ∑ h ∈ H, ∏ i, cubeSupply K n i ^ (h i * z) := by
  obtain ⟨κ, h1, h2, hq⟩ := square K
  obtain ⟨k, hk, hs⟩ := cube_dimension h1 h2 L k₀
  let n := 2 ^ k
  have hn : 0 < n := by unfold n; positivity
  let : Nonempty (Fin n) := Fin.pos_iff_nonempty.mp hn
  have hsize : 0 < n * Fintype.card (Cube n) := Nat.mul_pos hn Fintype.card_pos
  obtain ⟨b, C, H, hb, hcost, hcomp⟩ := Spectrum.powered_gadget
    add_le_add_right mul_le_mul_left zero_le one_le exists_bound
    (cubeDemand K n ^ L) (cubeSupply K n) _ hsize
    (fun _ => one_le_wire K _ _ (by decide)) (by
      intro g hg
      let f := NNReal.toRealHom.comp g
      have hf : Monotone f := fun _ _ h => hg h
      obtain ⟨i, hi⟩ := cube_gap K (by linarith) k L hs f hf (hq _ f hf)
      exact ⟨i, NNReal.coe_lt_coe.mp hi⟩) (cube_connected K n)
  refine ⟨n, hk.trans (Nat.lt_two_pow_self (n := k)).le, L * b, C, H, Nat.mul_pos hL hb, ?_, ?_⟩
  · intro h hh
    exact (Nat.mul_lt_mul_of_pos_left (hcost h hh) hL).trans_eq (by ring)
  · intro z hz
    rw [Nat.mul_assoc, pow_mul]
    exact hcomp z hz
```

### Spectral/Packets.lean

Lines 16–28

```lean
private theorem expand_profile {S I J : Type} [Fintype S] [DecidableEq S] [Nonempty S]
    [Fintype I] [Fintype J] (A : I → S × S) (hA : ∀ i, (A i).1 ≠ (A i).2)
    (B : J → S × S) (w C : ℕ) (H : Finset (J → ℕ)) (z : ℕ)
    (h : (∏ i, wire K (A i).1 (A i).2 2) ^ w ≤
      (C : Tensor K S) * ∑ q ∈ H, ∏ j, wire K (B j).1 (B j).2 2 ^ (q j * z)) :
    Packet K (fun e : I × Fin w => A e.1)
      (fun i : Fin C × ↥H => fun e : Σ j : J, Fin (i.2.1 j * z) => B e.1) := by
  classical
  apply packet_of_restriction _ (fun e => hA e.1)
  simp only [Fintype.prod_prod_type, Fintype.prod_sigma, Finset.prod_const, Finset.card_univ,
    Fintype.card_fin, Finset.prod_pow, Fintype.sum_prod_type, Finset.sum_const, nsmul_eq_mul]
  rw [Finset.sum_coe_sort H (fun q => ∏ j, wire K (B j).1 (B j).2 2 ^ (q j * z))]
  exact h
```

Lines 43–57

```lean
theorem cube_packets (L k₀ : ℕ) (hL : 0 < L) :
    ∃ n ≥ k₀, ∃ (b C : ℕ) (H : Finset (Cube n × Fin n → ℕ)), 0 < b ∧
      (∀ h ∈ H, L * ∑ i, h i < n * Fintype.card (Cube n) * b) ∧ ∀ z : ℕ, 0 < z →
        Packet K (fun e : Cube n × Fin (b * z) => (e.1, Cube.flip Finset.univ e.1))
          (fun i : Fin C × ↥H => fun e : Σ j : Cube n × Fin n, Fin (i.2.1 j * z) =>
            (e.1.1, Cube.flip {e.1.2} e.1.1)) := by
  obtain ⟨n, hn, b, C, H, hb, hcost, h⟩ := cube_site K L (k₀ + 1) hL
  refine ⟨n, by omega, b, C, H, hb, hcost, fun z hz => ?_⟩
  have ha (x : Cube n) : x ≠ Cube.flip Finset.univ x := by
    intro he
    have := congrFun he (⟨0, by omega⟩ : Fin n)
    simp [Cube.flip] at this
  refine expand_profile K (fun x : Cube n => (x, Cube.flip Finset.univ x)) ha
    (fun i : Cube n × Fin n => (i.1, Cube.flip {i.2} i.1)) (b * z) C H z ?_
  exact h z hz
```

### Diagonal/AlphabetCode.lean

Lines 471–519

```lean
theorem expand_packet {I : Type} {A : I → Type} [Fintype I] [∀ s, Fintype (A s)]
    (edges : (s : I) → A s → Cube n × Fin n)
    (packet : Spectral.Packet F (fun p : Cube n × Fin h => (p.1, Cube.flip Finset.univ p.1))
      (fun s e => ((edges s e).1, Cube.flip {(edges s e).2} (edges s e).1))) :
    ∃ E : Expansion F n h, Fintype.card E.I = Fintype.card I ∧
      ∀ s, ∃ t : I, ∑ i, E.cross s i = Fintype.card (A t) := by
  obtain ⟨f, hf, full⟩ := packet
  let κ (s : I) (i : Fin n) (e : A s) := region (edges s e) i
  let φ (s : I) (i : Fin n) (x y : Block n h) (b : Bool) (u : A s → Bool) :=
    half i (fun v => f s v x y u) b
  let E : Expansion F n h := {
    I := I
    weight := fun s x y => ∑ u, ∏ v, f s v x y u
    cross := fun s i => Fintype.card {e // κ s i e = none}
    side := fun s i b x y z => Cut.part (κ s i) (φ s i x y) b z
    full := fun x y => (full x y).symm
    split := by
      intro s i x y
      simp_rw [prod_half i]
      apply Cut.split (κ s i) (φ s i x y)
      intro b u u' he
      apply half_congr
      intro v hv
      refine hf _ _ _ _ _ _ _ _ (fun _ _ => rfl) (fun _ _ => rfl) (fun e hh => he e ?_)
      by_cases hc : (edges s e).2 = i
      · exact .inl (by simp [κ, region, hc])
      · right
        rcases hh with h₁ | h₂
        · change (edges s e).1 = v at h₁
          simp [κ, region, hc, h₁, hv]
        · have h₃ := (congrFun h₂ i).trans hv
          simpa [κ, region, Cube.flip, hc, Ne.symm hc] using h₃
    depends := by
      intro s i b x x' y y' z hx hy
      apply Finset.sum_congr rfl
      intro u _
      apply half_congr
      intro v hv
      apply hf
      · intro e he
        exact hx e ((congrFun he i).trans hv)
      · intro e he
        apply hy
        have h₁ := (congrFun he i).trans hv
        simp [Cube.flip] at h₁
        grind
      · intro e _; rfl
  }
  exact ⟨E, rfl, fun s => ⟨s, region_card (edges s)⟩⟩
```

Lines 524–540

```lean
theorem expansion (F : Type) [Field F] (L n₀ : ℕ) (hL : 0 < L) :
    ∃ n ≥ n₀, ∃ h : ℕ, ∃ E : Expansion F n h, 0 < h ∧
      L * Nat.clog 2 (Fintype.card E.I) ≤ h ∧
      ∀ s, L * ∑ i, E.cross s i < n * Fintype.card (Cube n) * h := by
  obtain ⟨n, hn, b, C, H, hb, hc, hp⟩ := Spectral.Tensor.cube_packets F L n₀ hL
  let z := L * Nat.clog 2 (Fintype.card (Fin C × ↥H)) + 1
  have hz : 0 < z := Nat.zero_lt_succ _
  let edges (s : Fin C × ↥H) (e : Σ j : Cube n × Fin n, Fin (s.2.val j * z)) := e.1
  obtain ⟨E, he, hcost⟩ := expand_packet edges (hp z hz)
  refine ⟨n, hn, b * z, E, Nat.mul_pos hb hz, ?_, ?_⟩
  · rw [he]
    exact (Nat.le_add_right _ _).trans (Nat.le_mul_of_pos_left z hb)
  · intro s
    obtain ⟨t, ht⟩ := hcost s
    rw [ht, Fintype.card_sigma]
    simp only [Fintype.card_fin, ← Finset.sum_mul, ← mul_assoc]
    exact Nat.mul_lt_mul_of_pos_right (hc t.2.val t.2.property) hz
```

Lines 584–602

```lean
private theorem plan_bound (hn : 2 ≤ n) (N C L : ℕ) (s : Packet n r d → A.I)
    (hN : 0 < N) (hL : 0 < L) (hsave : 32 * 256 ^ 2 * C ^ 2 ≤ L)
    (hw : Fintype.card (Position n r d h) ≤ C * r * N)
    (hc : ∀ s, L * ∑ i, A.cross s i ≤ n * Fintype.card (Cube n) * h)
    (hl : L * Nat.clog 2 (Fintype.card A.I) ≤ h) (u : Worker r) :
    128 * ∑ p, (A.plan hn N s).slotsAt h u p ≤ N := by
  let π := A.plan hn N s
  obtain ⟨i, hi, hm⟩ := A.coordinate (by omega) L hc hl s
  have he : ∑ p, A.cross (s p) π.cut = ∑ p, A.cross (s p) i := by
    refine le_antisymm ?_ (hi _)
    exact (first_spec hn (fun j => ∑ p, A.cross (s p) j) i).elim
      (fun h => h.le) (fun h => h.1.le)
  have hb : L * ∑ p, π.q p ≤ 2 * Fintype.card (Position n r d h) := by
    change L * (∑ p, (A.cross (s p) π.cut + Nat.clog 2 (Fintype.card A.I))) ≤ _
    rw [Finset.sum_add_distrib, he, ← Finset.sum_add_distrib]
    simpa only [Position, Fintype.card_prod, Fintype.card_fin, mul_assoc] using hm
  have := π.capacity hN hL hsave hw (other_ne hn π.cut) rfl hb u
  change 128 * ∑ p, π.slotsAt h u p ≤ N
  omega
```

Lines 609–618

```lean
theorem bounded_expansion (F : Type) [Field F] (C : ℕ) :
    ∃ n, ∃ hn : 2 ≤ n, ∃ h, ∃ A : Expansion F n h, 0 < h ∧
      ∀ (r d N : ℕ), 0 < N → Fintype.card (Position n r d h) ≤ C * r * N →
        ∀ (s : Packet n r d → A.I) (u : Worker r),
          128 * ∑ p, (A.plan hn N s).slotsAt h u p ≤ N := by
  let L := 32 * 256 ^ 2 * C ^ 2 + 1
  have hl : 0 < L := Nat.zero_lt_succ _
  obtain ⟨n, hn, h, A, hh, hs, hc⟩ := expansion F L 2 hl
  exact ⟨n, hn, h, A, hh, fun _ _ N hN hw s u =>
    A.plan_bound hn N C L s hN hl (Nat.le_succ _) hw (fun t => (hc t).le) hs u⟩
```

### Diagonal/Broadcast.lean

Lines 431–464

```lean
theorem prepare (A : Expansion F n h) (hn : 2 ≤ n) (hh : 0 < h) (hr : 0 < r) (hC : 4 ≤ C)
    (hA : ∀ d', Fintype.card (Position n r d' h) ≤ 52 * C * r * N →
      ∀ (s : Packet n r d' → A.I) (u : Worker r), 128 * ∑ p, (A.plan hn N s).slotsAt h u p ≤ N)
    (fa : Ports N l b → Label n r) (hfa : ∀ g, multiplicity fa g * r ^ n ≤ 4 * N)
    (hlit : ∀ i g, multiplicity (fun t : Fin l => fa (.inl t) i) g * r ≤ 2 * l)
    (hl : 1024 * l ≤ N) (hlink : 64 * (3 * r) ^ 2 ≤ N) (hB : B * d ≤ 4 * r)
    (G : MachineDiagram F N l b d) (hc : G.cost ≤ C * N * d) (hd : ∀ v, G.arity v ≤ D)
    (hsize : (D + 1) * (r ^ n) ^ 4 ≤ C * N * d) (hpad : (r ^ n) ^ 2 * 2 ^ n * h ≤ C * N * d) :
    ∃ S : Scheme A hn G r B, S.colors = fa ∧
      S.layout.d = rounds n r h (12 * (C * N * d)) ∧
      Fintype.card (Position n r S.layout.d h) ≤ 13 * C * N * d ∧
      ∀ z, multiplicity S.layout.color z * r ^ n ≤ 16 * C * N * d := by
  let f := fun s => fa (address N l b (port (tape N l b) d s).2)
  let g : slots (copies .unit d) → Label n r := fun s => (port .unit d s).2.elim
  have hf (z) : multiplicity f z * r ^ n ≤ C * N * d :=
    (exterior_copies fa hfa d z).trans (by gcongr)
  have hg (z) : multiplicity g z * r ^ n ≤ C * N * d := by
    have : IsEmpty (slots (copies .unit d)) := ⟨fun s => (port .unit d s).2.elim⟩
    simp [multiplicity]
  obtain ⟨L, heq, hpos, hf', _, hv⟩ := layout G f g hc hd hr hh hf hg hsize hpad
  have hall : Fintype.card (Position n r (B * L.d) h) ≤ 52 * C * r * N := by
    calc
      _ = B * Fintype.card (Position n r L.d h) := by
        simp only [Position, Packet, Fintype.card_prod, Fintype.card_fin]
        ring
      _ ≤ B * (13 * (C * N * d)) := Nat.mul_le_mul_left B hpos
      _ = 13 * C * N * (B * d) := by ring
      _ ≤ 13 * C * N * (4 * r) := Nat.mul_le_mul_left _ hB
      _ = _ := by ring
  let S : Scheme A hn G r B :=
    ⟨hlink, L, fa, hf', fun s => A.report_budget hn N s (hA _ hall s),
      fun R u => Pin.budget A hn R fa hlit hl u⟩
  exact ⟨S, rfl, heq, by simpa only [S, mul_assoc] using hpos,
    fun z => by simpa only [S, mul_assoc] using hv z⟩
```

### Diagonal/Replay.lean

Lines 110–113

```lean
def branch (code : Option Q ↪ (Fin b → Bool)) (run : Q → (Fin (N / 4) → Bool) → F)
    (z : Fin b → Bool) (y : Fin (N / 4) → Bool) : F :=
  Function.extend code (fun q => q.elim (if y = (fun _ => false) then 1 else 0) (fun t => run t y))
    (fun _ => 0) z
```

Lines 115–131

```lean
theorem branch_sum (hR : 0 < R) (code : Option Q ↪ (Fin b → Bool))
    (run : Fin R → Q → (Fin (N / 4) → Bool) → F) :
    (∑ z, ∑ y, ∏ u, branch code (run u) z (Links.view u y)) =
      1 + ∑ q, ∑ y, ∏ u, run u q (Links.view u y) := by
  rw [← sum_subset (subset_univ (univ.map code)), sum_map]
  · simp only [branch, code.injective.extend_apply]
    rw [Fintype.sum_option]
    simp only [Option.elim_none, Option.elim_some, prod_ite_zero, mem_univ, forall_const,
      prod_const_one, zero_iff, sum_ite_eq', ite_true]
  · intro z _ h
    have hz : ¬∃ q, code q = z := by
      rintro ⟨q, rfl⟩
      exact h (mem_map.mpr ⟨q, mem_univ q, rfl⟩)
    apply sum_eq_zero
    intro y _
    apply prod_eq_zero (mem_univ (⟨0, hR⟩ : Fin R))
    exact Function.extend_apply' _ _ _ hz
```

Lines 830–832

```lean
def simulate (serial : S.Serial R) (u : Fin (3 * r)) (ξ : Fin l → Bool) (bits : Fin (N / 4) → Bool) : F :=
  A.extraProcessor hn N R.1 (Worker.index r) R.2.2 S.links
    (Pin.src (S.toSite R)) (Pin.dst (S.toSite R)) (S.pins R.1) (S.joint pad R ξ) serial u bits
```

Lines 853–864

```lean
private theorem route (R : S.Report) (serial : S.Serial R) :
    (∑ y, ∑ a : A.Assignment hn N R.1, ∏ u, S.joint pad R (ξ (Worker.index r u)) u y a) =
      ∑ t, ∏ i, S.simulate pad R serial i (ξ i) (Links.view i t) := by
  let g := fun u => S.joint pad R (ξ (Worker.index r u)) u
  have he := A.physical_extra hn N R.1 (Worker.index r) R.2.2 S.links
    (Pin.src (S.toSite R)) (Pin.dst (S.toSite R)) (S.pins R.1) g
    (fun u => S.joint_local pad R _ u) serial
  refine he.trans (sum_congr rfl (fun t _ => prod_congr rfl (fun i _ => ?_)))
  unfold simulate
  apply congrFun
  apply A.extraProcessor_congr
  simp only [g, Equiv.apply_symm_apply]
```

Lines 866–877

```lean
private theorem unpinned (serial : (R : S.Report) → S.Serial R) :
    (∑ z, ∑ y, ∏ j : Fin B, G.eval (inputs pad ξ z y j) (unitCopies d)) =
      ∑ R, ∑ t, ∏ i, S.simulate pad R (serial R) i (ξ i) (Links.view i t) := by
  rw [← Pin.sum_masks ξ (fun s z y => ∏ j : Fin B, G.eval (inputs pad s z y j) (unitCopies d))]
  have he (y : Pin N l b r → Bool) :=
    S.layout.virtual A hn N (A.budget hn N) S.packets (Pin.input pad y) (fun _ => unitCopies d)
  simp only [Pin.input] at he
  simp_rw [he, mul_sum, ← prod_mul_distrib]
  rw [sum_comm]
  apply sum_congr rfl
  intro R _
  exact S.route pad ξ R (serial R)
```

Lines 881–884

```lean
theorem replay (U : Val (tape N l b) → Val .unit → F) (serial : (R : S.Report) → S.Serial R)
    (hR : 3 * r ≤ B * d) (hx : U pad () = 1) (hG : G.eval = power U d) :
    scalar U ξ = ∑ R, ∑ t, ∏ i, S.simulate pad R (serial R) i (ξ i) (Links.view i t) :=
  (substitute U ξ pad G hR hx hG).trans (S.unpinned pad ξ serial)
```

### Diagonal/Broadcast.lean

Lines 481–483

```lean
def Programs (code : Option S.Report ↪ (Fin b → Bool)) : Prop :=
  ∃ serial : (R : S.Report) → S.Serial R, ∀ i : Fin (3 * r), ∃ ξ : Fin l → Bool,
    ∀ z y, U (feed ξ z y) () = branch code (fun R => S.simulate pad R (serial R) i ξ) z y
```

Lines 487–498

```lean
theorem diagonal (code : Option S.Report ↪ (Fin b → Bool))
    (hr : 0 < r) (hR : 3 * r ≤ B * d) (hx : U pad () = 1)
    (hU : S.Programs pad U code) : G.eval ≠ power U d := by
  intro hG
  obtain ⟨serial, hprog⟩ := hU
  choose ξ hξ using hprog
  have h := S.replay pad ξ U serial hR hx hG
  unfold scalar at h
  simp_rw [hξ] at h
  rw [branch_sum (by omega : 0 < 3 * r) code
    (fun i R => S.simulate pad R (serial R) i (ξ i))] at h
  exact one_ne_zero (add_eq_right.mp h)
```

Lines 705–721

```lean
theorem repeat_trim {N C e : ℕ} (hd : 1 ≤ d) (he : d ≤ 2 * e) (ha : wd a ≤ N)
    (h : Wide F D L (copies a d) (copies .unit d) (power U d) m) (hm : m ≤ C * N * d) :
    ∃ v ≤ (3 * C + 12) * N * e, Wide F (max D 3) L (copies a e) (copies .unit e) (power U e) v := by
  let k := e / d + 1
  have hcover : e ≤ d * k ∧ d * k ≤ 3 * e := by
    have := Nat.mod_lt e hd
    have := Nat.mod_add_div e d
    simp only [k, mul_add, mul_one]
    omega
  have hh := join (repeat_wide (h.mono (le_max_left D 3) (Set.Subset.refl _)) k)
  obtain ⟨v, hv, hnet⟩ := trim pad hp (le_max_right _ _) hcover.1 hh
  refine ⟨v, hv.trans ?_, hnet⟩
  calc
    _ ≤ k * (C * N * d) + 4 * (d * k) * N := by gcongr
    _ = (C + 4) * N * (d * k) := by ring
    _ ≤ (C + 4) * N * (3 * e) := Nat.mul_le_mul_left _ hcover.2
    _ = _ := by ring
```

Lines 725–743

```lean
theorem batch_floor (w n : ℕ) :
    let m := 2 ^ (w / (32 * n))
    let r := 2 ^ (w / (16 * n))
    0 < m ∧ m ^ 2 ≤ r ∧ ∀ d : ℕ, d ^ 2 ≤ r → d ≤ 2 * m := by
  have h : w / (32 * n) = w / (16 * n) / 2 := by
    rw [Nat.div_div_eq_div_mul, show 16 * n * 2 = 32 * n by ring]
  rw [h]
  generalize w / (16 * n) = j
  intro m r
  refine ⟨by unfold m; positivity, ?_, ?_⟩
  · unfold m r
    rw [← pow_mul]
    exact Nat.pow_le_pow_right (by decide) (by omega)
  · have hr : r ≤ 2 * m ^ 2 := by
      unfold m r
      rw [← pow_mul, ← pow_succ']
      exact Nat.pow_le_pow_right (by decide) (by omega)
    intro d hd
    nlinarith
```

### Processor/WindowSize.lean

Lines 226–263

```lean
theorem window (n h C Q : ℕ) (hn : 2 ≤ n) (hh : 0 < h) (hC : 4 ≤ C)
    (A : Expansion F n h) :
    ∀ᶠ w : ℕ in atTop,
      let N := 2 ^ w
      let r := 2 ^ (w / (16 * n))
      let d := 2 ^ (w / (32 * n))
      let l := N / 2048
      let B := (3 * r + d - 1) / d
      ∀ fa : Machine.Ports N l l → Label n r,
        let P := Format.prescribed d h (rounds n r h (12 * (C * N * d))) fa
        ∀ M : Bank (Fin n), M.budget F ≤ N * w ^ Q → Nonempty (Controller P A hn M B) →
          ∃ W, Window P A M B C Q w W := by
  obtain ⟨δ, _, hs⟩ := Diagonal.scales C n h (Fintype.card A.I) hC hn
  filter_upwards [hs, eventually_ge_atTop (threshold A C), eventually_ge_atTop (32 * n),
    poly_bound (Q + 4) 1 (by decide)] with w hs hT hw he
  intro N r d l B fa P M hM hS
  obtain ⟨hd, hd', _⟩ := batch_floor w n
  have sc : Scales C n h (Fintype.card A.I) N r l l d 0 B :=
    hs.2 d 0 hd hd' (by simp)
  have hr : 3 ≤ r := by
    have h₁ : 2 ≤ w / (16 * n) := (Nat.le_div_iff_mul_le (by omega)).mpr (by omega)
    exact (show 3 ≤ 2 ^ 2 by decide).trans (Nat.pow_le_pow_right (by decide) h₁)
  have hp : r ^ (16 * n) ≤ N := by
    unfold r N
    rw [← pow_mul]
    exact Nat.pow_le_pow_right (by decide) (Nat.div_mul_le_self w (16 * n))
  have hw' : 2 ≤ w := by omega
  let j := Nat.clog 2 (w ^ (Q + 4))
  have hj : j ≤ w := Nat.clog_le_of_le_pow (by simpa only [Nat.div_one] using he)
  have hlog := Nat.clog_le_of_le_pow (Nat.lt_pow_succ_log_self (by decide : 1 < 2) (w ^ (Q + 4))).le
  have hhigh : 2 ^ j ≤ w ^ (Q + 10) := calc
    _ ≤ 2 ^ (Nat.log 2 (w ^ (Q + 4)) + 1) := Nat.pow_le_pow_right (by decide) hlog
    _ ≤ w ^ (Q + 4) * w := by rw [pow_succ]; exact Nat.mul_le_mul (Nat.pow_log_le_self _ (by positivity)) hw'
    _ ≤ _ := by rw [← pow_succ]; exact Nat.pow_le_pow_right (by omega) (by omega)
  refine ⟨w + j, finite_window P A hn M hS.some hh rfl hr hp (Nat.div_le_self _ _)
    (Window.position_bound sc.padding) sc.copies hM hT (by omega) (by omega) ?_ ?_⟩
  · rw [pow_add]; exact Nat.mul_le_mul_left N hhigh
  · rw [pow_add (2 : ℕ)]; exact Nat.mul_le_mul_left N (Nat.le_pow_clog (by decide) _)
```

### Processor/DispatchSize.lean

Lines 672–699

```lean
theorem dispatch
    (hg : ∀ k, Weight F (fun x : A.CountBit r (B * P.count) → Bool => if A.guard hn k x then 1 else 0) g)
    (hf : ∀ k, ∃ f, Weight F f m ∧ C.Handles k f) :
    ∃ U v, C.Agrees U ∧ Net F (Machine.tape N l l) .unit U v ∧
      v ≤ dispatchCost n * (g + m + (l + l + N / 4) + 1) := by
  choose f hm hsite using hf
  let I := Machine.Ports N l l
  let addr : Header n ⊕ I → I := Sum.elim (fun i => .inr (.inl (C.head i))) id
  have hs := Weight.switch (fun s => C.step f (parse s)) (fun s => C.step_net f hg hm (parse s))
  have hrun : Weight F (C.runAt f) (dispatchCost n * (g + m + (l + l + N / 4) + 1)) := by
    apply (hs.pull addr).mono
    have hi : Fintype.card I = l + l + N / 4 := by simp [I, Machine.Ports, Nat.add_assoc]
    have hh : 1 + Nat.clog 2 n ≤ l := by have := C.fits; omega
    simp only [Header, Fintype.card_sum, Fintype.card_unit, Fintype.card_fin]
    unfold dispatchCost
    ring_nf
    omega
  obtain ⟨v, hv, net⟩ := hrun (Machine.tape N l l) (Machine.address N l l).symm
  refine ⟨_, v, ?_, net, hv⟩
  have he (ξ z : Fin l → Bool) (y : Fin (N / 4) → Bool) :
      bits (Machine.tape N l l) (Machine.feed ξ z y) ∘ (Machine.address N l l).symm =
        Sum.elim ξ (Sum.elim z y) := by
    funext j
    rw [Function.comp_apply, Machine.bits_feed, Equiv.apply_symm_apply]
  constructor
  · simp only [pad, he, runAt_zero, ite_true]
  · intro i c z y
    simpa only [he] using C.runAt_eq f hsite i c z y
```

Lines 720–736

```lean
def ControlBound (F : Type) [Field F] [Fintype F] : Prop :=
  ∀ n h C Q : ℕ, ∀ hn : 2 ≤ n, 0 < h → 4 ≤ C → ∀ A : Expansion F n h,
    ∃ k : ℕ, ∀ᶠ w : ℕ in atTop,
      let N := 2 ^ w
      let r := 2 ^ (w / (16 * n))
      let d := 2 ^ (w / (32 * n))
      let l := N / 2048
      let B := (3 * r + d - 1) / d
      ∀ fa : Machine.Ports N l l → Label n r,
        (∀ z, multiplicity fa z * r ^ n ≤ 4 * N) →
        (∀ i g, multiplicity (fun t : Fin l => fa (.inl t) i) g * r ≤ 2 * l) →
        let P := Format.prescribed d h (rounds n r h (12 * (C * N * d))) fa
        ∀ M : Bank (Fin n), M.budget F ≤ N * w ^ Q → r ^ n ≤ C * N * d →
          (∀ i, 100 * C * N * d ≤ 2 ^ M.width i * r) →
          Nonempty (Controller P A hn M B) →
          ∃ (S : Controller P A hn M B) (U : Val (tape N l l) → Val .unit → F) (v : ℕ),
            S.Agrees U ∧ Net F (tape N l l) .unit U v ∧ v ≤ N * w ^ k
```

Lines 738–792

```lean
theorem programmable_of_control (hctrl : ControlBound F) (K q : ℕ) :
    Programmable F (Cheap F K q) := by
  intro C hC n h hn hh A
  obtain ⟨k, hk⟩ := hctrl n h C (q + 8) hn hh hC A
  obtain ⟨δ, _, hs⟩ := Diagonal.scales C n h (Fintype.card A.I) hC hn
  refine ⟨k, ?_⟩
  filter_upwards [hk, hs, scales F n (1000000 * C) K q hn, eventually_ge_atTop (64 * n)] with w ctrl sc hb hw
  let N := 2 ^ w
  let r := 2 ^ (w / (16 * n))
  let d := 2 ^ (w / (32 * n))
  let l := N / 2048
  let B := (3 * r + d - 1) / d
  have ⟨hd, hd', _⟩ := batch_floor w n
  have s : Diagonal.Scales C n h (Fintype.card A.I) N r l l d 0 B :=
    sc.2 d 0 hd hd' (by simp)
  have hr : 16 ≤ r := by
    have ht : 4 ≤ w / (16 * n) := (Nat.le_div_iff_mul_le (by omega)).mpr (by omega)
    exact Nat.pow_le_pow_right (n := 2) (by decide) ht
  obtain ⟨M, hm, hbit, hspec⟩ := hb
  have hl : 1 + 3 * r + Fintype.card (M.Bits F) ≤ l := by
    have hsq : r ^ 2 ≤ l := (Nat.pow_le_pow_right (by omega) hn).trans s.labels
    have hlin : 2 * (1 + 3 * r) ≤ l := by nlinarith
    change Fintype.card (M.Bits F) ≤ l / 2 at hbit
    omega
  have hsize : r ^ n ≤ C * N * d :=
    s.labels.trans ((by omega : l ≤ N).trans
      ((Nat.le_mul_of_pos_left N (by omega : 0 < C)).trans (Nat.le_mul_of_pos_right _ hd)))
  have hwidth : ∀ i, 1000000 * C * N * d * rate K q (N + 2) ≤ 2 ^ M.width i * r := by
    intro i
    rw [rate, show 1 + (N + 2) = N + 3 by omega]
    exact (hspec i).1
  have hrepeat (i : Fin n) : B ≤ 2 ^ M.repetitions i :=
    (Nat.mul_le_mul_right_iff hd).mp (s.copies.trans (by rw [mul_comm _ d]; exact (hspec i).2))
  intro fa hf hf'
  let P := Format.prescribed d h (rounds n r h (12 * (C * N * d))) fa
  by_cases hex : Nonempty (Controller P A hn M B)
  · have hroom (i : Fin n) : 100 * C * N * d ≤ 2 ^ M.width i * r := by
      calc
        _ ≤ 1000000 * C * N * d := by gcongr; decide
        _ ≤ 1000000 * C * N * d * rate K q (N + 2) := Nat.le_mul_of_pos_right _ (by unfold rate; omega)
        _ ≤ _ := hwidth i
    obtain ⟨S, U, v, hu, hnet, hv⟩ := ctrl fa hf hf' M hm hsize hroom hex
    refine ⟨U, S.pad, v, hu.1, hnet, hv, ?_⟩
    intro _ D hD G _ had T hcol ht hc hcode
    exact S.programs T (Format.ofScheme_eq T hcol ht) had
      (bill_bound T M hc hD s.pos hsize hwidth) hrepeat hcode U hu
  · obtain ⟨v, hv, hnet⟩ := DNet.discard F (tape N l l)
    let pad := (bits (tape N l l)).symm (fun _ => false)
    refine ⟨fun _ _ => 1, pad, v, rfl, hnet, ?_, ?_⟩
    · rw [wd_word] at hv
      exact (hv.trans s.tape).trans (Nat.le_mul_of_pos_right _ (pow_pos (by omega : 0 < w) _))
    · intro _ D hD G _ _ T hcol ht hc hcode
      have he := controller_exists T M hcode hl (bill_bound T M hc hD s.pos hsize hwidth)
      rw [Format.ofScheme_eq T hcol ht] at he
      exact (hex he).elim
```

### Processor/Resistance.lean

Lines 19–22

```lean
private theorem fit_window {N w W Q s : ℕ} (hs : W ≤ 2 * w) (hb : 2 ^ W ≤ N * w ^ (Q + 10))
    (hw : 200 ≤ w) : 2 ^ W * (100 + W) ^ (Q + s) ≤ N * w ^ ((Q + 10) + 2 * (Q + s)) := by
  rw [pow_add w, pow_mul, ← mul_assoc]
  exact Nat.mul_le_mul hb (Nat.pow_le_pow_left (by nlinarith : 100 + W ≤ w ^ 2) _)
```

Lines 24–58

```lean
theorem control (F : Type) [Field F] [Fintype F] : ControlBound F := by
  classical
  intro n h K Q hn hh hK A
  let t := (Q + 10) + 2 * (Q + 660)
  refine ⟨max t 5 + 1, ?_⟩
  filter_upwards [window n h K Q hn hh hK A, dispatch_scale A hn K t,
    eventually_ge_atTop 200] with w getW dispatch hw
  intro N r d l B fa hf _ P M hm hlabel hroom hex
  obtain ⟨W, H⟩ := getW fa M hm hex
  have hr : 0 < r := lt_of_lt_of_le (by decide) H.three
  let _ : NeZero r := ⟨hr.ne'⟩
  let _ : NeZero (3 * r) := ⟨by omega⟩
  have hn4 : 0 < N / 4 := by have := hex.some.links; have := Nat.le_pow (a := 3 * r) (by decide : 0 < 2); omega
  let _ : NeZero (N / 4) := ⟨hn4.ne'⟩
  let e := Fintype.equivFin (Name n r)
  let _ : LinearOrder (Name n r) := LinearOrder.lift' e e.injective
  have he : P.which = (Format.prescribed d h P.count P.colors).which := rfl
  have hg : ∀ z, Diagonal.multiplicity P.colors z * r ^ n ≤ 4 * N := hf
  have hc (i) : Nonempty (Cells.Domain P ↪ Fin (2 ^ M.width i)) := by
    apply Function.Embedding.nonempty_of_card_le
    rw [Fintype.card_fin]
    apply (Nat.mul_le_mul_right_iff hr).mp
    exact (Cells.domain_size P hlabel (by omega) H.positions).trans (hroom i)
  let a : Cells.Address P M := fun i => (hc i).some
  let C := Cells.install hn he hg a hex.some
  have hwW : w ≤ W := by
    by_contra! hlt
    have hh := Nat.pow_lt_pow_right (by decide : 1 < 2) hlt
    exact hh.not_ge H.limits.bits
  have hM : M.budget F ≤ 2 ^ W * (100 + W) ^ Q :=
    hm.trans (Nat.mul_le_mul H.limits.bits (Nat.pow_le_pow_left (by omega) _))
  refine ⟨C.numberInOrder, dispatch fa M C.numberInOrder ?_⟩
  intro i
  obtain ⟨f, hf, he⟩ := C.handles H he hg a rfl hM i
  exact ⟨f, hf.mono (fit_window H.small H.cost hw), he⟩
```

Lines 62–65

```lean
theorem resists (F : Type) [Field F] [Finite F] (K q : ℕ) :
    Diagonal.Resists F (Cheap F K q) := by
  let _ := Fintype.ofFinite F
  exact resists_of_programmable (programmable_of_control (control F) K q)
```

### Diagonal/Broadcast.lean

Lines 813–852

```lean
theorem resists_of_programmable {L : Set (Table F)} (hU : Programmable F L) : Resists F L := by
  intro C₀ hC₀
  let C := 3 * C₀ + 12
  have hC : 4 ≤ C := by omega
  obtain ⟨n, hn, h, A, hh, ha⟩ := bounded_expansion F (52 * C)
  obtain ⟨k, hk⟩ := hU C hC n h hn hh A
  obtain ⟨δ, hδ, hscale⟩ := scales C n h (Fintype.card A.I) hC hn
  refine ⟨k, δ, hδ, ?_⟩
  filter_upwards [hk, hscale] with w impl scale
  let N := 2 ^ w
  let r := 2 ^ (w / (16 * n))
  let l := N / 2048
  have hN : 0 < N := by unfold N; positivity
  have hr : 0 < r := by unfold r; positivity
  obtain ⟨hδd, hs⟩ := scale
  have s := hs 1 0 (by decide) (by exact hr) (by exact hN.le)
  obtain ⟨fa, hfa, hfl⟩ := exterior s.pos s.labels s.tape
  obtain ⟨U, pad, v, hpad, hnet, hv, hp⟩ := impl fa hfa hfl
  have hwd : wd (tape N l l) ≤ N := by rw [wd_word]; exact s.tape
  refine ⟨tape N l l, U, v, hwd, hnet, hv, ?_⟩
  intro d D m hd hsmall hD hwide
  by_contra! hm
  let d₀ := 2 ^ (w / (32 * n))
  obtain ⟨hd₀, hsmall₀, hcmp⟩ := batch_floor w n
  have hmax : (max D 3) ^ 2 ≤ N := by
    rcases le_total 3 D with h | h
    · rw [max_eq_left h]; exact hD
    · rw [max_eq_right h]
      have hli : 64 * (3 * r) ^ 2 ≤ N := s.links
      have : 0 < (3 * r) ^ 2 := by positivity
      omega
  obtain ⟨v', hm', hfixed⟩ := repeat_trim pad hpad hd (hcmp d (hδd d hsmall)) hwd hwide hm
  obtain ⟨G, hG, had, hcost⟩ := hfixed.diagram
  have sc := hs d₀ (max D 3) hd₀ hsmall₀ hmax
  obtain ⟨S, hcolors, hshape, hcomp⟩ := prepare A hn hh hr hC (fun d' => ha r d' N hN)
    fa hfa hfl sc.literal sc.links sc.copies G (hcost.le.trans hm') (fun v => (had v).1)
    sc.degree sc.padding
  have hb := S.fits_of_compact hh sc hcomp
  have prog := hp (max D 3) hmax G (hcost.le.trans hm') had S hcolors hshape hcomp hb
  exact S.diagonal pad U (S.encodeReport hb) hr sc.cover hpad prog hG
```
