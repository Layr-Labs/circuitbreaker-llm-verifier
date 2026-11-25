# From cheap tagged tests to the two implementations at a common width

**Where the cheap service comes from.** Negating eventual tagged hardness yields, for each fixed ε>0, arbitrarily large w where a single network for the exact square tagged permutation costs at most 2^w·w·w^ε. This is the negation with the quantifiers in order; the same network is reusable at the chosen width. The high-level compiler uses a sorting primitive for lists of 2^j records, w≤j≤2w, with w to 4w bits per record and a comparison having the specified Fast implementation. Its only other nontrivalent tables are in Cheap(F,50000,2), with arity bounded both by 2^((ceil log₂w)^3) and by 2^floor(w/2).

**Why there is a strict gain per group of subrounds.** Consider six sites A,B,C,D,M,N. The two demands are AB and CD; the eight permitted supplies are AC, AM, CM, DN, BC, BN, CN and DM. For a spectral distance put u=(AB+CD)/4. If every supply length were at most u, paths ACB and CMD would force AB=CD=2u. Then M and N would be midpoints of CD. Uniqueness of midpoints identifies them up to distance zero; hence MB≤u. Both C and M would now be midpoints of AB, giving CM=0 and CD≤u. This forces u=0, which conflicts with the connected supply carrying a nontrivial bit. For each price there is thus a strict comparison with a supply. The powered-gadget argument converts these strict comparisons into one finite set of shapes of integer total size strictly less than 4b per block b, uniformly for all amplification factors. The finite maximum gives a constant contraction ratio ρ<1 in the subround recipe. When tagged hardness fails, choose_rotor picks a reset graph and a fixed t>0 so its reset overhead is absorbed by ρ^t, even if exceptional tag lengths are sparse. Its construction chooses p with 2ρ^p<1 and then an exceptional size sufficiently large to cover a whole block of t subrounds. The exact inequality is 250(1+c+ceil(log₂w₀))ρ^t≤1/2 for the reset parameters c,w₀. One such block with its reset returns to degree eight. Iteration of this fixed block J times gives the factor 2^(−J).

**Where the saving comes from.** The encoded transfer runs J rounds on replicated diagrams. Its exact interior/exterior allowance is (C+1)(2^(H−J)kM+2^H kE), for original size parameter M, exterior width E, block height H≥J and batch multiplicity d=2^H k. Division by d gives (C+1)(2^(−J)M+E). This is the source of the power saving; shrinking the exterior term geometrically as well would be incorrect. The finite-field recipe obtains the J-round transfer using the subround and reset rules. The first assembles preparation, marking, issuance, candidate selection and sweep with explicit bit-width and cost bounds. The reset chooses a bounded slice and layout under the tag-cost and height assumptions. The completed sequence starts from the true input diagram, controls both the load and the number of exposed endpoints, and preserves the product of the original tensors. The original diagram may have arbitrary incidence; it enters as the argument G of encoded_of_rules. The zero-vertex case is handled separately by one copy.

**Turning transfer into tables.** Host rows, bin rows and exterior rows are counted together. If T denotes the allowance, the setup has enough bins because floor(512T/(v w))+1 bins give between 512T and 513T total capacity in units v w when v w≤T. The row estimate is rows·w≤16500T; padding up to a power of two and to at least 2^w costs at most 100000T in width-weighted row work. The transfer terminates in controls and an exterior screen which evaluate the original power by the preserved endpoint data. Batch replacement of label tables yields a wide sorting implementation of at most 10^35 T. All port weights of long sorts enter this accounting. To fund short records choose Q=c+2, R=Q+2c+10, V=R+c+5 and only after the record-cost exponent e is fixed set B=100(1+b+V+e), with b the fixed recipe radix. The quadratic reserve in H keeps the extra addressing space affordable: adding J to R adds J² height bits. The parameter choice gives 2^(BJ)≤w≤2^(2BJ), hence 2^(−J)≤w^(−1/(2B)). Batch multiplicity is at most (2^w)^min(κ,1/20) for each fixed κ>0. All constants multiplying the cost are independent of w, M and the input diagram. The resulting single compiled diagram costs Kd(Mw^(−a)+E), with a=1/(2B)>0.

**Why it gives two outputs in compatible models.** At exceptional widths choose ε=a/16. Short tables have ordinary overhead bounded by a constant times (2+(ceil log₂w)^3)² per incident port, eventually at most a constant times w^ε. They can instead be kept literally as Cheap macros of allowed arity. The full sorts have two services, using the tagged-test network for their ordinary implementation and the short-block wide implementation for the other. The source proof of sort_service below includes the distinct size charges. Substitution into the same graph adds the costs over its vertices, because the same exterior and contraction indices are kept in both cases. With C'=4+10¹²K, take w large enough that C'+1≤w^ε. Adding dE exterior bits on the ordinary side charges at most d(Mw^(−a/2)+E w^(a/8)). On the wide side the cost is at most C'd(Mw^(−a/2)+E), the allowed macros are exactly Cheap(F,50000,2), and D=2^floor(w/2) satisfies D²≤2^w. These are the two descent inequalities with α=a/2 used in network-proof.md. The eventually valid compiler and all sufficiently-large side inequalities intersect the arbitrarily large exceptional widths. No single sequence of widths is required to work simultaneously for infinitely many κ; the resistance argument fixes just one.

## Checked source steps

Each excerpt below is a complete named declaration from the pinned source. The surrounding namespaces and imports are identified by the file path. The excerpts are for reading alongside the derivation. The built originals and their imported lemmas supply the elaboration environment. Attribution and the applicable Apache 2.0 license are in verification.md.

### BatchDescent/Dictionary.lean

Lines 378–384

```lean
theorem frequent_tags (bad : ¬ Tagged.Hard F) {ε : ℝ} (hε : 0 < ε) :
    ∃ᶠ w : ℕ in atTop, ∃ m, Net F (Tagged.inputs w (word w)) .unit (Tagged.tag F w (word w)) m ∧
      (m : ℝ) ≤ (w : ℝ) ^ ε * (2 : ℝ) ^ w * w := by
  by_contra! h
  refine bad ⟨ε, hε, ?_⟩
  filter_upwards [h] with w hw m hm
  simpa only [mul_comm, mul_left_comm, mul_assoc] using (hw m hm).le
```

Lines 313–316

```lean
def Service (F : Type) [Field F] (w : ℕ) (C : ℝ) (a b : Wire)
    (k : Val a → Val b → F) (budget : ℕ) : Prop :=
  (∃ v, Net F a b k v ∧ (v : ℝ) ≤ C * budget) ∧
    ∃ s ≤ budget, Wide F (2 ^ (w / 2)) (Processor.Cheap F 50000 2) a b k s
```

Lines 318–321

```lean
theorem arity_cap (w : ℕ) : (2 ^ (w / 2)) ^ 2 ≤ 2 ^ w := by
  rw [← pow_mul]
  apply Nat.pow_le_pow_right (by decide)
  omega
```

Lines 402–404

```lean
def SortTables (F : Type) [Field F] (w : ℕ) : Set (Table F) :=
  {T | ∃ (j : ℕ) (a : Wire) (_ : LinearOrder (Val a)), Fast F a ∧
    w ≤ j ∧ j ≤ 2 * w ∧ w ≤ wd a ∧ wd a ≤ 4 * w ∧ T = sorter F a j}
```

Lines 408–410

```lean
def Allowed (F : Type) [Field F] (w : ℕ) : Set (Table F) :=
  SortTables F w ∪ {T | T ∈ Processor.Cheap F 50000 2 ∧
    Processor.arity T ≤ 2 ^ ((Nat.clog 2 w) ^ 3) ∧ Processor.arity T ≤ 2 ^ (w / 2)}
```

Lines 342–376

```lean
theorem sort_service {a : Wire} [LinearOrder (Val a)] (hf : Fast F a) {w j m : ℕ} {C : ℝ} (hw : 64 ≤ w) (hj : w ≤ j) (hj' : j ≤ 2 * w)
    (hb : 1 ≤ wd a) (hb' : wd a ≤ 4 * w) (hC : 1 ≤ C)
    (tag : Net F (Tagged.inputs w (word w)) .unit (Tagged.tag F w (word w)) m)
    (hc : (m : ℝ) ≤ C * (2 : ℝ) ^ w * w) :
    Service F w C (tree j a) (tree j a) (sorter F a j).2.2 (10 ^ 12 * 2 ^ j * w) := by
  have hp (k : ℕ) (hk : k ≤ j) : 2 ^ (j - k) * 2 ^ k = (2 ^ j : ℕ) := by
    rw [← pow_add, Nat.sub_add_cancel hk]
  constructor
  · have ht := sort_trivalent F a hf w (j - w) m 2 (by omega) hb' tag (by omega)
    rw [Nat.add_sub_of_le hj] at ht
    obtain ⟨v, hv, hv'⟩ := ht
    refine ⟨v, hv', ?_⟩
    have hp' : (2 : ℝ) ^ (j - w) * 2 ^ w = 2 ^ j := by exact_mod_cast hp w hj
    have hs : (w + 10000 * (wd a + 1) + 1 : ℕ) ≤ C * (50000 * (w : ℝ)) := by
      have : w + 10000 * (wd a + 1) + 1 ≤ 50000 * w := by omega
      exact (show (w + 10000 * (wd a + 1) + 1 : ℕ) ≤ 50000 * (w : ℝ) by exact_mod_cast this).trans
        (le_mul_of_one_le_left (by positivity) hC)
    calc
      _ ≤ (4 : ℝ) ^ 2 * (2 ^ (j - w) * (4 * m + 10000 * 2 ^ w * (w + 10000 * (wd a + 1) + 1 : ℕ))) := by exact_mod_cast hv
      _ ≤ 4 ^ 2 * (2 ^ (j - w) * (4 * (C * 2 ^ w * w) + 10000 * 2 ^ w * (C * (50000 * w)))) := by gcongr
      _ = C * (8000000064 * 2 ^ j * w) := by rw [← hp']; ring
      _ ≤ _ := by push_cast; gcongr; norm_num
  · have hle : w / 4 ≤ j := by omega
    have ht := sort_wide F a hf (w / 4) (j - w / 4) (2 ^ (w / 2)) 6 hb (by omega)
      (cheap_arity hw hb') (by omega)
    rw [Nat.add_sub_of_le hle] at ht
    obtain ⟨v, hv, hv'⟩ := ht
    refine ⟨v, hv.trans ?_, hv'⟩
    have hs : 1 + 2 * 2 ^ (w / 4) * wd a ≤ 9 * 2 ^ (w / 4) * w := by
      have : 0 < 2 ^ (w / 4) := by positivity
      nlinarith
    calc
      _ ≤ 4 ^ 6 * (2 ^ (j - w / 4) * (9 * 2 ^ (w / 4) * w)) := by gcongr
      _ = (4 ^ 6 * 9) * 2 ^ j * w := by rw [← hp _ hle]; ring
      _ ≤ _ := by gcongr; norm_num
```

Lines 414–425

```lean
theorem service_of_diagram (G : Diagram F a b) (cost : G.V → ℕ)
    (h : ∀ v, Service F w C (G.table v).1 (G.table v).2.1 (G.table v).2.2 (cost v)) :
    Service F w C a b G.eval (∑ v, cost v) := by
  classical
  constructor
  · choose v hv hcost using fun i => (h i).1
    refine ⟨∑ i, v i, Diagram.realize (Net F) Net.empty Net.par Net.pre Net.post Net.trace G v hv, ?_⟩
    simp only [Nat.cast_sum, mul_sum]
    exact sum_le_sum fun i _ => hcost i
  · choose s hs hw using fun i => (h i).2
    exact ⟨∑ i, s i, sum_le_sum (fun i _ => hs i),
      Diagram.realize _ Wide.empty Wide.par Wide.pre Wide.post Wide.trace G s hw⟩
```

Lines 429–471

```lean
theorem expand {u v m : ℕ} {k : Val a → Val b → F} (hw : 64 ≤ w) (hC : 1 ≤ C)
    (tag : Net F (Tagged.inputs w (word w)) .unit (Tagged.tag F w (word w)) m)
    (hc : (m : ℝ) ≤ C * (2 : ℝ) ^ w * w)
    (hl : ((2 + (Nat.clog 2 w) ^ 3) ^ 2 : ℕ) ≤ C)
    (h : Wide F u (Allowed F w) a b k v) : Service F w C a b k (10 ^ 12 * v) := by
  have table (T : Table F) (ht : Processor.arity T ≤ 3 ∨ T ∈ Allowed F w) :
      Service F w C T.1 T.2.1 T.2.2 (10 ^ 12 * (1 + Processor.arity T)) := by
    rcases ht with ht | ⟨j, b, _, hf, hj, hj', hb, hb', rfl⟩ | ⟨hcheap, hsmall, hD⟩
    · have hk := Net.table (F := F) T.1 T.2.1 T.2.2 ht
      have hu : 3 ≤ 2 ^ (w / 2) := by
        have := Nat.pow_le_pow_right (n := 2) (by decide) (show 2 ≤ w / 2 by omega)
        omega
      obtain ⟨s, hs, h⟩ := hk.includeWide (L := Processor.Cheap F 50000 2) hu
      have h' : Service F w C T.1 T.2.1 T.2.2 4 := ⟨⟨1, hk, by norm_num; linarith⟩, s, hs, h⟩
      exact h'.mono (by omega) (by linarith)
    · refine (sort_service hf hw hj hj' (by omega) hb' hC tag hc).mono ?_ (by linarith)
      rw [arity_sorter, mul_assoc (10 ^ 12)]
      apply Nat.mul_le_mul_left
      have := Nat.mul_le_mul_left (2 ^ j) hb
      rw [mul_assoc 2]
      omega
    · let n := Processor.arity T
      have hlog : Nat.clog 2 (1 + n) ≤ (Nat.clog 2 w) ^ 3 + 1 := by
        apply (Nat.clog_le_iff_le_pow (by decide)).mpr
        rw [pow_succ]
        have := Nat.two_pow_pos ((Nat.clog 2 w) ^ 3)
        omega
      obtain ⟨v, hv, hnet⟩ := show T ∈ Processor.Cheap F 50000 2 from hcheap
      refine ⟨⟨v, hnet, ?_⟩, 1 + n, by omega, Wide.table _ _ _ hD (.inr hcheap)⟩
      have hv' : v ≤ 10 ^ 12 * (1 + n) * ((2 + (Nat.clog 2 w) ^ 3) ^ 2) := calc
        _ ≤ (1 + n) * 50000 * (1 + Nat.clog 2 (1 + n)) ^ 2 := hv
        _ ≤ _ := by
          rw [mul_comm (1 + n) 50000]
          apply Nat.mul_le_mul
          · exact Nat.mul_le_mul_right _ (by norm_num)
          · exact Nat.pow_le_pow_left (by omega) 2
      rw [mul_comm C]
      exact (Nat.cast_le.mpr hv').trans (by rw [Nat.cast_mul]; gcongr)
  obtain ⟨G, he, had, hs⟩ := h.diagram
  have hv := service_of_diagram G (fun i => 10 ^ 12 * (1 + G.arity i)) (fun i => table _ (had i).2)
  rw [he, ← mul_sum] at hv
  rw [← hs]
  exact hv
```

### Spectral/CubeGadget.lean

Lines 17–18

```lean
def sixPair : Fin 8 → Fin 6 × Fin 6 :=
  ![(0, 2), (0, 4), (2, 4), (3, 5), (1, 2), (1, 5), (2, 5), (3, 4)]
```

Lines 48–91

```lean
private theorem six_gap (f : Tensor K (Fin 6) →+* ℝ) (hf : Monotone f) :
    ∃ i, f (sixDemand K) < f (sixSupply K i) ^ 4 := by
  by_contra! bad
  let d := dist f
  let l := (d 0 1 + d 2 3) / 4
  have edge_le (i : Fin 8) : d (sixPair i).1 (sixPair i).2 ≤ l := by
    have h := bad i
    simp only [sixDemand, sixSupply, f.map_mul, price_wire f hf _ _ (n := 2) (by decide),
      Nat.cast_ofNat, ← Real.rpow_add (by norm_num : (0 : ℝ) < 2),
      ← Real.rpow_mul_natCast (by norm_num : (0 : ℝ) ≤ 2)] at h
    have := (Real.rpow_le_rpow_left_iff (by norm_num : (1 : ℝ) < 2)).mp h
    unfold l d
    linarith
  have tri : ∀ a b c, d a c ≤ d a b + d b c := dist_triangle f hf
  have sym : ∀ a b, d a b = d b a := dist_comm f
  have ac : d 0 2 ≤ l := edge_le 0
  have am : d 0 4 ≤ l := edge_le 1
  have cm : d 2 4 ≤ l := edge_le 2
  have nd : d 5 3 ≤ l := (sym 5 3).le.trans (edge_le 3)
  have cb : d 2 1 ≤ l := (sym 2 1).le.trans (edge_le 4)
  have nb : d 5 1 ≤ l := (sym 5 1).le.trans (edge_le 5)
  have cn : d 2 5 ≤ l := edge_le 6
  have md : d 4 3 ≤ l := (sym 4 3).le.trans (edge_le 7)
  have ab : d 0 1 = 2 * l := by
    have : 4 * l = d 0 1 + d 2 3 := by unfold l; ring
    linarith [tri 0 2 1, tri 2 4 3]
  have cd : d 2 3 = 2 * l := by unfold l at ab ⊢; linarith
  have mid {a b c x} (hab : d a b = 2 * l) (h₁ : d a c ≤ l) (h₂ : d c b ≤ l)
      (h₃ : d a x ≤ l) (h₄ : d x b ≤ l) : d c x = 0 := by
    have h (z : Fin 6) (ha : d a z ≤ l) (hb : d z b ≤ l) :
        d a z = d a b / 2 ∧ d z b = d a b / 2 := by
      constructor <;> linarith [tri a z b]
    exact unique_midpoint f hf (h c h₁ h₂).1 (h c h₁ h₂).2 (h x h₃ h₄).1 (h x h₃ h₄).2
  have mn := mid cd cm md cn nd
  have mb : d 4 1 ≤ l := by linarith [tri 4 5 1]
  have h := mid ab ac cb am mb
  have hl : l = 0 := by linarith [tri 2 4 3]
  have hu (i : Fin 8) : f (sixSupply K i) = 1 := by
    have he := le_antisymm ((edge_le i).trans_eq hl) (dist_nonneg f hf _ _)
    unfold d at he
    rw [sixSupply, price_wire f hf _ _ (by decide), he, Real.rpow_zero]
  have := hf (six_connected K)
  simp only [map_ofNat, map_prod, hu, Finset.prod_const_one] at this
  norm_num at this
```

Lines 95–103

```lean
theorem six_site :
    ∃ (b C : ℕ) (H : Finset (Fin 8 → ℕ)), 0 < b ∧ (∀ h ∈ H, ∑ i, h i < 4 * b) ∧
      ∀ z : ℕ, 0 < z → sixDemand K ^ (b * z) ≤
        (C : Tensor K (Fin 6)) * ∑ h ∈ H, ∏ i, sixSupply K i ^ (h i * z) := by
  refine Spectrum.powered_gadget add_le_add_right mul_le_mul_left zero_le one_le exists_bound
    _ _ _ (by decide) (fun _ => one_le_wire K _ _ (by decide)) ?_ (six_connected K)
  intro g hg
  obtain ⟨i, hi⟩ := six_gap K (NNReal.toRealHom.comp g) (fun a b h => hg h)
  exact ⟨i, NNReal.coe_lt_coe.mp hi⟩
```

### Spectrum/Gadget.lean

Lines 148–175

```lean
theorem powered_gadget [Nontrivial S] [Nonempty I]
    (hadd : ∀ a b c : S, a ≤ b → a + c ≤ b + c)
    (hmul : ∀ a b c : S, a ≤ b → c * a ≤ c * b)
    (hzero : ∀ a : S, 0 ≤ a) (hunit : ∀ a : S, a ≠ 0 → 1 ≤ a)
    (hbound : ∀ a : S, ∃ n : ℕ, a ≤ (n : S)) (D : S) (W : I → S) (d : ℕ) (hd : 0 < d)
    (hW : ∀ i, (1 : S) ≤ W i)
    (hgap : ∀ g : S →+* ℝ≥0, Monotone g → ∃ i, g D < g (W i) ^ d)
    (hconn : (2 : S) ≤ ∏ i, W i) :
    ∃ (b C : ℕ) (H : Finset (I → ℕ)), 0 < b ∧ (∀ h ∈ H, ∑ i, h i < d * b) ∧
      ∀ z : ℕ, 0 < z → D ^ (b * z) ≤ (C : S) * ∑ h ∈ H, ∏ i, W i ^ (h i * z) := by
  let : AddLeftMono S := ⟨fun c a b h => by simpa only [add_comm] using hadd a b c h⟩
  let : MulLeftMono S := ⟨fun c a b h => hmul a b c h⟩
  obtain ⟨b, p, hb, hp, hprice⟩ :=
    uniform_margin hbound D W d hd hW hgap (Fintype.card I * (Fintype.card I + 1))
  have hy : (∑ i, W i ^ p) ≠ 0 := by
    obtain ⟨i⟩ := ‹Nonempty I›
    have h : (1 : S) ≤ W i ^ p := by simpa using pow_le_pow_left' (hW i) p
    exact ne_of_gt ((lt_of_le_of_ne (hzero 1) zero_ne_one).trans_le
      (h.trans (single_le_sum (f := fun i => W i ^ p) (fun _ _ => hzero _) (mem_univ i))))
  obtain ⟨C, hC⟩ := strictComparison hadd hmul (x := D ^ b) hzero hunit hbound hy (by
    intro g hg
    obtain ⟨i, hi⟩ := hprice g hg
    simp only [map_pow, map_sum]
    exact hi.trans_le (single_le_sum (f := fun i => g (W i) ^ p) (fun _ _ => zero_le) (mem_univ i)))
  obtain ⟨H, hs, he⟩ := profile_bound hzero W hW hconn p
  refine ⟨b, C, H, hb, fun h hh => (hs h hh).trans_lt hp, fun z hz => ?_⟩
  rw [pow_mul]
  exact (hC z).trans (hmul _ _ _ (he z hz))
```

### Spectral/Packets.lean

Lines 32–41

```lean
theorem six_packets :
    ∃ (b C : ℕ) (H : Finset (Fin 8 → ℕ)), 0 < b ∧
      (∀ h ∈ H, ∑ i, h i < 4 * b) ∧ ∀ z : ℕ, 0 < z →
        Packet K (fun e : Fin 2 × Fin (b * z) => sixDemandPair e.1)
          (fun i : Fin C × ↥H => fun e : Σ j : Fin 8, Fin (i.2.1 j * z) => sixPair e.1) := by
  obtain ⟨b, C, H, hb, hcost, h⟩ := six_site K
  refine ⟨b, C, H, hb, hcost, fun z hz => ?_⟩
  apply expand_profile K sixDemandPair (by decide) sixPair
  rw [Fin.prod_univ_two]
  exact h z hz
```

### Packing/FieldLift.lean

Lines 404–412

```lean
structure Gadget (F : Type) [Field F] where
  b : ℕ
  copies : ℕ
  shapes : Finset (Fin 8 → ℕ)
  pos : 0 < b
  cost : ∀ f ∈ shapes, ∑ j, f j < 4 * b
  packet : ∀ z : ℕ, 0 < z →
    Packet F (fun e : Fin 2 × Fin (b * z) => Tensor.sixDemandPair e.1)
      (fun i : Fin copies × ↥shapes => fun e : Σ j : Fin 8, Fin (i.2.1 j * z) => Tensor.sixPair e.1)
```

Lines 414–416

```lean
theorem gadget_exists (F : Type) [Field F] : Nonempty (Gadget F) :=
  let ⟨b, c, H, hb, hc, h⟩ := Tensor.six_packets F
  ⟨⟨b, c, H, hb, hc, h⟩⟩
```

### Packing/Beat.lean

Lines 202–274

```lean
theorem choose_rotor {F : Type} [Field F] (tag : Radix.WideTag F ∨ ¬Tagged.Hard F) (D k : ℕ)
    (ρ B : ℝ) (hρ : 0 ≤ ρ) (hρ' : ρ < 1) (hB : 0 ≤ B) :
    ∃ (t y w : ℕ) (R : Rotor F (2 ^ y) w) (c : ℝ),
      0 < t ∧ 0 < w ∧ 1 ≤ c ∧
      (Fintype.card R.V : ℝ) + Fintype.card R.E ≤ c * (2 : ℝ) ^ y * w ∧
      (∀ Δ ≤ D ^ (k ^ t), 4 * Δ * 2 ^ (64 * Nat.clog 2 (2 * Δ)) ≤ 2 ^ y) ∧
      B * (1 + c + Nat.clog 2 w) * ρ ^ t ≤ 1 / 2 := by
  let b := 131 * (D + 1) * (k + 1)
  have hb : 0 < b := by unfold b; positivity
  obtain ⟨p, hp, hsmall⟩ := ((eventually_gt_atTop 0).and
    ((tendsto_pow_atTop_nhds_zero_of_lt_one hρ hρ').eventually_lt_const
      (by norm_num : (0 : ℝ) < 1 / 4))).exists
  let ζ : ℝ := 1 / ((b : ℝ) * p)
  have hζ : 0 < ζ := by unfold ζ; positivity
  obtain ⟨A, C, htag⟩ := frequent_tag tag hζ
  let K := B * (2 + 6 * |C| + b * p * A)
  have tend := tendsto_pow_atTop_nhds_zero_of_lt_one (show 0 ≤ 2 * ρ ^ p by positivity)
    (show 2 * ρ ^ p < 1 by linarith)
  have lim := tend.const_mul K
  rw [mul_zero] at lim
  obtain ⟨i₀, hi₀⟩ := eventually_atTop.mp ((eventually_gt_atTop 0).and
    (lim.eventually_lt_const (by norm_num : (0 : ℝ) < 1 / 2)))
  let u := 2 ^ (b * p)
  have hu : 1 < u := one_lt_pow₀ (by decide) (Nat.mul_pos hb hp).ne'
  obtain ⟨y, m, hy, hlo, hnet, hs⟩ := htag (u ^ i₀)
  let i := Nat.log u y
  obtain ⟨hi, hdec⟩ := hi₀ i (Nat.le_log_of_pow_le hu hlo)
  let t := p * i
  have ht : 0 < t := Nat.mul_pos hp hi
  have low : 2 ^ (b * t) ≤ y := by
    rw [show b * t = (b * p) * i from (mul_assoc _ _ _).symm, pow_mul]
    exact Nat.pow_log_le_self u hy.ne'
  have up : y ≤ 2 ^ (b * p * (i + 1)) := by
    rw [pow_mul]
    exact (Nat.lt_pow_of_log_lt hu (Nat.lt_succ_self i)).le
  have scale : (y : ℝ) ^ ζ ≤ 2 * (2 : ℝ) ^ i := calc
    _ ≤ ((2 : ℝ) ^ (b * p * (i + 1))) ^ ζ :=
      Real.rpow_le_rpow (by positivity) (by exact_mod_cast up) hζ.le
    _ = (2 : ℝ) ^ (i + 1) := by
      rw [← Real.rpow_natCast_mul (by norm_num), ← Real.rpow_natCast 2 (i + 1)]
      congr 1
      simp only [ζ, Nat.cast_mul]
      field_simp
    _ = _ := pow_succ' _ _
  have clog : Nat.clog 2 (y ^ A) ≤ b * p * A * 2 ^ i := by
    apply Nat.clog_le_of_le_pow
    refine (Nat.pow_le_pow_left up A).trans ?_
    rw [← pow_mul, mul_right_comm (b * p)]
    exact Nat.pow_le_pow_right (by decide) (Nat.mul_le_mul_left _ Nat.lt_two_pow_self)
  obtain ⟨R, hv, he⟩ := Rotor.of_tag hnet
  let c := 1 + 6 * |C| * (2 : ℝ) ^ i
  refine ⟨t, y, y ^ A, R, c, ht, by positivity, le_add_of_nonneg_right (by positivity), ?_, ?_, ?_⟩
  · rw [hv]
    have hem : (Fintype.card R.E : ℝ) ≤ 2 * m := by
      exact_mod_cast (show Fintype.card R.E ≤ 2 * m by omega)
    calc
      _ ≤ 3 * (m : ℝ) := by linarith
      _ ≤ 3 * (|C| * (2 : ℝ) ^ y * (y ^ A : ℕ) * (2 * (2 : ℝ) ^ i)) := by
        apply mul_le_mul_of_nonneg_left (hs.trans ?_) (by norm_num)
        gcongr
        exact le_abs_self _
      _ ≤ _ := by unfold c; nlinarith [show 0 ≤ (2 : ℝ) ^ y * (y ^ A : ℕ) by positivity]
  · intro Δ hΔ
    exact (enough_rows D k t ht Δ hΔ).trans (Nat.pow_le_pow_right (by decide) low)
  · have ha : (Nat.clog 2 (y ^ A) : ℝ) ≤ (b : ℝ) * p * A * (2 : ℝ) ^ i := by exact_mod_cast clog
    have hfact : 1 + c + Nat.clog 2 (y ^ A) ≤ (2 + 6 * |C| + (b : ℝ) * p * A) * (2 : ℝ) ^ i := by
      unfold c
      linarith [one_le_pow₀ (by norm_num : (1 : ℝ) ≤ 2) (n := i)]
    calc
      _ ≤ B * ((2 + 6 * |C| + (b : ℝ) * p * A) * (2 : ℝ) ^ i) * ρ ^ t :=
        mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hfact hB) (by positivity)
      _ = K * (2 * ρ ^ p) ^ i := by simp only [t, K, mul_pow, pow_mul, mul_assoc]
      _ ≤ _ := hdec.le
```

### BatchDescent/Restore/Elect.lean

Lines 548–588

```lean
theorem sequence_recipe (ht : ¬Tagged.Hard F) (hgroup : SubroundRule F) (hreset : ResetRule F) :
    ∃ (P : Recipe F) (K d : ℕ), 0 < K ∧ 3 ≤ d ∧ 4 ^ (P.n * (K + 2) + K + 2) ≤ P.H ∧
      ∀ j : ℕ, 0 < j → Sequence F K 8 8 d d (P.b ^ j) (P.n * j) ((1 / 2) ^ j) P.C := by
  obtain ⟨P⟩ := gadget_exists F
  obtain ⟨t, y, w, R, c, ht, hw, hc, hsize, en, hdecay⟩ :=
    choose_rotor (.inr ht) 8 degreePower P.ratio 250 P.ratio_bounds.1 P.ratio_bounds.2 (by norm_num)
  let D := 8 ^ (degreePower ^ t)
  have : 0 < degreePower := by decide
  have hD : 3 ≤ D := (show 3 ≤ 8 by decide).trans (Nat.le_self_pow (by positivity) _)
  have enough : 4 * coreBins D ≤ 2 ^ y := by simpa only [coreBins, mul_assoc] using en D le_rfl
  have hR : (Fintype.card R.E : ℝ) ≤ c * (2 : ℝ) ^ y * w :=
    (le_add_of_nonneg_left (Nat.cast_nonneg _)).trans hsize
  let r := 2 * resetWidth c w
  have hr : 0 ≤ r := by unfold r resetWidth; positivity
  have bound : r ≤ 250 * (1 + c + Nat.clog 2 w) := by
    unfold r resetWidth
    linarith [Nat.cast_nonneg (α := ℝ) (Nat.clog 2 w)]
  have decrease := (mul_le_mul_of_nonneg_right bound (pow_nonneg P.ratio_bounds.1 t)).trans hdecay
  let b := P.b ^ t * 1
  let n := t + 1
  let e := degreePower ^ t * (900 + 2 ^ y + 1)
  let C := r * t + r
  obtain ⟨X, g⟩ := Reduction.subrounds P 8 (by decide) t ht
  obtain ⟨Y, q⟩ := Reduction.reset R (by positivity) hw hD enough c (by linarith) (by simpa using hR)
  let H₀ := Reduction.fee X (P.ratio ^ t) t * Y
  have step : Reduction F 8 3 b n e H₀ (1 / 2) C :=
    ((g F).comp (q F)).relax le_rfl decrease le_rfl
  obtain ⟨d', reset⟩ := hreset y w D R c hw hD enough hc hR
  let d := max 3 d'
  obtain ⟨v, k, hg⟩ := Sequence.contractions hgroup P 8 d (by decide) t ht
  obtain ⟨k', hk, hq⟩ := reset v
  obtain ⟨l, hl⟩ := Sequence.single hq (by decide) hk hr hr
  let K := max k l + 1
  have seq : Sequence F K 8 8 d d b n (1 / 2) C :=
    ((hg.enlarge (show k ≤ K by unfold K; omega)).comp (hl.enlarge (by unfold K; omega))).relax
      (by decide) (le_max_right _ _) decrease le_rfl
  let H := max (Reduction.fee H₀ (1 / 2) C) (4 ^ (n * (K + 2) + K + 2))
  let recipe : Recipe F := ⟨b, n, e, H, 2 * C, fun j hj =>
    (step.iterate (by decide) (by unfold e; positivity) j hj).relax
      (Nat.pow_le_pow_left (le_max_left _ _) j) le_rfl le_rfl⟩
  exact ⟨recipe, K, d, by unfold K; omega, le_max_left _ _, le_max_right _ _, seq.iterate step.loss⟩
```

Lines 859–865

```lean
theorem subroundRule (F : Type) [Field F] : SubroundRule F := by
  intro P Δ hd d
  apply rule_of_bounds P Δ d hd
  · exact (((prepareWidth_bound P d).add (markWidth_bound P d)).add
      (issueWidth_bound P Δ d)).add ((short_const 2).mul (next_bits P Δ d))
  · have h₁ := (((input_bound P d).scale 100).add (prepareCost_bound P d)).add (headedCost_bound P Δ d)
    exact (((h₁.add (secondCost_bound P Δ d)).add ((candidate_bound P d).scale 4)).add (sweepCost_bound P Δ d))
```

### BatchDescent/Restore/Rule.lean

Lines 666–672

```lean
theorem choose_layout : Nonempty (LocalPrint h s hk P hc cost) := by
  have enough (i : (frame S T u Δ R).state) : (poll h s hk i).Room P (fits h s hk i P hc cost)
      (next P S T u Δ t R (T.scheme.h - r)).scheme := room h s hk i P hc cost
  let A : selection h s hk P hc cost := History.Reset.Setup.select
    (poll h s hk) P (fun i => fits h s hk i P hc cost) (board h s hk)
    (book S T u t R) (meter S t R) (Nat.sub_le T.scheme.h r) enough
  exact ⟨arrange A⟩
```

Lines 679–686

```lean
theorem resetRule (F : Type) [Field F] [Fintype F] : ResetRule F := by
  intro y w Δ P κ hw hΔ hy hc hcost
  let _ : NeZero w := ⟨hw.ne'⟩
  have cost : (Fintype.card P.E : ℝ) ≤ κ * (2 ^ y : ℕ) * w := by simpa only [Nat.cast_pow, Nat.cast_ofNat] using hcost
  obtain ⟨t, s⟩ := exists_slice Δ (2 ^ y) w hy hΔ
  refine ⟨2, fun d => phase_of_local P s hc cost ?_⟩
  intro x y S T R w r u K σ h hk
  exact choose_layout h s hk P hc cost
```

### BatchDescent/Encoding.lean

Lines 19–20

```lean
def allowance (C h J k M E : ℕ) : ℕ :=
  (C + 1) * (2 ^ (h - J) * k * M + 2 ^ h * k * E)
```

Lines 22–28

```lean
theorem allowance_value (C h J k M E : ℕ) (hj : J ≤ h) :
    (allowance C h J k M E : ℝ) =
      (C + 1) * (2 ^ h * k : ℕ) * ((1 / 2 : ℝ) ^ J * M + E) := by
  unfold allowance
  push_cast
  rw [pow_sub₀ (2 : ℝ) (by norm_num) hj, one_div, inv_pow]
  ring
```

Lines 840–853

```lean
private theorem complete {x y : Wire} {S : Setup F x y} {P : Recipe F}
    {K d R B J M w : ℕ} (h : Fit P K d R B J M w S) (hd : 3 ≤ d)
    (q : Sequence F K 8 8 d d (P.b ^ J) (P.n * J) ((1 / 2) ^ J) P.C) :
    Nonempty (Completion S (limits P K d J M w S) (Stack.initialTemplate S w M).sweep
      (Stack.initialBudget S M) S.U S.z 8 8 d 0 (P.n * J) ((1 / 2) ^ J) P.C) := by
  let L := limits P K d J M w S
  have hz : S.z ∣ S.U := by rw [h.size]; exact dvd_mul_left _ _
  have hl : Stack.length F + 1 ≤ L.A := by
    have := h.factors
    have hx : 2 ^ J ≤ L.A := h.lanes.1.trans h.lanes.2
    omega
  have ht := Stack.first_inside S L M d h.enough hd rfl rfl hz hl
  have hσ : 0 ≤ L.σ := by change 0 ≤ (1 / 2 : ℝ) ^ J * M * S.count; positivity
  exact q.run h.enough hσ ht (by rw [zero_add]; exact le_rfl) h.size h.zpos (dvd_refl _) (dvd_refl _) le_rfl
```

Lines 855–891

```lean
theorem encoded_of_rules (ht : ¬Tagged.Hard F) (hgroup : SubroundRule F) (hreset : ResetRule F) :
    Encoded F := by
  obtain ⟨P, K, d, hK, hd, hH, seq⟩ := sequence_recipe ht hgroup hreset
  obtain ⟨C, hC⟩ := exists_nat_gt P.C
  let l := lane P K d
  refine ⟨P, 3 + 2 * C, 3 * l, l, by omega, ?_⟩
  intro Q R V _ hR _
  refine ⟨5 * result P K d R, ?_⟩
  intro B hB
  have hg : label P K d + 10 ≤ B := by unfold result at hB; omega
  filter_upwards [large P K d R B] with J hJ
  intro w z T hQ hV hh hz hv hlo hhi hwidth _ x y G M hn htri hsize hM hM'
  let D := P.dimensions w J z (R + J) Q V (2 ^ (l * J)) T hQ hV hh hz hv
  have hj := hJ.beats
  have hw : 100 ≤ w := by
    have hb : J ≤ B * J := Nat.le_mul_of_pos_left _ (by omega)
    have := (dyad (hj.trans hb)).trans hlo
    omega
  obtain ⟨hm, hp⟩ := Stack.hostCount hw hM hM'
  obtain ⟨e⟩ := Stack.name_exists (x := x) (y := y)
    (by omega : wd x + wd y ≤ M) ((le_max_right _ _).trans hm)
  let S := D.setup G e (cap P K J M D.copies) (P.group J)
  have hu : S.U ≤ w := by
    change P.b ^ J * D.z ≤ w
    rw [Recipe.dimensions_z]
    refine le_trans ?_ hwidth
    rw [← mul_assoc 16, mul_assoc (16 * 2 ^ (Q * J))]
    exact Nat.le_mul_of_pos_left _ (by positivity)
  have h : Fit P K d R B J M w S :=
    ⟨hJ, hK, hH, hg, rfl, rfl, D.zpos, hu, rfl, by change 2 * w ≤ _; omega,
      rfl, rfl, htri, hsize, hm, hp, hlo, hhi⟩
  obtain ⟨t⟩ := complete h hd (seq J (by omega))
  have hr : 2 ^ (l * J) * D.q ≤ 2 ^ (R * J) := by
    change 2 ^ (l * J) * 2 ^ (Q * J) ≤ _
    rw [← pow_add, ← add_mul]
    exact mono J (by omega)
  exact collected h hC.le hn hr t
```

### BatchDescent/Compilation.lean

Lines 136–138

```lean
theorem reserve_height (R J : ℕ) : P.blocks (R + J) J = P.blocks R J + J ^ 2 := by
  unfold blocks
  ring
```

Lines 146–166

```lean
theorem reserve_parameters (R Q V B J : ℕ) (hB : 5 + Q + P.b + V ≤ B)
    {κ : ℝ} (hk : 0 < κ) :
    ∀ᶠ w : ℕ in atTop, ∃ j z : ℕ, J ≤ j ∧ 1 ≤ j ∧ 2 ^ (B * j) ≤ w ∧ w ≤ 2 ^ (2 * B * j) ∧
      0 < z ∧ 2 ^ (V * j) ∣ z ∧ 16 * (2 ^ (Q * j) * P.b ^ j) * z ≤ w ∧ w ≤ 32 * (2 ^ (Q * j) * P.b ^ j) * z ∧
      0 < P.batch (R + j) j z ∧ (P.batch (R + j) j z : ℝ) ≤ ((2 : ℝ) ^ w : ℝ) ^ κ := by
  obtain ⟨N, hN⟩ := eventually_atTop.mp (replication (half_pos hk) 1)
  filter_upwards [P.parameters R Q V B (max J N) hB (half_pos hk)] with w hw
  obtain ⟨j, z, hj, hj', hl, hu, hz, hv, hn, he, hs, hp⟩ := hw
  have hh : 2 ^ j ≤ w := (Nat.pow_le_pow_right (by decide)
    (Nat.le_mul_of_pos_left _ (by omega : 0 < B))).trans hl
  have ht := hN j (by omega) w hh
  have hsize : P.batch (R + j) j z = 2 ^ (j ^ 2) * P.batch R j z := by
    unfold batch
    rw [P.reserve_height, pow_add]
    ring
  refine ⟨j, z, by omega, hj', hl, hu, hz, hv, hn, he, ?_, ?_⟩
  · rw [hsize]; positivity
  · rw [hsize, Nat.cast_mul, Nat.cast_pow, Nat.cast_ofNat]
    rw [one_mul] at ht
    refine (mul_le_mul ht hp (by positivity) (by positivity)).trans_eq ?_
    rw [← Real.rpow_add (by positivity), add_halves]
```

Lines 460–471

```lean
theorem row_count : D.rows S x y * w ≤ 16500 * T := by
  have h₁ := B.host
  have h₂ := B.bin_budget.2
  have h₃ : (wd x + wd y) * D.n * w ≤ 64 * T := by
    calc
      _ ≤ (wd x + wd y) * D.n * (32 * (D.q * D.b) * D.z) := Nat.mul_le_mul_left _ B.wide
      _ = 32 * (wd x + wd y) * (D.n * D.q) * (D.b * D.z) := by ring
      _ = 64 * ((wd x + wd y) * D.copies) := by rw [D.cover]; unfold copies; ring
      _ ≤ _ := Nat.mul_le_mul_left _ B.exterior
  unfold rows
  ring_nf at h₁ h₂ h₃ ⊢
  omega
```

Lines 658–694

```lean
theorem build (budget : D.Budget b.scheme x y w T)
    (hT : T = allowance C D.h J (D.b * D.z) M (wd x + wd y)) (hJ : J ≤ D.h)
    (hu : 2 ^ (c * J) * D.b ≤ D.u) (hv : D.PortCount b.scheme + 2 + 2 * 2 ^ b.scheme.h ≤ D.v)
    (hs : D.shortBudget M e J ≤ T) (hl : D.longBudget M c R J w ≤ T)
    (he : D.e = 2 * w) (ha : D.a = 3 * w) (hh : 8 * D.shortWidth b.scheme (Fintype.card b.Z) ≤ w)
    (hext : D.h + 4 ≤ D.b * D.z)
    (hsmall : 2 * 2 ^ (D.h + 2) * (D.h + 1 + D.b * D.z + 2) ≤ min (2 ^ ((Nat.clog 2 w) ^ 3)) (2 ^ (w / 2)))
    (htop : 10 ^ 38 * T ≤ 2 ^ (5 * w / 4)) :
    ∃ u s, Wide F u (Allowed F w) (copies x D.copies) (copies y D.copies)
      (MatrixExponent.power G.eval D.copies) s ∧ s ≤ 10 ^ 35 * T := by
  let cap := max b.scale (2 * 2 ^ (2 * w) * (4 * w))
  let L := Allowed F w ∪ (b.labels : Set (Table F))
  let j := D.height b.scheme x y w
  obtain ⟨hc, hr⟩ := budget.room he ha hh hext (show 2 * 2 ^ (2 * w) * (4 * w) ≤ cap from le_max_right ..) hsmall
  have done := terminal b.stop b.view hc hr (b.endpoints budget hT hJ hu hv)
  have large : Allowed F w ⊆ L := Set.subset_union_left
  have prog := b.transfer.widen (le_max_left b.scale (2 * 2 ^ (2 * w) * (4 * w))) (Set.Subset.refl L)
  obtain ⟨n, hn, net⟩ := prog.public_end hc.gates (done.widen le_rfl large) b.first b.code b.encoded
  have hb := b.cost
  have hex : wd (copies x D.copies) + wd (copies y D.copies) ≤ T := by
    rw [wd_copies, wd_copies, ← mul_add, mul_comm]
    exact budget.exterior
  have width : b.stop.public.port.width = b.stop.port.width +
      (wd (copies x D.copies) + wd (copies y D.copies)) := by
    rw [show b.stop.public.port = b.stop.port.pair (.wire (.pair _ _)) from rfl, Format.width_pair]
    exact congrArg _ (card_slots (.pair _ _))
  have hcost := budget.charge
  rw [width] at hn
  have h₁ : n ≤ 10 ^ 26 * T := by omega
  have h₃ := budget.min
  have hn' : n ≤ 2 ^ (2 * w - 1) := (h₁.trans ((by omega : 10 ^ 26 * T ≤ 10 ^ 38 * T).trans htop)).trans
    (Nat.pow_le_pow_right (by decide) (by omega))
  rw [b.power] at net
  obtain ⟨u, s, hs, hsize⟩ := batch_labels h₃ b.labels b.few b.small net hn'
  have hw := budget.base
  rw [mul_comm (2 ^ w)] at hw
  exact ⟨u, s, hs, by omega⟩
```

Lines 771–781

```lean
private theorem final_cost (B J w C h k M E s : ℕ) (hB : 0 < B) (hw : 1 ≤ w)
    (hu : w ≤ 2 ^ (2 * B * J)) (hJ : J ≤ h) (hs : s ≤ 10 ^ 35 * allowance C h J k M E) :
    (s : ℝ) ≤ (10 ^ 35 * (C + 1) : ℕ) * (2 ^ h * k : ℕ) *
      ((M : ℝ) * (w : ℝ) ^ (-(1 / (2 * B : ℝ))) + E) := by
  have h := Nat.cast_le (α := ℝ) |>.mpr hs
  rw [Nat.cast_mul, allowance_value _ _ _ _ _ _ hJ] at h
  refine h.trans ?_
  push_cast
  rw [mul_comm (M : ℝ), ← mul_assoc, ← mul_assoc]
  gcongr
  exact saving hB hw hu
```

Lines 783–857

```lean
theorem compilation_of_encoded (h : Encoded F) : Compilation F := by
  obtain ⟨P, C, c, d, hd, h⟩ := h
  let Q := c + 2
  let R := Q + 2 * c + 10
  let V := R + c + 5
  obtain ⟨e, H⟩ := h Q R V le_rfl le_rfl le_rfl
  let B := 100 * (1 + P.b + V + e)
  have hQR : Q ≤ R := by unfold R; omega
  have hVQ : Q ≤ V := by unfold V; omega
  have hB : 5 + Q + P.b + V ≤ B := by unfold B; omega
  have hb : 0 < B := by omega
  obtain ⟨L, hL⟩ := eventually_atTop.mp <| (H B le_rfl).and <|
    (P.reserve_limits R Q B (by omega)).and (P.shallow R)
  refine ⟨10 ^ 35 * (C + 1), 1 / (2 * (B : ℝ)), by positivity, fun κ hk => ?_⟩
  let δ := min κ (1 / 20 : ℝ)
  filter_upwards [P.reserve_parameters R Q V B L hB (show 0 < δ by unfold δ; positivity),
    allowance_small C, eventually_ge_atTop 1] with w hw hall hw'
  obtain ⟨J, z, hL', hJ, hl, hu, hz, hdiv, hlo, hhi, hp, hdim⟩ := hw
  obtain ⟨run, limits, shallow⟩ := hL J hL'
  intro x y k v M net h₀ low high
  obtain ⟨G, rfl, hg, hcard⟩ := net.diagram
  by_cases hG : Nonempty G.V
  swap
  · let _ := not_nonempty_iff.mp hG
    have hv : v = 0 := hcard.symm.trans Fintype.card_of_isEmpty
    obtain ⟨n, hn, hnet⟩ := (net.power 1).includeWide (L := Allowed F w) (le_rfl : 3 ≤ 3)
    have hzero : n = 0 := by simpa [hv] using hn
    refine ⟨1, 3, n, le_rfl, ?_, hnet, ?_⟩
    · rw [Nat.cast_one]
      apply Real.one_le_rpow _ hk.le
      exact_mod_cast Nat.one_le_pow w 2 (by decide)
    · rw [hzero, Nat.cast_zero]
      positivity
  let b := P.blocks (R + J) J
  let T := allowance C b J (P.b ^ J * z) M (wd x + wd y)
  have hq : 0 < Q * J := Nat.mul_pos (by unfold Q; omega) hJ
  have hv : 0 < V * J := Nat.mul_pos (by unfold V; omega) hJ
  have hr : R * J ≤ b := by
    unfold b Recipe.blocks
    exact (Nat.mul_le_mul_right J (Nat.le_add_right ..)).trans (Nat.le_add_right ..)
  have hr' : J ≤ b := (Nat.le_mul_of_pos_left _ (by unfold R; omega)).trans hr
  have hh : Q * J ≤ b + 1 := (Nat.mul_le_mul_right _ hQR).trans (hr.trans (Nat.le_succ b))
  let D := P.dimensions w J z (R + J) Q V (2 ^ (d * J)) T hq hv hh hz hdiv
  have hZ : D.z = z := P.dimensions_z ..
  have copies : D.copies = P.batch (R + J) J z := P.dimensions_copies ..
  have hT : T = allowance C D.h J (D.b * D.z) M (wd x + wd y) := by rw [hZ]; rfl
  obtain ⟨a⟩ := run w z T hq hv hh hz hdiv hl hu hlo hhi x y G M hG hg (by omega) low high
  change Transfer G D M C c e R J w at a
  have hn : (D.copies : ℝ) ≤ (2 ^ w : ℕ) ^ δ := by
    simpa only [copies, Nat.cast_pow, Nat.cast_ofNat] using hdim
  have big : 1 ≤ (↑(2 ^ w : ℕ) : ℝ) := by exact_mod_cast Nat.one_le_pow w 2 (by decide)
  have top : 10 ^ 38 * T ≤ 2 ^ (5 * w / 4) := by
    rw [hT]
    exact hall _ _ _ _ _ hr' (by omega) high
      (hn.trans (Real.rpow_le_rpow_of_exponent_le big (min_le_right ..)))
  have hsmall : 8 * D.shortWidth a.scheme (Fintype.card a.Z) ≤ w :=
    a.short hJ hQR (le_rfl : R + c + 5 ≤ V) (le_rfl : 100 * (1 + P.b + V + e) ≤ B) shallow rfl rfl
      (Nat.sub_le ..) hl
  obtain ⟨budget, short, long, room⟩ := a.funded Q V B hJ le_rfl le_rfl le_rfl
    (by unfold B; omega) hr low hl (by rw [hZ]; exact hlo) (by rw [hZ]; exact hhi) rfl
    (P.dimensions_rows ..) (Nat.pow_le_pow_right (by decide) (Nat.mul_le_mul_right _ hd))
    T hT rfl top
  have hunit : 2 ^ (c * J) * D.b ≤ D.u := by
    apply Nat.mul_le_mul_right
    apply Nat.pow_le_pow_right (by decide)
    unfold Q
    rw [add_mul]
    omega
  obtain ⟨hdepth, hsort⟩ := limits z w hl hu hlo hhi
  obtain ⟨u, s, hwide, hcost⟩ := a.build budget hT hr' hunit room short long rfl rfl hsmall
    (by rw [hZ]; exact hdepth) (by rw [hZ]; exact hsort) top
  refine ⟨D.copies, u, s, by rw [copies]; exact hp,
    hn.trans (Real.rpow_le_rpow_of_exponent_le big (min_le_left ..)), hwide, ?_⟩
  rw [hT] at hcost
  exact final_cost B J w C D.h (D.b * D.z) M (wd x + wd y) s hb hw' hu hr' hcost
```

Lines 867–869

```lean
theorem encoded {F : Type} [Field F] [Finite F] (h : ¬ Tagged.Hard F) : Encoded F := by
  let _ := Fintype.ofFinite F
  exact encoded_of_rules h (subroundRule F) (resetRule F)
```

### BatchDescent/Plan.lean

Lines 59–69

```lean
theorem saving {B j w : ℕ} (hB : 0 < B) (hw : 1 ≤ w) (hu : w ≤ 2 ^ (2 * B * j)) :
    (1 / 2 : ℝ) ^ j ≤ (w : ℝ) ^ (-(1 / (2 * B : ℝ))) := by
  have hb : 0 < (B : ℝ) := Nat.cast_pos.mpr hB
  have hp : (w : ℝ) ^ (1 / (2 * B : ℝ)) ≤ (2 : ℝ) ^ j := calc
    _ ≤ ((2 : ℝ) ^ (2 * B * j : ℕ) : ℝ) ^ (1 / (2 * B : ℝ)) := by gcongr; exact_mod_cast hu
    _ = _ := by
      rw [← Real.rpow_natCast_mul (by norm_num)]
      have : ((2 * B * j : ℕ) : ℝ) * (1 / (2 * B : ℝ)) = j := by push_cast; field_simp
      rw [this, Real.rpow_natCast]
  rw [one_div_pow, Real.rpow_neg (Nat.cast_nonneg _), ← one_div]
  exact one_div_le_one_div_of_le (by positivity) hp
```

### BatchDescent/Encoded.lean

Lines 17–24

```lean
def Compilation (F : Type) [Field F] : Prop :=
  ∃ K : ℕ, ∃ a : ℝ, 0 < a ∧ ∀ κ : ℝ, 0 < κ →
    ∀ᶠ w : ℕ in atTop, ∀ (x y : Wire) (k : Val x → Val y → F) (v M : ℕ),
      Net F x y k v → v + wd x + wd y ≤ M → 2 ^ w ≤ M →
      (M : ℝ) ≤ (2 ^ w : ℕ) ^ (21 / 20 : ℝ) →
      ∃ d u s : ℕ, 1 ≤ d ∧ (d : ℝ) ≤ (2 ^ w : ℕ) ^ κ ∧
        Wide F u (Allowed F w) (copies x d) (copies y d) (power k d) s ∧
        (s : ℝ) ≤ K * d * ((M : ℝ) * (w : ℝ) ^ (-a) + (wd x + wd y : ℕ))
```

Lines 58–78

```lean
private theorem cost {K d M e s v w a : ℝ}
    (hd : 0 ≤ d) (hM : 0 ≤ M) (he : 0 ≤ e) (hw : 1 ≤ w) (ha : 0 < a)
    (hbig : K + 1 ≤ w ^ (a / 16)) (hc : s ≤ K * d * (M * w ^ (-a) + e)) (hv : v ≤ w ^ (a / 16) * s) :
    v + d * e ≤ d * (M * w ^ (-(a / 2)) + e * w ^ (a / 2 / 4)) := by
  have hu : 1 ≤ w ^ (a / 16) := Real.one_le_rpow hw (by linarith)
  have hw' : 0 < w := zero_lt_one.trans_le hw
  let u := w ^ (a / 16)
  let S := d * (M * w ^ (-a) + e)
  have hS : 0 ≤ S := by unfold S; positivity
  have hex : d * e ≤ S := mul_le_mul_of_nonneg_left (le_add_of_nonneg_left (by positivity)) hd
  have hs : s ≤ K * S := by rwa [mul_assoc] at hc
  have hprod : u * u = w ^ (a / 2 / 4) := by unfold u; rw [← Real.rpow_add hw']; congr 1; ring
  have hpow : w ^ (a / 2 / 4) * w ^ (-a) ≤ w ^ (-(a / 2)) := by
    rw [← Real.rpow_add hw']
    exact Real.rpow_le_rpow_of_exponent_le hw (by linarith)
  calc
    _ ≤ u * (K * S) + u * S := add_le_add (hv.trans (by gcongr)) (hex.trans (le_mul_of_one_le_left hS hu))
    _ = u * (K + 1) * S := by ring
    _ ≤ u * u * S := by gcongr
    _ = d * (M * (w ^ (a / 2 / 4) * w ^ (-a)) + e * w ^ (a / 2 / 4)) := by rw [hprod]; unfold S; ring
    _ ≤ _ := by gcongr
```

Lines 82–113

```lean
theorem descent {F : Type} [Field F] (compile : ¬Tagged.Hard F → Compilation F) :
    Diagonal.Descent F (Processor.Cheap F 50000 2) := by
  intro bad
  obtain ⟨K, a, ha, h⟩ := compile bad
  let C := 4 + 10 ^ 12 * K
  refine ⟨C, by unfold C; omega, a / 2, half_pos ha, fun κ hκ => ?_⟩
  have hp : 0 < a / 16 := by linarith
  have large := ((tendsto_rpow_atTop hp).comp tendsto_natCast_atTop_atTop).eventually_ge_atTop ((C : ℝ) + 1)
  refine ((frequent_tags bad hp).and_eventually ((eventually_ge_atTop 64).and
    (large.and ((short_cost hp).and (h κ hκ))))).mono ?_
  rintro w ⟨⟨m, tag, cheap⟩, hw, hlarge, hshort, g⟩
  dsimp only
  intro x y k v M hv hm hn hM
  obtain ⟨d, u, s, hd, hd', network, hs⟩ := g x y k v M hv hm hn hM
  let e := wd x + wd y
  have hw' : 1 ≤ (w : ℝ) := by exact_mod_cast (show 1 ≤ w by omega)
  have hC : 1 ≤ (w : ℝ) ^ (a / 16) := Real.one_le_rpow hw' hp.le
  obtain ⟨⟨r, ht, hr⟩, s', hs', hsnet⟩ := expand hw hC tag cheap hshort network
  have hb : (10 ^ 12 * s : ℕ) ≤ (C : ℝ) * d * ((M : ℝ) * (w : ℝ) ^ (-a) + e) := by
    have hl : (10 : ℝ) ^ 12 * K ≤ C := by unfold C; push_cast; linarith
    rw [Nat.cast_mul, Nat.cast_pow, Nat.cast_ofNat]
    calc
      _ ≤ 10 ^ 12 * ((K : ℝ) * d * ((M : ℝ) * (w : ℝ) ^ (-a) + e)) := by gcongr
      _ = (10 ^ 12 * (K : ℝ)) * d * ((M : ℝ) * (w : ℝ) ^ (-a) + e) := by ring
      _ ≤ _ := by gcongr
  constructor
  · refine ⟨d, r + d * e, hd, hd', ⟨r, ht, rfl⟩, ?_⟩
    rw [Nat.cast_add, Nat.cast_mul]
    exact cost (by positivity) (by positivity) (by positivity) hw' ha hlarge hb hr
  · refine ⟨d, 2 ^ (w / 2), s', hd, hd', arity_cap w, hsnet, ?_⟩
    refine ((Nat.cast_le.mpr hs').trans hb).trans ?_
    gcongr; linarith
```

### BatchDescent/Hard.lean

Lines 17–18

```lean
theorem hard (F : Type) [Field F] [Finite F] : Tagged.Hard F :=
  Processor.hard_of_descent (descent fun bad => compilation_of_encoded (encoded bad))
```

Lines 20–24

```lean
theorem manyList_full (F : Type) [Field F] [Finite F] :
    ∃ η : ℝ, 0 < η ∧ ∀ᶠ w : ℕ in atTop, ∀ t m : ℕ,
      1 ≤ t → (t : ℝ) ≤ (2 : ℝ) ^ ((w : ℝ) ^ (1 + η)) → Realizes F w t m →
      (t : ℝ) * listBits w * (w : ℝ) ^ η < m :=
  Tagged.manyList_full (hard F)
```
