import Std

/-! Local lemmas and a finite counterexample for cubic matching descent.
    No theorem in this file asserts a uniform numerical cutwidth bound. -/

namespace MatchingDescent

abbrev Cut (α : Type u) := α → Prop
def union (s t : Cut α) : Cut α := fun x => s x ∨ t x
def inter (s t : Cut α) : Cut α := fun x => s x ∧ t x
def subset (s t : Cut α) : Prop := ∀ x, s x → t x
def put (s : Cut α) (v : α) : Cut α := fun x => x = v ∨ s x
def Submodular (F : Cut α → Nat) : Prop :=
  ∀ s t, F (union s t) + F (inter s t) ≤ F s + F t

/-- A certified saving survives any intervening additions. -/
theorem marginal_persists {F : Cut α → Nat} (sub : Submodular F)
    {s t : Cut α} {v : α} {d : Nat} (st : subset s t) (fresh : ¬t v)
    (cert : F (put s v) + d ≤ F s) : F (put t v) + d ≤ F t := by
  have U : union (put s v) t = put t v := by
    funext x
    apply propext
    constructor
    · rintro (h | z)
      · exact h.elim Or.inl (fun a => Or.inr (st x a))
      · exact Or.inr z
    · rintro (h | h)
      · exact Or.inl (Or.inl h)
      · exact Or.inr h
  have I : inter (put s v) t = s := by
    funext x
    apply propext
    constructor
    · rintro ⟨h, b⟩
      exact h.elim (fun e => (fresh (e ▸ b)).elim) id
    · intro h
      exact ⟨Or.inr h, st x h⟩
  have step := sub (put s v) t
  rw [U, I] at step
  omega

def putList (s : Cut α) : List α → Cut α
  | [] => s
  | v :: rest => putList (put s v) rest

/-- Apply any list of distinct certified moves and keep every certified unit. -/
theorem batch_descent {F : Cut α → Nat} (sub : Submodular F)
    (s : Cut α) (gain : α → Nat) (moves : List α) (t : Cut α)
    (st : subset s t) (nd : moves.Nodup)
    (fresh : ∀ v ∈ moves, ¬t v)
    (cert : ∀ v ∈ moves, F (put s v) + gain v ≤ F s) :
    F (putList t moves) + (moves.map gain).sum ≤ F t := by
  induction moves generalizing t with
  | nil => simp [putList]
  | cons v rest ih =>
    obtain ⟨outside, nodup⟩ := List.nodup_cons.mp nd
    have step := marginal_persists sub st (fresh v (by simp)) (cert v (by simp))
    have tail := ih (put t v) (fun x hx => Or.inr (st x hx)) nodup
      (by
        intro w hw hw'
        rcases hw' with h | h
        · exact outside (h ▸ hw)
        · exact fresh w (List.mem_cons_of_mem _ hw) h)
      (fun w hw => cert w (List.mem_cons_of_mem _ hw))
    simp only [putList, List.map_cons, List.sum_cons]
    omega

/-- Whole patch proposals may overlap. -/
theorem overlapping_patch {F : Cut α → Nat} (sub : Submodular F)
    (z c : Cut α) (opt : F z ≤ F (inter z c)) : F (union z c) ≤ F c := by
  have h := sub z c
  omega

/-- Projecting away finitely attained auxiliary optimum values preserves
    joint submodularity. This is the step used for the edge-bit cover energy. -/
theorem submodular_projection {G : Cut α → Cut β → Nat} {F : Cut α → Nat}
    (joint : ∀ s t x y, G (union s t) (union x y) + G (inter s t) (inter x y)
                        ≤ G s x + G t y)
    (attains : ∀ s, ∃ x, G s x = F s) (lower : ∀ s x, F s ≤ G s x) :
    Submodular F := by
  intro s t
  obtain ⟨x, hx⟩ := attains s
  obtain ⟨y, hy⟩ := attains t
  have h := joint s t x y
  have hu := lower (union s t) (union x y)
  have hi := lower (inter s t) (inter x y)
  omega

def complement (s : Cut α) : Cut α := fun x => ¬s x
def Hits (patch s : Cut α) : Prop := ∃ x, patch x ∧ s x
noncomputable def hitScore (patch s : Cut α) : Nat :=
  open Classical in if Hits patch s then 1 else 0
noncomputable def mixedScore (patch s : Cut α) : Nat :=
  open Classical in if Hits patch s ∧ Hits patch (complement s) then 1 else 0

theorem hit_submodular (p : Cut α) : Submodular (hitScore p) := by
  classical
  intro s t
  have hu : Hits p (union s t) ↔ Hits p s ∨ Hits p t := by
    unfold Hits union
    grind
  have hi : Hits p (inter s t) → Hits p s ∧ Hits p t := by
    unfold Hits inter
    grind
  unfold hitScore
  by_cases a : Hits p s <;> by_cases b : Hits p t <;>
    by_cases c : Hits p (inter s t) <;> simp_all

theorem complement_submodular {F : Cut α → Nat} (h : Submodular F) :
    Submodular (fun s => F (complement s)) := by
  classical
  intro s t
  have u : complement (union s t) = inter (complement s) (complement t) := by
    funext x; apply propext; simp [complement, union, inter]
  have i : complement (inter s t) = union (complement s) (complement t) := by
    funext x; apply propext; simp [complement, union, inter]; grind
  have both := h (complement s) (complement t)
  dsimp only
  rw [u, i]
  omega

theorem mixed_submodular (p : Cut α) : Submodular (mixedScore p) := by
  classical
  by_cases present : ∃ x, p x
  · have charge (s : Cut α) :
        mixedScore p s + 1 = hitScore p s + hitScore p (complement s) := by
      obtain ⟨v, hp⟩ := present
      have either : Hits p s ∨ Hits p (complement s) := by
        by_cases hv : s v
        · exact Or.inl ⟨v, hp, hv⟩
        · exact Or.inr ⟨v, hp, hv⟩
      by_cases a : Hits p s <;> by_cases b : Hits p (complement s) <;>
        simp_all [mixedScore, hitScore]
    intro s t
    have plus := hit_submodular p s t
    have minus := complement_submodular (hit_submodular p) s t
    dsimp only at minus
    have a := charge s
    have b := charge t
    have c := charge (union s t)
    have d := charge (inter s t)
    omega
  · have none (s : Cut α) : ¬Hits p s := fun ⟨x, h, _⟩ => present ⟨x,h⟩
    intro s t; simp [mixedScore, none]

noncomputable def totalCharge (patches : List (Cut α)) (s : Cut α) : Nat :=
  (patches.map (fun p => mixedScore p s)).sum

theorem charge_submodular (patches : List (Cut α)) : Submodular (totalCharge patches) := by
  induction patches with
  | nil => intro s t; simp [totalCharge]
  | cons p ps ih =>
    intro s t
    have a := mixed_submodular p s t
    have b := ih s t
    change mixedScore p (union s t) + totalCharge ps (union s t) +
           (mixedScore p (inter s t) + totalCharge ps (inter s t)) ≤
           mixedScore p s + totalCharge ps s + (mixedScore p t + totalCharge ps t)
    omega

theorem minimum_exists (f : α → Nat) (initial : α) : ∃ z, ∀ w, f z ≤ f w := by
  classical
  suffices least : ∀ k, (∃ x, f x = k) → ∃ z, ∀ w, f z ≤ f w from
    least (f initial) ⟨initial, rfl⟩
  intro k
  induction k using Nat.strongRecOn with
  | ind k ih =>
    rintro ⟨v, hv⟩
    by_cases better : ∃ w, f w < k
    · obtain ⟨w, hw⟩ := better
      exact ih (f w) hw ⟨w, rfl⟩
    · exact ⟨v, fun w => by have hw : ¬ f w < k := fun h => better ⟨w,h⟩; omega⟩

def joinBits (s : Cut α) (x : Cut β) : Cut (Sum α β)
  | .inl v => s v
  | .inr e => x e

noncomputable def coverScore (patches : List (Cut (Sum α β))) (s : Cut α) : Nat :=
  let cost := fun x : Cut β => totalCharge patches (joinBits s x)
  cost (Classical.choose (minimum_exists cost (fun _ => False)))

/-- For vertex/incident-edge patches this is exactly the minimum over edge bits
    of the number of nonmonochromatic vertex stars. -/
theorem crossing_cover_submodular (patches : List (Cut (Sum α β))) :
    Submodular (coverScore patches) := by
  let G := fun s x => totalCharge patches (joinBits s x)
  apply submodular_projection (G := G)
  · intro s t x y
    have u : union (joinBits s x) (joinBits t y) = joinBits (union s t) (union x y) := by
      funext q; cases q <;> rfl
    have i : inter (joinBits s x) (joinBits t y) = joinBits (inter s t) (inter x y) := by
      funext q; cases q <;> rfl
    have both := charge_submodular patches (joinBits s x) (joinBits t y)
    rw [u, i] at both
    exact both
  · intro s
    exact ⟨Classical.choose (minimum_exists (G s) (fun _ => False)), rfl⟩
  · intro s x
    exact Classical.choose_spec (minimum_exists (G s) (fun _ => False)) x

/-- Edges have two specified ends. Taking different edges with the same ends
    or taking equal ends represents parallel edges or loops. -/
def incident (left right : β → α) (v : α) (e : β) : Prop :=
  left e = v ∨ right e = v

def incidenceStar (left right : β → α) (v : α) : Cut (Sum α β)
  | .inl u => u = v
  | .inr e => incident left right v e

def charged (left right : β → α) (s : Cut α) (x : Cut β) : Cut α := fun v =>
  Hits (incidenceStar left right v) (joinBits s x) ∧
  Hits (incidenceStar left right v) (complement (joinBits s x))

def IsCrossingCover (left right : β → α) (s cover : Cut α) : Prop :=
  ∀ e, (s (left e) ↔ s (right e)) ∨ cover (left e) ∨ cover (right e)

/-- A star is charged precisely when an incident edge bit differs from its vertex bit. -/
theorem charged_iff {left right : β → α} {s : Cut α} {x : Cut β} {v : α} :
    charged left right s x v ↔ ∃ e, incident left right v e ∧ ¬ (x e ↔ s v) := by
  classical
  constructor
  · intro h
    apply Classical.byContradiction
    intro impossible
    have safe (e : β) (inc : incident left right v e) : x e ↔ s v :=
      Classical.byContradiction (fun h => impossible ⟨e, inc, h⟩)
    have same (q : Sum α β) (inStar : incidenceStar left right v q) :
        joinBits s x q ↔ s v := by
      cases q with
      | inl u => change u = v at inStar; subst u; rfl
      | inr e => exact safe e inStar
    obtain ⟨⟨a, ha, yes⟩, ⟨b, hb, no⟩⟩ := h
    exact no ((same b hb).mpr ((same a ha).mp yes))
  · rintro ⟨e, inc, different⟩
    by_cases yes : s v
    · have no : ¬ x e := fun xe => different ⟨fun _ => yes, fun _ => xe⟩
      exact ⟨⟨.inl v, rfl, yes⟩, ⟨.inr e, inc, no⟩⟩
    · have xe : x e := by
        apply Classical.byContradiction
        intro h
        exact different (by simp [h, yes])
      exact ⟨⟨.inr e, inc, xe⟩, ⟨.inl v, rfl, yes⟩⟩

/-- A vertex set covers the crossings if and only if one can choose edge bits
    that charge only vertices in that set. This includes multigraphs. -/
theorem cover_iff_edge_bits (left right : β → α) (s cover : Cut α) :
    IsCrossingCover left right s cover ↔
      ∃ x : Cut β, subset (charged left right s x) cover := by
  classical
  constructor
  · intro hc
    let x : Cut β := fun e =>
      (¬ cover (left e) ∧ s (left e)) ∨ (¬ cover (right e) ∧ s (right e))
    refine ⟨x, ?_⟩
    intro v hv
    obtain ⟨e, inc, different⟩ := charged_iff.mp hv
    apply Classical.byContradiction
    intro uncovered
    have covered := hc e
    apply different
    change ((¬ cover (left e) ∧ s (left e)) ∨ (¬ cover (right e) ∧ s (right e))) ↔ s v
    rcases inc with h | h <;> subst v <;> grind
  · rintro ⟨x, fits⟩ e
    by_cases a : cover (left e)
    · exact Or.inr (Or.inl a)
    by_cases b : cover (right e)
    · exact Or.inr (Or.inr b)
    have agree (v : α) (inc : incident left right v e) (outside : ¬cover v) : x e ↔ s v := by
      apply Classical.byContradiction
      intro h
      exact outside (fits v (charged_iff.mpr ⟨e, inc, h⟩))
    exact Or.inl ((agree (left e) (Or.inl rfl) a).symm.trans (agree (right e) (Or.inr rfl) b))

noncomputable def sizeOn (vs : List α) (cover : Cut α) : Nat :=
  open Classical in (vs.map (fun v => if cover v then 1 else 0)).sum

theorem sizeOn_mono (vs : List α) {c d : Cut α} (inc : subset c d) :
    sizeOn vs c ≤ sizeOn vs d := by
  classical
  induction vs with
  | nil => simp [sizeOn]
  | cons v vs ih =>
    have hv := inc v
    change (if c v then 1 else 0) + sizeOn vs c ≤ (if d v then 1 else 0) + sizeOn vs d
    by_cases a : c v <;> by_cases b : d v <;> simp_all <;> omega

theorem total_is_charged_size (vs : List α) (left right : β → α) (s : Cut α) (x : Cut β) :
    totalCharge (vs.map (incidenceStar left right)) (joinBits s x) =
      sizeOn vs (charged left right s x) := by
  classical
  induction vs with
  | nil => rfl
  | cons v vs ih =>
    change mixedScore _ _ + totalCharge (vs.map (incidenceStar left right)) (joinBits s x) =
      (if charged left right s x v then 1 else 0) + sizeOn vs (charged left right s x)
    rw [ih]
    by_cases h : charged left right s x v <;> unfold charged at h <;>
      simp [mixedScore, charged, h]

/-- The minimized incidence charge is attained by a crossing vertex cover
    and is at most the size of every crossing vertex cover. -/
theorem incidence_score_minimum_cover (vs : List α) (left right : β → α) (s : Cut α) :
    let F := coverScore (vs.map (incidenceStar left right)) s
    (∃ c, IsCrossingCover left right s c ∧ sizeOn vs c = F) ∧
    ∀ c, IsCrossingCover left right s c → F ≤ sizeOn vs c := by
  classical
  let patches := vs.map (incidenceStar left right)
  let cost := fun x : Cut β => totalCharge patches (joinBits s x)
  let best := Classical.choose (minimum_exists cost (fun _ => False))
  have least : ∀ x, cost best ≤ cost x :=
    Classical.choose_spec (minimum_exists cost (fun _ => False))
  change (∃ c, IsCrossingCover left right s c ∧ sizeOn vs c = cost best) ∧
         ∀ c, IsCrossingCover left right s c → cost best ≤ sizeOn vs c
  constructor
  · exact ⟨charged left right s best, (cover_iff_edge_bits ..).mpr ⟨best, fun _ h => h⟩,
           (total_is_charged_size vs left right s best).symm⟩
  · intro c hc
    obtain ⟨x, inc⟩ := (cover_iff_edge_bits ..).mp hc
    apply Nat.le_trans (least x)
    change totalCharge (vs.map (incidenceStar left right)) (joinBits s x) ≤ sizeOn vs c
    rw [total_is_charged_size]
    exact sizeOn_mono vs inc

theorem in_putList (s : Cut α) (vs : List α) (w : α) :
    putList s vs w ↔ s w ∨ w ∈ vs := by
  induction vs generalizing s with
  | nil => simp [putList]
  | cons v vs ih =>
    rw [putList, ih, List.mem_cons]
    simp only [put]
    grind

def oldPrefix {n : Nat} (k : Fin n) : Cut (Fin n) := fun w => w.val ≤ k.val
def lexPrefix {n : Nat} (pivot : Fin n → Fin n) (last : Fin n) : Cut (Fin n) := fun w =>
  (pivot w).val < (pivot last).val ∨ ((pivot w).val = (pivot last).val ∧ w.val ≤ last.val)

/-- Exact strict-prefix statement, including halfway through a tied pivot group.
    Each input vertex is its rank in the old total order. -/
theorem every_advance_prefix {n : Nat} (F : Cut (Fin n) → Nat) (sub : Submodular F)
    (pivot : Fin n → Fin n) (adv : ∀ v, (pivot v).val ≤ v.val)
    (cert : ∀ k v, k.val < v.val → (pivot v).val ≤ k.val →
              F (put (oldPrefix k) v) ≤ F (oldPrefix k)) (last : Fin n) :
    F (lexPrefix pivot last) ≤ F (oldPrefix (pivot last)) := by
  classical
  let old := oldPrefix (pivot last)
  let new := lexPrefix pivot last
  have covers : subset old new := by
    intro w hw
    change w.val ≤ (pivot last).val at hw
    change (pivot w).val < (pivot last).val ∨ ((pivot w).val = (pivot last).val ∧ w.val ≤ last.val)
    have hw' := adv w
    have hl := adv last
    omega
  let extras := (List.finRange n).filter (fun v => decide (new v ∧ ¬old v))
  have membership (v : Fin n) : v ∈ extras ↔ new v ∧ ¬old v := by
    simp [extras, List.mem_finRange]
  have nd : extras.Nodup := List.Pairwise.filter _ (List.nodup_finRange n)
  have moves := batch_descent sub old (fun _ => 0) extras old (fun _ x => x) nd
    (fun v hv => ((membership v).mp hv).2)
    (by
      intro v hv
      obtain ⟨hn, ho⟩ := (membership v).mp hv
      change ¬ v.val ≤ (pivot last).val at ho
      change (pivot v).val < (pivot last).val ∨
             ((pivot v).val = (pivot last).val ∧ v.val ≤ last.val) at hn
      exact cert (pivot last) v (by omega) (by omega))
  have same : putList old extras = new := by
    funext w
    apply propext
    rw [in_putList, membership]
    have inc := covers w
    grind
  rw [same] at moves
  change F new ≤ F old
  omega

/-- The all-prefix result specialized to the incidence-cover objective.
    Its submodularity has been discharged in this file. -/
theorem incidence_cover_every_prefix {n : Nat} (patches : List (Cut (Sum (Fin n) β)))
    (pivot : Fin n → Fin n) (adv : ∀ v, (pivot v).val ≤ v.val)
    (cert : ∀ k v, k.val < v.val → (pivot v).val ≤ k.val →
      coverScore patches (put (oldPrefix k) v) ≤ coverScore patches (oldPrefix k))
    (last : Fin n) :
    coverScore patches (lexPrefix pivot last) ≤ coverScore patches (oldPrefix (pivot last)) :=
  every_advance_prefix _ (crossing_cover_submodular patches) pivot adv cert last

structure Availability where
  early : Bool
  late : Bool
deriving DecidableEq, Repr

def sound (b x : Availability) : Prop :=
  (b.early = true → x.early = true) ∧ (x.late = true → b.late = true)

def children (inside : Bool) (a b : Availability) : Availability :=
  if inside then ⟨!(a.late || b.late), false⟩ else ⟨false, !(a.early || b.early)⟩

inductive ExposureTree where
  | frontier (inside exposed : Bool)
  | branch (inside : Bool) (left right : ExposureTree)

def exact : ExposureTree → Availability
  | .frontier i f => ⟨i && f, !i && f⟩
  | .branch i l r => children i (exact l) (exact r)

def bound : ExposureTree → Availability
  | .frontier i _ => ⟨false, !i⟩
  | .branch i l r => children i (bound l) (bound r)

theorem children_sound (i : Bool) (a b x y : Availability)
    (ax : sound a x) (by' : sound b y) : sound (children i a b) (children i x y) := by
  cases i <;> simp_all [sound, children] <;> grind

/-- Sound for all boundary exposures, even when those boundary values are dependent. -/
theorem exposure_sound (t : ExposureTree) : sound (bound t) (exact t) := by
  induction t with
  | frontier i f => cases i <;> cases f <;> simp [sound, bound, exact]
  | branch i l r hl hr => exact children_sound i _ _ _ _ hl hr

/- A bipartite matching recursively skips a left vertex or consumes one
   adjacent, available right vertex. Distinct listed vertices give the
   ordinary cardinality problem for a cut. -/
inductive Feasible (adj : Nat → Nat → Bool) : List Nat → List Nat → Nat → Prop where
  | empty (rs) : Feasible adj [] rs 0
  | skip : Feasible adj ls rs k → Feasible adj (v :: ls) rs k
  | take : w ∈ rs → adj v w = true → Feasible adj ls (rs.erase w) k →
           Feasible adj (v :: ls) rs (k+1)

def sizes (adj : Nat → Nat → Bool) : List Nat → List Nat → List Nat
  | [], _ => [0]
  | v :: ls, rs =>
    sizes adj ls rs ++ (rs.filter (adj v)).flatMap
      (fun w => (sizes adj ls (rs.erase w)).map Nat.succ)

theorem sizes_correct (adj : Nat → Nat → Bool) (ls rs : List Nat) (k : Nat) :
    k ∈ sizes adj ls rs ↔ Feasible adj ls rs k := by
  induction ls generalizing rs k with
  | nil =>
    constructor
    · intro h
      have : k = 0 := by simpa [sizes] using h
      subst k
      exact .empty rs
    · intro h
      cases h
      simp [sizes]
  | cons v ls ih =>
    simp only [sizes, List.mem_append, List.mem_flatMap, List.mem_filter, List.mem_map]
    constructor
    · rintro (s | ⟨w, ⟨member, edge⟩, a, valid, rfl⟩)
      · exact .skip ((ih rs k).mp s)
      · exact .take member edge ((ih _ a).mp valid)
    · intro h
      cases h with
      | skip old => exact Or.inl ((ih _ _).mpr old)
      | take member edge old =>
        exact Or.inr ⟨_, ⟨member, edge⟩, _, (ih _ _).mpr old, rfl⟩

theorem zero_possible (adj : Nat → Nat → Bool) (ls rs : List Nat) :
    0 ∈ sizes adj ls rs := by
  induction ls with
  | nil => simp [sizes]
  | cons v ls ih => exact List.mem_append_left _ ih

def matchingNumber (adj : Nat → Nat → Bool) (ls rs : List Nat) : Nat :=
  (sizes adj ls rs).max (List.ne_nil_of_mem (zero_possible adj ls rs))

theorem maximum_spec (adj : Nat → Nat → Bool) (ls rs : List Nat) :
    Feasible adj ls rs (matchingNumber adj ls rs) ∧
    ∀ k, Feasible adj ls rs k → k ≤ matchingNumber adj ls rs := by
  constructor
  · exact (sizes_correct ..).mp (List.max_mem _)
  · intro k hk
    exact List.le_max_of_mem ((sizes_correct ..).mpr hk)

/-- Giving a matching more available right vertices keeps it feasible. -/
theorem enlarge_right {adj : Nat → Nat → Bool} {ls r : List Nat} {k : Nat}
    (witness : Feasible adj ls r k) {s : List Nat} (grow : r.Sublist s) :
    Feasible adj ls s k := by
  induction witness generalizing s with
  | empty => exact .empty _
  | skip _ ih => exact .skip (ih grow)
  | take member edge _ ih =>
    exact .take (grow.subset member) edge (ih (List.Sublist.erase _ grow))

theorem new_left_upper {adj : Nat → Nat → Bool} {v : Nat} {ls rs : List Nat} {k : Nat}
    (witness : Feasible adj (v :: ls) rs k) : k ≤ matchingNumber adj ls rs + 1 := by
  cases witness with
  | skip old =>
    have h := (maximum_spec adj ls rs).2 _ old
    omega
  | take member edge old =>
    have h := (maximum_spec adj ls rs).2 _ (enlarge_right old List.erase_sublist)
    omega

/-- Adding one left vertex improves the optimum by one exactly when some
    neighbor can be left unused by an old maximum matching. -/
theorem new_vertex_exposure (adj : Nat → Nat → Bool) (v : Nat) (ls rs : List Nat) :
    matchingNumber adj (v :: ls) rs = matchingNumber adj ls rs + 1 ↔
    ∃ w ∈ rs, adj v w = true ∧ Feasible adj ls (rs.erase w) (matchingNumber adj ls rs) := by
  constructor
  · intro improved
    have witness := (maximum_spec adj (v :: ls) rs).1
    rw [improved] at witness
    cases witness with
    | skip old =>
      have h := (maximum_spec adj ls rs).2 _ old
      omega
    | take member edge old => exact ⟨_, member, edge, old⟩
  · rintro ⟨w, member, edge, old⟩
    have lo := (maximum_spec adj (v :: ls) rs).2 _ (.take member edge old)
    have hi := new_left_upper (maximum_spec adj (v :: ls) rs).1
    omega

def cube : Nat → List Nat
  | 0 => [1,3,4] | 1 => [0,2,5] | 2 => [1,3,6] | 3 => [0,2,7]
  | 4 => [0,5,7] | 5 => [1,4,6] | 6 => [2,5,7] | 7 => [3,4,6]
  | _ => []

/-- Nonbacktracking messages for one cut on the simple cube. -/
def cubeMessage (b : Nat) : Nat → Nat → Nat → Availability
  | 0, _, v => ⟨false, !(b.testBit v)⟩
  | d+1, p, v =>
    let msgs := ((cube v).filter (· != p)).map (fun w => cubeMessage b d v w)
    if b.testBit v then ⟨!(msgs.any Availability.late), false⟩
    else ⟨false, !(msgs.any Availability.early)⟩

def cubeAdvance (b radius : Nat) : Nat :=
  (List.range 8).foldl (fun bits v =>
    let msgs := (cube v).map (fun w => cubeMessage b (radius-1) v w)
    if b.testBit v || msgs.any Availability.early || !(msgs.any Availability.late)
    then bits ||| (1 <<< v) else bits) 0

def cubeTwoPasses (b : Nat) : Nat :=
  255 ^^^ cubeAdvance (255 ^^^ cubeAdvance b 3) 3

def cubeRetime (b : Nat) : Nat := cubeAdvance (cubeTwoPasses b) 4

def cubeWidth (b : Nat) : Nat :=
  matchingNumber (fun u v => (cube u).contains v)
    ((List.range 8).filter b.testBit) ((List.range 8).filter (fun v => !(b.testBit v)))

def H (b : Nat) : Nat := cubeWidth (cubeRetime b)

set_option maxRecDepth 16384
set_option maxHeartbeats 4000000

theorem cube_four_images : [1, 3, 5, 7].map cubeRetime = [0, 0, 0, 39] := by decide
theorem cube_four_values : [1, 3, 5, 7].map H = [0, 0, 0, 3] := by decide

def MaskSubmodular (f : Nat → Nat) : Prop :=
  ∀ a b, f (a &&& b) + f (a ||| b) ≤ f a + f b

theorem retimed_cube_not_submodular : ¬ MaskSubmodular H := by
  intro h
  have bad := h 3 5
  have strict : H 3 + H 5 < H (3 &&& 5) + H (3 ||| 5) := by decide
  omega

theorem cube_short_pass_values :
    [7, 15, 23, 31].map (fun b => cubeWidth (cubeAdvance b 3)) = [3, 0, 2, 0] ∧
    [1, 3, 5, 7].map (fun b => cubeWidth (cubeTwoPasses b)) = [0, 0, 0, 3] := by decide

/-- Ending the rule after its first or its second pass also breaks submodularity. -/
theorem shorter_passes_not_submodular :
    ¬ MaskSubmodular (fun b => cubeWidth (cubeAdvance b 3)) ∧
    ¬ MaskSubmodular (fun b => cubeWidth (cubeTwoPasses b)) := by
  constructor
  · intro h
    exact (by decide : ¬ (cubeWidth (cubeAdvance (15 &&& 23) 3) +
      cubeWidth (cubeAdvance (15 ||| 23) 3) ≤
      cubeWidth (cubeAdvance 15 3) + cubeWidth (cubeAdvance 23 3))) (h 15 23)
  · intro h
    exact (by decide : ¬ (cubeWidth (cubeTwoPasses (3 &&& 5)) +
      cubeWidth (cubeTwoPasses (3 ||| 5)) ≤
      cubeWidth (cubeTwoPasses 3) + cubeWidth (cubeTwoPasses 5))) (h 3 5)

def blockMask (b : Nat) := (b &&& 7) ||| (if b.testBit 3 then 248 else 0)
theorem cube_gaussian_sign_table :
    (List.range 16).map (fun b => H (blockMask b)) =
    [0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0] := by decide

def independentAxis (v : Fin 8) : Nat := if v.val < 3 then v.val else 3
def correlatedAxis (v : Fin 8) : Nat := if v.val = 2 then 1 else independentAxis v
def covarianceEntry (assign : Fin 8 → Nat) (i j : Fin 8) : Nat :=
  if assign i = assign j then 1 else 0

/-- Each coordinate is one of four independent standard normal axes.
    Identifying axes 1 and 2 only increases covariance entries. -/
theorem cube_covariance_order :
    ∀ i j : Fin 8, covarianceEntry independentAxis i j ≤ covarianceEntry correlatedAxis i j := by
  decide

def correlatedBlockMask (b : Nat) : Nat :=
  blockMask ((b &&& 11) ||| (if b.testBit 1 then 4 else 0))

/-- All sixteen sign assignments of the four independent normals have equal
    mass at zero. These are the two expectation numerators with denominator 16. -/
theorem cube_expectation_numerators :
    ((List.range 16).map (fun b => H (blockMask b))).sum = 3 ∧
    ((List.range 16).map (fun b => H (correlatedBlockMask b))).sum = 6 := by decide

/-- Cross multiplication for 3/16 < 3*(3/4)^3*(1-3/4). -/
theorem noncentral_probability_is_larger :
    3 * 4^4 < 16 * (3 * 3^3 * (4 - 3)) := by decide

theorem coefficient_test (m n : Nat) : 49*(2*m) < 10*(2*m+n) ↔ 39*m < 5*n := by omega

def bipSign (side : Bool) : Int := if side then 1 else -1

/-- The local (3I+A) filter kills the bipartition sign vector on any cubic
    graph, pointwise. All three neighbor bits are opposite the center. -/
theorem kills_bipartition (v a b c : Bool) (ha : a = !v) (hb : b = !v) (hc : c = !v) :
    3*bipSign v + bipSign a + bipSign b + bipSign c = 0 := by
  subst a; subst b; subst c
  cases v <;> decide

end MatchingDescent
