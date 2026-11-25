import OAI.LinearAlgebra.MatrixExponent.BatchDescent.Hard
import OAI.LinearAlgebra.MatrixExponent.Diagram.Present
import Mathlib.Data.ZMod.Basic

namespace CircuitBreakerList
open Filter OAI.LinearAlgebra.MatrixExponent

/-- The external input to the gate-to-network argument, specialized to one list and F₂. -/
theorem singleList :
    ∃ σ : ℝ, 0 < σ ∧ ∀ᶠ w : ℕ in atTop, ∀ v : ℕ,
      Net (ZMod 2) (lists w 1) .unit (fun x _ => distinctPower (ZMod 2) w 1 x) v →
      (listBits w : ℝ) * (w : ℝ) ^ σ < (v + listBits w : ℕ) := by
  obtain ⟨σ, hσ, h⟩ := BatchDescent.manyList_full (ZMod 2)
  refine ⟨σ, hσ, ?_⟩
  filter_upwards [h] with w hw
  intro v hv
  have hsize : Realizes (ZMod 2) w 1 (v + listBits w) := ⟨v, hv, by simp⟩
  have ht : (1 : ℝ) ≤ (2 : ℝ) ^ ((w : ℝ) ^ (1 + σ)) :=
    Real.one_le_rpow (by norm_num) (Real.rpow_nonneg (Nat.cast_nonneg _) _)
  simpa using hw 1 (v + listBits w) (by decide) (by simpa using ht) hsize

/-- The same bound on the explicit graph used for the gate and copy construction. -/
theorem singleListDiagram :
    ∃ σ : ℝ, 0 < σ ∧ ∀ᶠ w : ℕ in atTop,
      ∀ G : Diagram (ZMod 2) (lists w 1) .unit, G.Trivalent →
      G.eval = (fun x _ => distinctPower (ZMod 2) w 1 x) →
      (listBits w : ℝ) * (w : ℝ) ^ σ < (Fintype.card G.V + listBits w : ℕ) := by
  obtain ⟨σ, hσ, h⟩ := singleList
  refine ⟨σ, hσ, ?_⟩
  filter_upwards [h] with w hw
  intro G hG he
  exact hw _ (he ▸ G.net hG)

#print axioms singleList
#print axioms singleListDiagram
end CircuitBreakerList
