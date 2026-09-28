import Std.Tactic

/-! Integer resource consequences. The operational/entropy hypotheses are
premises, not conclusions of this file. No quantum state semantics are formalized.
-/
namespace CausalResources

def slots (n r : Nat) : Nat := n * r / 2
def exports (r : Nat) : Nat → Nat
  | 0 => 0
  | n + 1 => exports r n + (slots (n + 1) r - slots n r)

theorem slots_monotone (n r : Nat) : slots n r ≤ slots (n + 1) r := by
  unfold slots
  rw [Nat.add_mul, Nat.one_mul]
  omega

theorem export_schedule (n r : Nat) : exports r n = slots n r := by
  induction n with
  | zero => simp [exports, slots]
  | succ n ih =>
      simp only [exports, ih]
      have h := slots_monotone n r
      omega

theorem budget_inversion (purity Q excess : Nat)
    (h : purity ≤ 14 * Q + 4 * excess + 3) :
    purity - 4 * excess - 3 ≤ 14 * Q := by omega

theorem gadget_dimension : 1 + 4 * 3 + 3 * (2 + 12 + 12 + 5 + 6) - 4 * 5 = 104 := by decide
theorem matched_rank : 26 + 6 = 32 := by decide
theorem matched_system_dimension : 2 * 104 = 208 := by decide
theorem dyadic_rank : 32 = 2 ^ 5 := by decide
theorem dephasing_colors : 3 * 174 * 173 + 1 < 2 ^ 17 := by decide
theorem half_rate_after_pairing (n : Nat) : 2 * (n / 2) ≤ n := by omega

#print axioms export_schedule
-- Exact scalar ledger after the analytic weighted-norm estimates in the paper.
-- These do not formalize the matrix inequalities or entropy estimates.
theorem prefix_charge : 50 + 11 * 112 < 1400 := by decide
theorem polar_scale : 16 * 104 < 41 * 41 := by decide
theorem extraction_scale : 4000 * 3 * 41 < 4 * 2^17 := by decide
theorem raw_commutator_charge : 25 + 48 + 24 ≤ 200 := by decide
theorem product_commutator_charge : 4 * (200 + 24) = 896 := by decide
theorem rounded_commutator_charge :
    896 + 4 * 464 + 4 * 12 < 20000 ∧
    448 + 2 * (16 + 464) + 2 * (6 + 12) < 20000 := by decide
theorem displacement_margin : 20000 * 21^2 < 400 * (30000 - 1) := by decide
theorem sector_margin : 16 * (499^2 - 3 - 3000) > 3000000 := by decide
theorem entropy_scalar_margin : 1500000 < 2^21 ∧ 24 * 2 * 1000 < 3 * 1000000 := by decide
theorem centrality_scalar_margin : 19 < 20 ∧ 40 < 3 * 4^2 := by decide

#print axioms rounded_commutator_charge
#print axioms sector_margin
#print axioms budget_inversion
end CausalResources
