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
#print axioms budget_inversion
end CausalResources
