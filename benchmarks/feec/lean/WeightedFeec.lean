import Mathlib.LinearAlgebra.FiniteDimensional
import Mathlib.LinearAlgebra.Dimension.Constructions
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Tactic

set_option warningAsError true

/-!
# Weighted derivative coefficient structure

For a first-cell polynomial `N(s) = ∑ j=1..p, c_j s^j`, the normalized
weighted derivative is `s^(α+1) (s^(-α) N)' = s N' - α N` on `s > 0`.
The linear map below is its coefficient map. Its kernel/rank statements
concern this finite-dimensional coefficient model. They do not assume or
prove a global FEEC exact sequence, commuting projection, or MHD convergence.
-/

noncomputable section

namespace WeightedFeec

/-- Polynomial degrees run from one through `p`: the first-cell left trace is zero. -/
def degree {p : ℕ} (i : Fin p) : ℝ := (i.val + 1 : ℕ)

theorem degree_injective (p : ℕ) : Function.Injective (@degree p) := by
  intro i j h
  dsimp [degree] at h
  have hn : i.val + 1 = j.val + 1 := by exact_mod_cast h
  apply Fin.ext
  omega

/-- The coefficient map of `s N' - α N`. -/
def normalizedDerivative (p : ℕ) (α : ℝ) : (Fin p → ℝ) →ₗ[ℝ] (Fin p → ℝ) where
  toFun c i := (degree i - α) * c i
  map_add' c d := by
    funext i
    simp only [Pi.add_apply, mul_add]
  map_smul' a c := by
    funext i
    simp only [Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
    ring

@[simp] theorem normalizedDerivative_apply {p : ℕ} (α : ℝ) (c : Fin p → ℝ)
    (i : Fin p) : normalizedDerivative p α c i = (degree i - α) * c i := rfl

/-- Exactly the coefficients at a resonant degree can survive in the kernel. -/
theorem kernel_iff {p : ℕ} (α : ℝ) (c : Fin p → ℝ) :
    normalizedDerivative p α c = 0 ↔ ∀ i, degree i ≠ α → c i = 0 := by
  constructor
  · intro h i hi
    have heq := congrFun h i
    change (degree i - α) * c i = 0 at heq
    exact (mul_eq_zero.mp heq).resolve_left (sub_ne_zero.mpr hi)
  · intro h
    funext i
    change (degree i - α) * c i = 0
    by_cases hi : degree i = α
    · simp [hi]
    · simp [h i hi]

theorem injective_of_nonresonant {p : ℕ} (α : ℝ)
    (hα : ∀ i : Fin p, degree i ≠ α) :
    Function.Injective (normalizedDerivative p α) := by
  intro c d h
  funext i
  have hi := congrFun h i
  exact mul_left_cancel₀ (sub_ne_zero.mpr (hα i)) hi

/-- A resonant monomial is a genuine nonzero kernel vector. -/
theorem resonant_kernel_vector {p : ℕ} (k : Fin p) :
    normalizedDerivative p (degree k) (Pi.single k 1) = 0 ∧
      (Pi.single k 1 : Fin p → ℝ) ≠ 0 := by
  constructor
  · rw [kernel_iff]
    intro i hi
    have hik : i ≠ k := by
      intro h
      exact hi (congrArg degree h)
    simp [Pi.single_eq_of_ne hik]
  · intro h
    have heq := congrFun h k
    norm_num at heq

theorem injective_iff_nonresonant {p : ℕ} (α : ℝ) :
    Function.Injective (normalizedDerivative p α) ↔
      ∀ i : Fin p, degree i ≠ α := by
  constructor
  · intro h k hk
    have hv := resonant_kernel_vector k
    have hz : normalizedDerivative p α (Pi.single k 1) = normalizedDerivative p α 0 := by
      simpa [← hk] using hv.1
    exact hv.2 (h hz)
  · exact injective_of_nonresonant α

/-- At a resonant degree the whole kernel, rather than only one witness, is a line. -/
theorem resonant_kernel {p : ℕ} (k : Fin p) :
    LinearMap.ker (normalizedDerivative p (degree k)) =
      Submodule.span ℝ {Pi.single k 1} := by
  ext c
  rw [LinearMap.mem_ker, kernel_iff, Submodule.mem_span_singleton]
  constructor
  · intro hc
    refine ⟨c k, ?_⟩
    funext i
    by_cases hik : i = k
    · subst i
      simp
    · have hdeg : degree i ≠ degree k := fun h => hik (degree_injective p h)
      simp [Pi.single_eq_of_ne hik, hc i hdeg]
  · rintro ⟨a, rfl⟩ i hi
    have hik : i ≠ k := fun h => hi (congrArg degree h)
    simp [Pi.single_eq_of_ne hik]

theorem resonant_nullity {p : ℕ} (k : Fin p) :
    Module.finrank ℝ (LinearMap.ker (normalizedDerivative p (degree k))) = 1 := by
  rw [resonant_kernel]
  exact finrank_span_singleton (resonant_kernel_vector k).2

theorem nonresonant_rank {p : ℕ} (α : ℝ) (hα : ∀ i : Fin p, degree i ≠ α) :
    Module.finrank ℝ (LinearMap.range (normalizedDerivative p α)) = p := by
  rw [LinearMap.finrank_range_of_inj (injective_of_nonresonant α hα)]
  simp

theorem resonant_rank {p : ℕ} (k : Fin p) :
    Module.finrank ℝ (LinearMap.range (normalizedDerivative p (degree k))) = p - 1 := by
  have hr := (normalizedDerivative p (degree k)).finrank_range_add_finrank_ker
  rw [resonant_nullity] at hr
  simpa using Nat.eq_sub_of_add_eq hr

/-- All physical GLISS preset weights are in this injective range. -/
theorem physical_injective {p : ℕ} {α : ℝ} (hα : α ≤ 1 / 2) :
    Function.Injective (normalizedDerivative p α) := by
  apply injective_of_nonresonant α
  intro i hi
  have hnonneg : (0 : ℝ) ≤ i.val := Nat.cast_nonneg _
  simp only [degree, Nat.cast_add, Nat.cast_one] at hi
  linarith

theorem physical_rank {p : ℕ} {α : ℝ} (hα : α ≤ 1 / 2) :
    Module.finrank ℝ (LinearMap.range (normalizedDerivative p α)) = p := by
  rw [LinearMap.finrank_range_of_inj (physical_injective hα)]
  simp

/-- The polynomial represented by the first-cell coefficient model. -/
def firstPolynomial {p : ℕ} (c : Fin p → ℝ) : Polynomial ℝ :=
  ∑ i, Polynomial.monomial (i.val + 1) (c i)

/-- Formal polynomial differentiation gives precisely the implemented coefficient map. -/
theorem polynomial_derivative_bridge {p : ℕ} (α : ℝ) (c : Fin p → ℝ) :
    Polynomial.X * (firstPolynomial c).derivative - Polynomial.C α * firstPolynomial c =
      firstPolynomial (normalizedDerivative p α c) := by
  classical
  simp only [firstPolynomial, map_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro i _
  simp only [Polynomial.derivative_monomial, Nat.add_sub_cancel,
    Polynomial.X_mul_monomial, Polynomial.C_mul_monomial,
    ← Polynomial.monomial_sub, normalizedDerivative_apply, degree]
  congr 1
  ring

/-- Distinct positive monomial degrees make the polynomial representation faithful. -/
theorem firstPolynomial_coefficient {p : ℕ} (c : Fin p → ℝ) (i : Fin p) :
    (firstPolynomial c).coeff (i.val + 1) = c i := by
  classical
  simp only [firstPolynomial, Polynomial.finset_sum_coeff]
  rw [Finset.sum_eq_single i]
  · simp
  · intro j _ hji
    have hne : j.val + 1 ≠ i.val + 1 := by
      intro h
      apply hji
      apply Fin.ext
      omega
    simp [Polynomial.coeff_monomial, hne]
  · simp

theorem firstPolynomial_injective (p : ℕ) : Function.Injective (@firstPolynomial p) := by
  intro c d h
  funext i
  have hi := congrArg (fun f : Polynomial ℝ => f.coeff (i.val + 1)) h
  simpa only [firstPolynomial_coefficient] using hi

theorem polynomial_kernel_iff {p : ℕ} (α : ℝ) (c : Fin p → ℝ) :
    Polynomial.X * (firstPolynomial c).derivative - Polynomial.C α * firstPolynomial c = 0 ↔
      normalizedDerivative p α c = 0 := by
  rw [polynomial_derivative_bridge]
  constructor
  · intro h
    apply firstPolynomial_injective p
    simpa [firstPolynomial] using h
  · intro h
    rw [h]
    simp [firstPolynomial]

/-- On an interior cell `N(s) = a + b*s` need not have zero constant term. -/
def affineNormalizedDerivative (α a b s : ℝ) : ℝ := -α * a + (1 - α) * b * s

/-- The affine-cell obstruction uses the same formal differential operator. -/
theorem affine_polynomial_derivative_bridge (α a b s : ℝ) :
    (Polynomial.X * (Polynomial.C a + Polynomial.C b * Polynomial.X).derivative -
      Polynomial.C α * (Polynomial.C a + Polynomial.C b * Polynomial.X)).eval s =
      affineNormalizedDerivative α a b s := by
  simp [affineNormalizedDerivative]
  ring

/-- The old weighted-P0 ansatz gives normalized value `e*s`. Even agreement
at two distinct positive interior points forces `α*a = 0`. -/
theorem degree_one_interior_obstruction {α a b e x y : ℝ}
    (hxy : x ≠ y)
    (hx : affineNormalizedDerivative α a b x = e * x)
    (hy : affineNormalizedDerivative α a b y = e * y) :
    α * a = 0 := by
  dsimp [affineNormalizedDerivative] at hx hy
  have hproduct : ((1 - α) * b - e) * (x - y) = 0 := by
    nlinarith [hx, hy]
  have hslope : (1 - α) * b - e = 0 :=
    (mul_eq_zero.mp hproduct).resolve_right (sub_ne_zero.mpr hxy)
  have hslope' : (1 - α) * b = e := sub_eq_zero.mp hslope
  rw [hslope'] at hx
  linarith

theorem degree_one_interior_no_old_companion {α a b : ℝ} (hα : α ≠ 0) (ha : a ≠ 0)
    {x y : ℝ} (hxy : x ≠ y) :
    ¬∃ e : ℝ, affineNormalizedDerivative α a b x = e * x ∧
      affineNormalizedDerivative α a b y = e * y := by
  rintro ⟨e, hx, hy⟩
  exact (mul_ne_zero hα ha) (degree_one_interior_obstruction hxy hx hy)

theorem tent_sample_points {h : ℝ} (hpos : 0 < h) (hsmall : h < 1) :
    0 < h ∧ h < h + (1 - h) / 3 ∧
      h + (1 - h) / 3 < h + 2 * (1 - h) / 3 ∧ h + 2 * (1 - h) / 3 < 1 := by
  exact ⟨hpos, by linarith, by linarith, by linarith⟩

/-- The descending cell of a two-cell boundary-zero tent is not an axis-only artifact.
`tent_sample_points` places both test points strictly between `h` and one. -/
theorem degree_one_tent_obstruction {h : ℝ} (hpos : 0 < h) (hsmall : h < 1) :
    ¬∃ e : ℝ,
      affineNormalizedDerivative (1 / 2) (h / (1 - h)) (-h / (1 - h))
        (h + (1 - h) / 3) = e * (h + (1 - h) / 3) ∧
      affineNormalizedDerivative (1 / 2) (h / (1 - h)) (-h / (1 - h))
        (h + 2 * (1 - h) / 3) = e * (h + 2 * (1 - h) / 3) := by
  apply degree_one_interior_no_old_companion (by norm_num : (1 / 2 : ℝ) ≠ 0)
  · exact div_ne_zero (ne_of_gt hpos) (by linarith)
  · linarith

end WeightedFeec

#print axioms WeightedFeec.resonant_nullity
#print axioms WeightedFeec.resonant_rank
#print axioms WeightedFeec.physical_injective
#print axioms WeightedFeec.polynomial_kernel_iff
#print axioms WeightedFeec.degree_one_tent_obstruction
