import Mathlib
open Complex

namespace S3Zeta

/-- Spectral data of -Δ+1 on the round S³: eigenvalue k² with multiplicity k², k ≥ 1. -/
noncomputable def term (s : ℂ) (k : ℕ) : ℂ := (k : ℂ) ^ 2 * ((k : ℂ) ^ 2) ^ (-s)

lemma term_eq (s : ℂ) (k : ℕ) (hk : 1 ≤ k) :
    term s k = 1 / (k : ℂ) ^ (2 * s - 2) := by
  have hk0 : (k : ℂ) ≠ 0 := by exact_mod_cast (by omega : k ≠ 0)
  have hkpos : (0 : ℝ) < k := by exact_mod_cast hk
  have hlog : Complex.log ((k : ℂ)) = ((Real.log k : ℝ) : ℂ) := by
    have : ((k : ℂ)) = ((k : ℝ) : ℂ) := by simp
    rw [this, ← Complex.ofReal_log hkpos.le]
  have hsq : (k : ℂ) ^ 2 ≠ 0 := pow_ne_zero 2 hk0
  have hlog2 : Complex.log ((k : ℂ) ^ 2) = 2 * Complex.log (k : ℂ) := by
    have : ((k : ℂ) ^ 2) = ((k:ℝ)^2 : ℝ) := by simp
    rw [this, ← Complex.ofReal_log (by positivity), Real.log_pow, hlog]; push_cast; ring
  unfold term
  rw [Complex.cpow_def_of_ne_zero hsq, Complex.cpow_def_of_ne_zero hk0, hlog2]
  have h2 : (k : ℂ) ^ 2 = Complex.exp (2 * Complex.log (k : ℂ)) := by
    have := Complex.exp_nat_mul (Complex.log (k : ℂ)) 2
    rw [Complex.exp_log hk0] at this
    simpa using this.symm
  rw [h2, ← Complex.exp_add, one_div, ← Complex.exp_neg]
  congr 1
  ring

/-- On Re s > 3/2 the spectral zeta function of -Δ+1 on S³ equals ζ(2s-2). -/
theorem spectral_zeta_eq (s : ℂ) (hs : 3 / 2 < s.re) :
    ∑' k : ℕ, term s k = riemannZeta (2 * s - 2) := by
  have hre : 1 < (2 * s - 2).re := by simp; linarith
  rw [zeta_eq_tsum_one_div_nat_cpow hre]
  congr 1; ext k
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · have hne : (2 * s - 2) ≠ 0 := by
      intro h; have := hre; rw [h] at this; norm_num at this
    simp [term, Complex.zero_cpow hne]
  · exact term_eq s k hk

/-- RH (Mathlib's statement) is equivalent to: every nontrivial zero of the continued spectral
zeta function ζ(2s-2) lies on Re s = 5/4. -/
theorem RH_iff_spectral :
    RiemannHypothesis ↔
      ∀ s : ℂ, riemannZeta (2 * s - 2) = 0 → (¬ ∃ n : ℕ, 2 * s - 2 = -2 * (n + 1)) →
        2 * s - 2 ≠ 1 → s.re = 5 / 4 := by
  constructor
  · intro h s h0 hn h1
    have := h (2 * s - 2) h0 hn h1
    simp at this; linarith
  · intro h t h0 hn h1
    have := h ((t + 2) / 2) (by ring_nf; ring_nf at h0; simpa using h0)
      (by ring_nf; ring_nf at hn; simpa using hn) (by ring_nf; ring_nf at h1; simpa using h1)
    simp at this; linarith

end S3Zeta
