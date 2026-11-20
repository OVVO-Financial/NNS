-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.MeasurePM
public import NNSProofs.Moments
public import NNSProofs.CoMomentsMatrix
public import Mathlib.Probability.Distributions.Gaussian.Real

namespace NNS

open MeasureTheory ProbabilityTheory Finset
open scoped ENNReal

@[expose] public noncomputable section

namespace GapProjection

variable {n : ℕ} {μ ν : Measure ℝ}

theorem integral_id_empMeasure (x : Fin n → ℝ) : ∫ v, v ∂empMeasure x = mean x :=
  integral_empMeasure x fun v => v

theorem integral_sub_pow_empMeasure (r : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    ∫ v, (v - t) ^ r ∂empMeasure x = centralMoment r t x :=
  integral_empMeasure x fun v => (v - t) ^ r

theorem variance_empMeasure (hn : 0 < n) (x : Fin n → ℝ) :
    Var[fun v : ℝ => v; empMeasure x] = popVariance x := by
  have := isProbabilityMeasure_empMeasure hn x
  rw [variance_eq_mUPM_two_add_mLPM_two (memLp_empMeasure hn x 2), integral_id_empMeasure,
    mUPM_empMeasure, mLPM_empMeasure, popVariance_eq_PM, add_comm]

def symProj (r : ℕ) (p : ℝ × ℝ) : ℝ := p.1 + (-1) ^ r * p.2

theorem integral_sub_pow_eq_symProj [IsFiniteMeasure μ] {r : ℕ}
    (h : MemLp (fun x : ℝ => x) r μ) (t : ℝ) :
    ∫ x, (x - t) ^ r ∂μ = symProj r (mUPM r t μ, mLPM r t μ) :=
  (mUPM_add_neg_one_pow_mul_mLPM h t).symm

theorem integral_sub_pow_eq_of_pair_eq [IsFiniteMeasure μ] [IsFiniteMeasure ν] {r : ℕ}
    (hμ : MemLp (fun x : ℝ => x) r μ) (hν : MemLp (fun x : ℝ => x) r ν) {s t : ℝ}
    (h : (mUPM r s μ, mLPM r s μ) = (mUPM r t ν, mLPM r t ν)) :
    ∫ x, (x - s) ^ r ∂μ = ∫ x, (x - t) ^ r ∂ν := by
  rw [integral_sub_pow_eq_symProj hμ, integral_sub_pow_eq_symProj hν, h]

theorem variance_eq_symProj [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ) :
    Var[fun x : ℝ => x; μ] = symProj 2 (mUPM 2 (∫ x, x ∂μ) μ, mLPM 2 (∫ x, x ∂μ) μ) := by
  rw [variance_eq_mUPM_two_add_mLPM_two h]
  simp [symProj]

theorem variance_eq_of_pair_eq [IsProbabilityMeasure μ] [IsProbabilityMeasure ν]
    (hμ : MemLp (fun x : ℝ => x) 2 μ) (hν : MemLp (fun x : ℝ => x) 2 ν)
    (h : (mUPM 2 (∫ x, x ∂μ) μ, mLPM 2 (∫ x, x ∂μ) μ)
      = (mUPM 2 (∫ x, x ∂ν) ν, mLPM 2 (∫ x, x ∂ν) ν)) :
    Var[fun x : ℝ => x; μ] = Var[fun x : ℝ => x; ν] := by
  rw [variance_eq_symProj hμ, variance_eq_symProj hν, h]

lemma upperDev_add (r : ℕ) (t a : ℝ) : upperDev r t (t + a) = if 0 < a then a ^ r else 0 := by
  simp [upperDev]

lemma lowerDev_add (r : ℕ) (t a : ℝ) :
    lowerDev r t (t + a) = if a ≤ 0 then (-a) ^ r else 0 := by
  simp [lowerDev]

theorem exists_samples_moment_eq_pair_ne (r : ℕ) (t : ℝ) :
    ∃ x y : Fin 2 → ℝ, centralMoment r t x = centralMoment r t y ∧
      (UPM r t x, LPM r t x) ≠ (UPM r t y, LPM r t y) := by
  rcases Nat.even_or_odd r with hr | hr
  · refine ⟨![t + -1, t + 1], ![t + 1, t + 1], ?_, ?_⟩
    · simp [centralMoment, Fin.sum_univ_two, hr.neg_one_pow]
    · have h1 : LPM r t ![t + -1, t + 1] = 1 / 2 := by
        norm_num [LPM, Fin.sum_univ_two, lowerDev_add]
      have h2 : LPM r t ![t + 1, t + 1] = 0 := by
        norm_num [LPM, Fin.sum_univ_two, lowerDev_add]
      intro h
      have := congrArg Prod.snd h
      simp only [h1, h2] at this
      norm_num at this
  · have hr0 : r ≠ 0 := by
      rintro rfl
      exact absurd hr (by decide)
    have h2r : (1 : ℝ) < 2 ^ r := one_lt_pow₀ one_lt_two hr0
    refine ⟨![t + -1, t + 1], ![t + -2, t + 2], ?_, ?_⟩
    · simp [centralMoment, Fin.sum_univ_two, hr.neg_one_pow, hr.neg_pow]
    · have h1 : UPM r t ![t + -1, t + 1] = 1 / 2 := by
        norm_num [UPM, Fin.sum_univ_two, upperDev_add]
      have h2 : UPM r t ![t + -2, t + 2] = 2 ^ r / 2 := by
        norm_num [UPM, Fin.sum_univ_two, upperDev_add]
      intro h
      have := congrArg Prod.fst h
      simp only [h1, h2] at this
      linarith

theorem exists_laws_moment_eq_pair_ne (r : ℕ) (t : ℝ) :
    ∃ μ ν : Measure ℝ, IsProbabilityMeasure μ ∧ IsProbabilityMeasure ν ∧
      MemLp (fun x : ℝ => x) r μ ∧ MemLp (fun x : ℝ => x) r ν ∧
      ∫ x, (x - t) ^ r ∂μ = ∫ x, (x - t) ^ r ∂ν ∧
      (mUPM r t μ, mLPM r t μ) ≠ (mUPM r t ν, mLPM r t ν) := by
  obtain ⟨x, y, hm, hp⟩ := exists_samples_moment_eq_pair_ne r t
  refine ⟨empMeasure x, empMeasure y, isProbabilityMeasure_empMeasure two_pos x,
    isProbabilityMeasure_empMeasure two_pos y, memLp_empMeasure two_pos x r,
    memLp_empMeasure two_pos y r, ?_, ?_⟩
  · rw [integral_sub_pow_empMeasure, integral_sub_pow_empMeasure, hm]
  · simpa only [mUPM_empMeasure, mLPM_empMeasure] using hp

theorem pair_partition_strictly_refines (r : ℕ) (t : ℝ) :
    (∀ μ ν : Measure ℝ, IsProbabilityMeasure μ → IsProbabilityMeasure ν →
      MemLp (fun x : ℝ => x) r μ → MemLp (fun x : ℝ => x) r ν →
      (mUPM r t μ, mLPM r t μ) = (mUPM r t ν, mLPM r t ν) →
      ∫ x, (x - t) ^ r ∂μ = ∫ x, (x - t) ^ r ∂ν) ∧
    ¬ ∀ μ ν : Measure ℝ, IsProbabilityMeasure μ → IsProbabilityMeasure ν →
      MemLp (fun x : ℝ => x) r μ → MemLp (fun x : ℝ => x) r ν →
      ∫ x, (x - t) ^ r ∂μ = ∫ x, (x - t) ^ r ∂ν →
      (mUPM r t μ, mLPM r t μ) = (mUPM r t ν, mLPM r t ν) := by
  refine ⟨fun μ ν _ _ hμ hν h => integral_sub_pow_eq_of_pair_eq hμ hν h, fun hall => ?_⟩
  obtain ⟨μ, ν, hPμ, hPν, hμ, hν, hm, hp⟩ := exists_laws_moment_eq_pair_ne r t
  exact hp (hall μ ν hPμ hPν hμ hν hm)

theorem not_exists_recover_pair (r : ℕ) (t : ℝ) :
    ¬ ∃ g : ℝ → ℝ × ℝ, ∀ μ : Measure ℝ, IsProbabilityMeasure μ →
      MemLp (fun x : ℝ => x) r μ → g (∫ x, (x - t) ^ r ∂μ) = (mUPM r t μ, mLPM r t μ) := by
  rintro ⟨g, hg⟩
  obtain ⟨μ, ν, hPμ, hPν, hμ, hν, hm, hp⟩ := exists_laws_moment_eq_pair_ne r t
  exact hp (by rw [← hg μ hPμ hμ, ← hg ν hPν hν, hm])

def tableA : Fin 36 → ℝ := fun i => if i.val = 0 then 18 else if i.val ≤ 9 then -2 else 0

def tableB : Fin 5 → ℝ := ![5, -5, 0, 0, 0]

lemma sum_tableA (g : ℝ → ℝ) : ∑ i, g (tableA i) = g 18 + 9 * g (-2) + 26 * g 0 := by
  have h := Fin.sum_univ_eq_sum_range
    (fun k => g (if k = 0 then 18 else if k ≤ 9 then -2 else 0)) 36
  simp only [tableA]
  rw [h]
  simp [Finset.sum_range_succ]
  ring

lemma sum_tableB (g : ℝ → ℝ) : ∑ i, g (tableB i) = g 5 + g (-5) + 3 * g 0 := by
  simp [tableB, Fin.sum_univ_succ]
  ring

theorem mean_tableA : mean tableA = 0 := by
  have h := sum_tableA fun v => v
  simp only [mean, h]
  norm_num

theorem mean_tableB : mean tableB = 0 := by
  have h := sum_tableB fun v => v
  simp only [mean, h]
  norm_num

theorem table_A :
    UPM 2 (mean tableA) tableA = 9 ∧ LPM 2 (mean tableA) tableA = 1 ∧
      popVariance tableA = 10 := by
  have hU : UPM 2 (mean tableA) tableA = 9 := by
    rw [mean_tableA, UPM, sum_tableA (upperDev 2 0)]
    norm_num [upperDev]
  have hL : LPM 2 (mean tableA) tableA = 1 := by
    rw [mean_tableA, LPM, sum_tableA (lowerDev 2 0)]
    norm_num [lowerDev]
  refine ⟨hU, hL, ?_⟩
  rw [popVariance_eq_PM, hU, hL]
  norm_num

theorem table_B :
    UPM 2 (mean tableB) tableB = 5 ∧ LPM 2 (mean tableB) tableB = 5 ∧
      popVariance tableB = 10 := by
  have hU : UPM 2 (mean tableB) tableB = 5 := by
    rw [mean_tableB, UPM, sum_tableB (upperDev 2 0)]
    norm_num [upperDev]
  have hL : LPM 2 (mean tableB) tableB = 5 := by
    rw [mean_tableB, LPM, sum_tableB (lowerDev 2 0)]
    norm_num [lowerDev]
  refine ⟨hU, hL, ?_⟩
  rw [popVariance_eq_PM, hU, hL]
  norm_num

theorem table_laws :
    (mUPM 2 (∫ x, x ∂empMeasure tableA) (empMeasure tableA) = 9 ∧
      mLPM 2 (∫ x, x ∂empMeasure tableA) (empMeasure tableA) = 1 ∧
      Var[fun x : ℝ => x; empMeasure tableA] = 10) ∧
    (mUPM 2 (∫ x, x ∂empMeasure tableB) (empMeasure tableB) = 5 ∧
      mLPM 2 (∫ x, x ∂empMeasure tableB) (empMeasure tableB) = 5 ∧
      Var[fun x : ℝ => x; empMeasure tableB] = 10) := by
  simp only [integral_id_empMeasure, mUPM_empMeasure, mLPM_empMeasure,
    variance_empMeasure (by norm_num : 0 < 36), variance_empMeasure (by norm_num : 0 < 5)]
  exact ⟨table_A, table_B⟩

theorem exists_laws_variance_eq_pair_ne :
    ∃ μ ν : Measure ℝ, IsProbabilityMeasure μ ∧ IsProbabilityMeasure ν ∧
      MemLp (fun x : ℝ => x) 2 μ ∧ MemLp (fun x : ℝ => x) 2 ν ∧
      Var[fun x : ℝ => x; μ] = Var[fun x : ℝ => x; ν] ∧
      (mUPM 2 (∫ x, x ∂μ) μ, mLPM 2 (∫ x, x ∂μ) μ)
        ≠ (mUPM 2 (∫ x, x ∂ν) ν, mLPM 2 (∫ x, x ∂ν) ν) := by
  obtain ⟨⟨hAU, hAL, hAV⟩, hBU, hBL, hBV⟩ := table_laws
  refine ⟨empMeasure tableA, empMeasure tableB,
    isProbabilityMeasure_empMeasure (by norm_num) _,
    isProbabilityMeasure_empMeasure (by norm_num) _, memLp_empMeasure (by norm_num) _ 2,
    memLp_empMeasure (by norm_num) _ 2, by rw [hAV, hBV], ?_⟩
  rw [hAU, hAL, hBU, hBL]
  norm_num

theorem pair_partition_strictly_refines_variance :
    (∀ μ ν : Measure ℝ, IsProbabilityMeasure μ → IsProbabilityMeasure ν →
      MemLp (fun x : ℝ => x) 2 μ → MemLp (fun x : ℝ => x) 2 ν →
      (mUPM 2 (∫ x, x ∂μ) μ, mLPM 2 (∫ x, x ∂μ) μ)
        = (mUPM 2 (∫ x, x ∂ν) ν, mLPM 2 (∫ x, x ∂ν) ν) →
      Var[fun x : ℝ => x; μ] = Var[fun x : ℝ => x; ν]) ∧
    ¬ ∀ μ ν : Measure ℝ, IsProbabilityMeasure μ → IsProbabilityMeasure ν →
      MemLp (fun x : ℝ => x) 2 μ → MemLp (fun x : ℝ => x) 2 ν →
      Var[fun x : ℝ => x; μ] = Var[fun x : ℝ => x; ν] →
      (mUPM 2 (∫ x, x ∂μ) μ, mLPM 2 (∫ x, x ∂μ) μ)
        = (mUPM 2 (∫ x, x ∂ν) ν, mLPM 2 (∫ x, x ∂ν) ν) := by
  refine ⟨fun μ ν _ _ hμ hν h => variance_eq_of_pair_eq hμ hν h, fun hall => ?_⟩
  obtain ⟨μ, ν, hPμ, hPν, hμ, hν, hv, hp⟩ := exists_laws_variance_eq_pair_ne
  exact hp (hall μ ν hPμ hPν hμ hν hv)

theorem not_exists_recover_pair_of_variance :
    ¬ ∃ g : ℝ → ℝ × ℝ, ∀ μ : Measure ℝ, IsProbabilityMeasure μ → MemLp (fun x : ℝ => x) 2 μ →
      g Var[fun x : ℝ => x; μ] = (mUPM 2 (∫ x, x ∂μ) μ, mLPM 2 (∫ x, x ∂μ) μ) := by
  rintro ⟨g, hg⟩
  obtain ⟨μ, ν, hPμ, hPν, hμ, hν, hv, hp⟩ := exists_laws_variance_eq_pair_ne
  exact hp (by rw [← hg μ hPμ hμ, ← hg ν hPν hν, hv])

lemma upperDev_neg_neg {r : ℕ} (hr : r ≠ 0) (t v : ℝ) :
    upperDev r (-t) (-v) = lowerDev r t v := by
  simp only [upperDev, lowerDev, neg_lt_neg_iff, neg_sub_neg]
  rcases lt_trichotomy v t with h | h | h
  · simp [h, h.le]
  · simp [h, zero_pow hr]
  · simp [not_lt.mpr h.le, not_le.mpr h]

lemma lowerDev_neg_neg {r : ℕ} (hr : r ≠ 0) (t v : ℝ) :
    lowerDev r (-t) (-v) = upperDev r t v := by
  simp only [upperDev, lowerDev, neg_le_neg_iff, neg_sub_neg]
  rcases lt_trichotomy v t with h | h | h
  · simp [not_lt.mpr h.le, not_le.mpr h]
  · simp [h, zero_pow hr]
  · simp [h, h.le]

theorem mUPM_map_neg {r : ℕ} (hr : r ≠ 0) (t : ℝ) (μ : Measure ℝ) :
    mUPM r (-t) (μ.map fun x => -x) = mLPM r t μ := by
  rw [mUPM, integral_map measurable_neg.aemeasurable
    (upperDev_measurable r (-t)).aestronglyMeasurable]
  exact integral_congr_ae (.of_forall fun v => upperDev_neg_neg hr t v)

theorem mLPM_map_neg {r : ℕ} (hr : r ≠ 0) (t : ℝ) (μ : Measure ℝ) :
    mLPM r (-t) (μ.map fun x => -x) = mUPM r t μ := by
  rw [mLPM, integral_map measurable_neg.aemeasurable
    (lowerDev_measurable r (-t)).aestronglyMeasurable]
  exact integral_congr_ae (.of_forall fun v => lowerDev_neg_neg hr t v)

theorem integral_id_map_neg (μ : Measure ℝ) : ∫ x, x ∂μ.map (fun x => -x) = -∫ x, x ∂μ := by
  rw [integral_map (φ := fun x : ℝ => -x) (f := fun x : ℝ => x) (by fun_prop) (by fun_prop),
    integral_neg]

theorem variance_map_neg (μ : Measure ℝ) :
    Var[fun x : ℝ => x; μ.map fun x => -x] = Var[fun x : ℝ => x; μ] := by
  rw [variance_map aemeasurable_id' measurable_neg.aemeasurable]
  exact variance_fun_neg

theorem reflection_swaps_pair_at_mean (μ : Measure ℝ) :
    mUPM 2 (∫ x, x ∂μ.map fun x => -x) (μ.map fun x => -x) = mLPM 2 (∫ x, x ∂μ) μ ∧
    mLPM 2 (∫ x, x ∂μ.map fun x => -x) (μ.map fun x => -x) = mUPM 2 (∫ x, x ∂μ) μ ∧
    Var[fun x : ℝ => x; μ.map fun x => -x] = Var[fun x : ℝ => x; μ] := by
  rw [integral_id_map_neg]
  exact ⟨mUPM_map_neg two_ne_zero _ μ, mLPM_map_neg two_ne_zero _ μ, variance_map_neg μ⟩

theorem exists_laws_variance_eq_pair_swapped :
    ∃ μ ν : Measure ℝ, IsProbabilityMeasure μ ∧ IsProbabilityMeasure ν ∧
      MemLp (fun x : ℝ => x) 2 μ ∧ MemLp (fun x : ℝ => x) 2 ν ∧
      Var[fun x : ℝ => x; μ] = 10 ∧ Var[fun x : ℝ => x; ν] = 10 ∧
      mUPM 2 (∫ x, x ∂μ) μ = 9 ∧ mLPM 2 (∫ x, x ∂μ) μ = 1 ∧
      mUPM 2 (∫ x, x ∂ν) ν = 1 ∧ mLPM 2 (∫ x, x ∂ν) ν = 9 := by
  obtain ⟨⟨hAU, hAL, hAV⟩, -⟩ := table_laws
  obtain ⟨hU, hL, hV⟩ := reflection_swaps_pair_at_mean (empMeasure tableA)
  have hP := isProbabilityMeasure_empMeasure (by norm_num : 0 < 36) tableA
  have hm : MemLp (fun x : ℝ => x) 2 (empMeasure tableA) := memLp_empMeasure (by norm_num) _ 2
  refine ⟨empMeasure tableA, (empMeasure tableA).map fun x => -x, hP,
    (Measure.isProbabilityMeasure_map_iff measurable_neg.aemeasurable).mpr hP, hm, ?_, hAV,
    by rw [hV, hAV], hAU, hAL, by rw [hU, hAL], by rw [hL, hAU]⟩
  exact (memLp_map_measure_iff aestronglyMeasurable_id measurable_neg.aemeasurable).mpr hm.neg

section Correlation

variable {Ω : Type*} {mΩ : MeasurableSpace Ω} {P : Measure Ω} {X Y : Ω → ℝ}

def corr (X Y : Ω → ℝ) (P : Measure Ω) : ℝ :=
  cov[X, Y; P] / (Real.sqrt Var[X; P] * Real.sqrt Var[Y; P])

theorem covariance_sq_le_variance_mul_variance [IsFiniteMeasure P] (hX : MemLp X 2 P)
    (hY : MemLp Y 2 P) : cov[X, Y; P] ^ 2 ≤ Var[X; P] * Var[Y; P] := by
  have h : ∀ s : ℝ, 0 ≤ Var[Y; P] * (s * s) + 2 * cov[X, Y; P] * s + Var[X; P] := by
    intro s
    have hs := variance_nonneg (X + s • Y) P
    rw [variance_add hX (hY.const_smul s), covariance_smul_right, variance_smul] at hs
    linarith
  have hd := discrim_le_zero h
  rw [discrim] at hd
  nlinarith

theorem abs_corr_le_one [IsFiniteMeasure P] (hX : MemLp X 2 P) (hY : MemLp Y 2 P) :
    |corr X Y P| ≤ 1 := by
  have h := Real.abs_le_sqrt (covariance_sq_le_variance_mul_variance hX hY)
  rw [Real.sqrt_mul (variance_nonneg X P)] at h
  have hd : 0 ≤ Real.sqrt Var[X; P] * Real.sqrt Var[Y; P] :=
    mul_nonneg (Real.sqrt_nonneg _) (Real.sqrt_nonneg _)
  rw [corr, abs_div, abs_of_nonneg hd]
  exact div_le_one_of_le₀ h hd

theorem corr_mem_Icc [IsFiniteMeasure P] (hX : MemLp X 2 P) (hY : MemLp Y 2 P) :
    -1 ≤ corr X Y P ∧ corr X Y P ≤ 1 :=
  abs_le.mp (abs_corr_le_one hX hY)

end Correlation

def sampleCorr (x y : Fin n → ℝ) : ℝ :=
  covariance x y / (Real.sqrt (popVariance x) * Real.sqrt (popVariance y))

theorem sample_covariance_self (x : Fin n → ℝ) : covariance x x = popVariance x := by
  simp only [covariance, coMoment, popVariance, centralMoment, sq]

theorem sample_covariance_sq_le (x y : Fin n → ℝ) :
    covariance x y ^ 2 ≤ popVariance x * popVariance y :=
  sq_sum_div_le (fun i => x i - mean x) fun i => y i - mean y

theorem abs_sampleCorr_le_one (x y : Fin n → ℝ) : |sampleCorr x y| ≤ 1 := by
  have h := Real.abs_le_sqrt (sample_covariance_sq_le x y)
  rw [Real.sqrt_mul (popVariance_nonneg x)] at h
  have hd : 0 ≤ Real.sqrt (popVariance x) * Real.sqrt (popVariance y) :=
    mul_nonneg (Real.sqrt_nonneg _) (Real.sqrt_nonneg _)
  rw [sampleCorr, abs_div, abs_of_nonneg hd]
  exact div_le_one_of_le₀ h hd

section Matrix

variable {k : ℕ}

def corrMatrix (S : Matrix (Fin k) (Fin k) ℝ) : Matrix (Fin k) (Fin k) ℝ :=
  Matrix.of fun i j => S i j / Real.sqrt (S i i * S j j)

theorem corrMatrix_popCovMatrix (X : Matrix (Fin n) (Fin k) ℝ) (i j : Fin k) :
    corrMatrix (popCovMatrix X) i j = sampleCorr (colOf X i) (colOf X j) := by
  simp only [corrMatrix, popCovMatrix, Matrix.of_apply, sample_covariance_self, sampleCorr]
  rw [Real.sqrt_mul (popVariance_nonneg _)]

end Matrix

def unifFin (n : ℕ) : Measure (Fin n) := (n : ℝ≥0∞)⁻¹ • Measure.count

theorem integral_unifFin (f : Fin n → ℝ) : ∫ i, f i ∂unifFin n = (∑ i, f i) / n := by
  rw [unifFin, integral_smul_measure, integral_count]
  simp [div_eq_inv_mul]

theorem covariance_unifFin (x y : Fin n → ℝ) : cov[x, y; unifFin n] = covariance x y := by
  simp only [ProbabilityTheory.covariance, integral_unifFin]
  rfl

theorem variance_unifFin (x : Fin n → ℝ) : Var[x; unifFin n] = popVariance x := by
  rw [← covariance_self Measurable.of_discrete.aemeasurable, covariance_unifFin,
    sample_covariance_self]

theorem sampleCorr_eq_corr (x y : Fin n → ℝ) : sampleCorr x y = corr x y (unifFin n) := by
  rw [corr, covariance_unifFin, variance_unifFin, variance_unifFin, sampleCorr]

section Nonlinear

variable {Ω : Type*} {mΩ : MeasurableSpace Ω} {P : Measure Ω} {X : Ω → ℝ}

theorem integral_odd_eq_zero_of_symm (hX : AEMeasurable X P)
    (hsym : P.map (fun ω => -X ω) = P.map X) {f : ℝ → ℝ} (hf : Measurable f)
    (hodd : ∀ v, f (-v) = -f v) : ∫ ω, f (X ω) ∂P = 0 := by
  have h1 : ∫ ω, f (X ω) ∂P = ∫ v, f v ∂P.map X :=
    (integral_map hX hf.aestronglyMeasurable).symm
  have h2 : ∫ v, f v ∂P.map X = ∫ ω, f (-X ω) ∂P := by
    rw [← hsym, integral_map (φ := fun ω => -X ω) hX.neg hf.aestronglyMeasurable]
  have h3 : ∫ ω, f (-X ω) ∂P = -∫ ω, f (X ω) ∂P := by
    simp only [hodd]
    exact integral_neg _
  linarith

lemma memLp_sq_of_memLp_four (hX : MemLp X 4 P) : MemLp (fun ω => X ω ^ 2) 2 P := by
  have h4 : MemLp X ((4 : ℕ) : ℝ≥0∞) P := by exact_mod_cast hX
  rw [memLp_two_iff_integrable_sq (hX.aemeasurable.pow_const 2).aestronglyMeasurable]
  refine (h4.integrable_norm_pow (by norm_num)).congr (.of_forall fun ω => ?_)
  change ‖X ω‖ ^ 4 = (X ω ^ 2) ^ 2
  rw [Real.norm_eq_abs, ← sq_abs (X ω)]
  ring

theorem covariance_sq_eq_zero_of_symm [IsProbabilityMeasure P]
    (hsym : P.map (fun ω => -X ω) = P.map X) (hX : MemLp X 4 P) :
    cov[X, fun ω => X ω ^ 2; P] = 0 := by
  have hm := hX.aemeasurable
  have hE1 : ∫ ω, X ω ∂P = 0 :=
    integral_odd_eq_zero_of_symm hm hsym (f := fun v => v) measurable_id fun _ => rfl
  have hE3 : ∫ ω, X ω ^ 3 ∂P = 0 :=
    integral_odd_eq_zero_of_symm hm hsym (f := fun v => v ^ 3) (by fun_prop) fun v => by ring
  have hX2 : MemLp X 2 P := hX.mono_exponent (by norm_num)
  have h3 : (X * fun ω => X ω ^ 2) = fun ω => X ω ^ 3 := by
    ext ω
    simp only [Pi.mul_apply]
    ring
  rw [covariance_eq_sub hX2 (memLp_sq_of_memLp_four hX), hE1, zero_mul, sub_zero, h3, hE3]

theorem corr_sq_eq_zero_of_symm [IsProbabilityMeasure P]
    (hsym : P.map (fun ω => -X ω) = P.map X) (hX : MemLp X 4 P) :
    corr X (fun ω => X ω ^ 2) P = 0 := by
  rw [corr, covariance_sq_eq_zero_of_symm hsym hX, zero_div]

theorem variance_pos_of_variance_sq_pos [IsProbabilityMeasure P] (hX : MemLp X 2 P)
    (hv : 0 < Var[fun ω => X ω ^ 2; P]) : 0 < Var[X; P] := by
  by_contra h
  have h0 : Var[X; P] = 0 := le_antisymm (not_lt.mp h) (variance_nonneg X P)
  rw [variance, ENNReal.toReal_eq_zero_iff] at h0
  have he : evariance X P = 0 := h0.resolve_right (evariance_lt_top hX).ne
  have hae := (evariance_eq_zero_iff hX.aemeasurable).mp he
  have hsq : (fun ω => X ω ^ 2) =ᵐ[P] fun _ => (∫ ω, X ω ∂P) ^ 2 :=
    hae.mono fun ω hω => by simp only [hω]
  rw [variance_congr hsq, ← covariance_self aemeasurable_const, covariance_const_left] at hv
  exact lt_irrefl _ hv

theorem corr_sq_gaussianReal :
    corr (fun x : ℝ => x) (fun x => x ^ 2) (gaussianReal 0 1) = 0 := by
  refine corr_sq_eq_zero_of_symm ?_ (memLp_id_gaussianReal' 4 (by norm_num))
  rw [gaussianReal_map_neg, neg_zero]
  exact Measure.map_id'.symm

end Nonlinear

theorem mean_neg (x : Fin n → ℝ) : mean (-x) = -mean x := by
  simp [mean, neg_div]

theorem mean_smul (c : ℝ) (x : Fin n → ℝ) : mean (c • x) = c * mean x := by
  simp [mean, ← mul_sum, mul_div_assoc]

theorem CoLPM_neg_neg {d e : ℕ} (hd : d ≠ 0) (he : e ≠ 0) (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoLPM d e (-tx) (-ty) (-x) (-y) = CoUPM d e tx ty x y := by
  simp only [CoLPM, CoUPM, Pi.neg_apply, lowerDev_neg_neg hd, lowerDev_neg_neg he]

theorem CoUPM_neg_neg {d e : ℕ} (hd : d ≠ 0) (he : e ≠ 0) (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoUPM d e (-tx) (-ty) (-x) (-y) = CoLPM d e tx ty x y := by
  simp only [CoLPM, CoUPM, Pi.neg_apply, upperDev_neg_neg hd, upperDev_neg_neg he]

theorem sample_covariance_neg_neg (x y : Fin n → ℝ) : covariance (-x) (-y) = covariance x y := by
  simp only [covariance, coMoment, mean_neg, Pi.neg_apply]
  congr 1
  exact sum_congr rfl fun i _ => by ring

section Matrix

variable {k : ℕ}

theorem PMmatrix_neg (X : Matrix (Fin n) (Fin k) ℝ) :
    (PMmatrix 1 1 (colMeans (-X)) (-X) false false).clpm
        = (PMmatrix 1 1 (colMeans X) X false false).cupm ∧
      (PMmatrix 1 1 (colMeans (-X)) (-X) false false).cupm
        = (PMmatrix 1 1 (colMeans X) X false false).clpm ∧
      (PMmatrix 1 1 (colMeans (-X)) (-X) false false).cov
        = (PMmatrix 1 1 (colMeans X) X false false).cov := by
  have hcol : ∀ i, colOf (-X) i = -colOf X i := fun i => rfl
  refine ⟨?_, ?_, ?_⟩
  · ext i j
    simp only [PMmatrix, Bool.false_eq_true, ↓reduceIte, clpmM, cupmM, Matrix.of_apply,
      colMeans, hcol, mean_neg, CoLPM_neg_neg one_ne_zero one_ne_zero]
  · ext i j
    simp only [PMmatrix, Bool.false_eq_true, ↓reduceIte, clpmM, cupmM, Matrix.of_apply,
      colMeans, hcol, mean_neg, CoUPM_neg_neg one_ne_zero one_ne_zero]
  · rw [PMmatrix_cov_eq_popCov, PMmatrix_cov_eq_popCov]
    ext i j
    simp only [popCovMatrix, Matrix.of_apply, hcol, sample_covariance_neg_neg]

end Matrix

theorem sample_covariance_smul_smul (c d : ℝ) (x y : Fin n → ℝ) :
    covariance (c • x) (d • y) = c * d * covariance x y := by
  simp only [covariance, coMoment, mean_smul, Pi.smul_apply, smul_eq_mul]
  rw [← mul_div_assoc, mul_sum]
  congr 1
  exact sum_congr rfl fun i _ => by ring

section Matrix

variable {k : ℕ}

theorem corrMatrix_scale (S : Matrix (Fin k) (Fin k) ℝ) {d : Fin k → ℝ} (hd : ∀ i, 0 < d i) :
    corrMatrix (Matrix.of fun i j => d i * d j * S i j) = corrMatrix S := by
  ext i j
  have hij : 0 < d i * d j := mul_pos (hd i) (hd j)
  have h : d i * d i * S i i * (d j * d j * S j j) = (d i * d j) ^ 2 * (S i i * S j j) := by
    ring
  simp only [corrMatrix, Matrix.of_apply]
  rw [h, Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hij.le, mul_div_mul_left _ _ hij.ne']

theorem popCovMatrix_scale (X : Matrix (Fin n) (Fin k) ℝ) (d : Fin k → ℝ) :
    popCovMatrix (Matrix.of fun r j => d j * X r j)
      = Matrix.of fun i j => d i * d j * popCovMatrix X i j := by
  ext i j
  have hcol : ∀ i, colOf (Matrix.of fun r j => d j * X r j) i = d i • colOf X i := fun i => rfl
  simp only [popCovMatrix, Matrix.of_apply, hcol, sample_covariance_smul_smul]

theorem exists_popCovMatrix_eq_clpm_ne :
    ∃ X X' : Matrix (Fin 3) (Fin 2) ℝ, popCovMatrix X = popCovMatrix X' ∧
      (PMmatrix 1 1 (colMeans X) X false false).clpm
        ≠ (PMmatrix 1 1 (colMeans X') X' false false).clpm := by
  set X : Matrix (Fin 3) (Fin 2) ℝ := Matrix.of fun r _ => (![2, -1, -1] : Fin 3 → ℝ) r
  obtain ⟨h1, -, h3⟩ := PMmatrix_neg X
  refine ⟨X, -X, ?_, fun h => ?_⟩
  · rw [← PMmatrix_cov_eq_popCov, ← PMmatrix_cov_eq_popCov, h3]
  · rw [h1] at h
    have h00 := congrFun (congrFun h 0) 0
    have hcol : colOf X 0 = ![2, -1, -1] := rfl
    have hm : mean (![2, -1, -1] : Fin 3 → ℝ) = 0 := by norm_num [mean, Fin.sum_univ_three]
    simp only [PMmatrix, Bool.false_eq_true, ↓reduceIte, clpmM, cupmM, Matrix.of_apply,
      colMeans, hcol, hm, pmScale] at h00
    norm_num [CoLPM, CoUPM, lowerDev, upperDev, Fin.sum_univ_three] at h00

end Matrix

end GapProjection

end

end NNS
