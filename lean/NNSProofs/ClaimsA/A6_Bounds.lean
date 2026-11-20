-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.GapOperators
public import NNSProofs.GapQuantile

namespace NNS.ClaimsA

open MeasureTheory ProbabilityTheory Set

@[expose] public noncomputable section

variable {μ : Measure ℝ} [IsProbabilityMeasure μ] {g t α : ℝ}
variable {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P] {X : Ω → ℝ}

theorem claim_C089 (hX : MemLp X 2 P) (hg : g < ∫ ω, X ω ∂P) :
    P.real {ω | X ω ≤ g} ≤ variance X P / (variance X P + ((∫ ω, X ω ∂P) - g) ^ 2) :=
  rv_cantelli hX hg

theorem claim_C090 (hμ : MemLp id 2 μ) :
    (μ.map (fun v => 2 * pmean μ - v) = μ → ∀ g < pmean μ,
      μ.real (Iic g) ≤ pvar μ / (2 * (pmean μ - g) ^ 2)) ∧
    (∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧ MemLp id 2 ν ∧
      ν.map (fun v => 2 * pmean ν - v) ≠ ν ∧
      ∀ g < pmean ν, ν.real (Iic g) ≤ pvar ν / (2 * (pmean ν - g) ^ 2)) ∧
    pvar (empirical zSample) / (2 * (pmean (empirical zSample) - (-1)) ^ 2) <
      (empirical zSample).real (Iic (-1)) := by
  refine ⟨fun hs _ hg => halved_chebyshev_of_symm hμ hg hs, ?_, ?_⟩
  · let x : Fin 3 → ℝ := ![-1, -1, 2]
    let ν := empirical x
    have : IsProbabilityMeasure ν := isProbabilityMeasure_empirical (by decide) x
    have hm : mean x = 0 := by norm_num [mean, x, Fin.sum_univ_three]
    have hv : pvar ν = 2 := by
      change pvar (empirical x) = 2
      rw [pvar_empirical, popVariance, centralMoment, hm]
      norm_num [x, Fin.sum_univ_three]
    have hs : semivar ν = 2 / 3 := by
      change semivar (empirical x) = 2 / 3
      rw [semivar_empirical, hm]
      norm_num [LPM, lowerDev, x, Fin.sum_univ_three]
    refine ⟨ν, inferInstance, memLp_empirical x, ?_, fun _ hg => ?_⟩
    · intro hsym
      have heq := semivar_eq_half_pvar_of_symm (memLp_empirical x) hsym
      change semivar ν = pvar ν / 2 at heq
      linarith
    · exact halved_chebyshev_of_semivar_le (memLp_empirical x) hg (by rw [hs, hv]; norm_num)
  · rw [zSample_bounds.1, zSample_bounds.2.1]
    norm_num

theorem claim_C091 :
    let ν := empirical zSample
    pmean ν = 0 ∧ pvar ν = 1 / 3 ∧ ν.real (Iic (-1)) = 1 / 4 ∧
      pvar ν / (2 * (pmean ν - (-1)) ^ 2) = 1 / 6 ∧
      pvar ν / (pvar ν + (pmean ν - (-1)) ^ 2) = ν.real (Iic (-1)) ∧
      pvar ν / (2 * (pmean ν - (-1)) ^ 2) < ν.real (Iic (-1)) := by
  dsimp only
  rw [pmean_empirical, pvar_empirical, zSample_facts.1, zSample_facts.2.1, zSample_bounds.1]
  norm_num

theorem claim_C092 (_hmean : Integrable id μ) (hg : g < pmean μ) :
    μ (Iic g) ≤ (∫⁻ v, ENNReal.ofReal (max (pmean μ - v) 0 ^ 2) ∂μ) /
      ENNReal.ofReal ((pmean μ - g) ^ 2) ∧
    (Integrable (lowerDev 2 (pmean μ)) μ →
      μ.real (Iic g) ≤ (Real.sqrt (∫ v, max (pmean μ - v) 0 ^ 2 ∂μ) / (pmean μ - g)) ^ 2) := by
  constructor
  · simpa only [rLowerDev_eq_max (by norm_num : (2 : ℝ) ≠ 0), Real.rpow_two] using
      measure_Iic_le_lintegral_div μ (by norm_num : (0 : ℝ) ≤ 2) hg
  · intro h
    have hl : lpm 2 (pmean μ) μ = ∫ v, max (pmean μ - v) 0 ^ 2 ∂μ := by
      simp only [lpm, lowerDev_eq_max_pow (by decide : 2 ≠ 0)]
    rw [← hl, div_pow, Real.sq_sqrt (lpm_nonneg _ _ _)]
    exact real_Iic_le_lpm_div 2 hg h

theorem claim_C093 (hg : g < pmean μ) :
    (Integrable (lowerDev 2 (pmean μ)) μ →
      μ.real (Iic g) ≤ semivar μ / (pmean μ - g) ^ 2) ∧
    (MemLp id 2 μ → semivar μ / (pmean μ - g) ^ 2 ≤ pvar μ / (pmean μ - g) ^ 2) :=
  ⟨fun h => real_Iic_le_lpm_div 2 hg h, fun h => semivar_div_le_pvar_div h⟩

theorem claim_C094 (hμ : MemLp id 2 μ) (hs : μ.map (fun v => 2 * pmean μ - v) = μ) :
    semivar μ = pvar μ / 2 ∧ ∀ g < pmean μ,
      (Real.sqrt (semivar μ) / (pmean μ - g)) ^ 2 = pvar μ / (2 * (pmean μ - g) ^ 2) ∧
      μ.real (Iic g) ≤ pvar μ / (2 * (pmean μ - g) ^ 2) := by
  have heq := semivar_eq_half_pvar_of_symm hμ hs
  refine ⟨heq, fun b hb => ⟨?_, halved_chebyshev_of_symm hμ hb hs⟩⟩
  rw [div_pow, Real.sq_sqrt (semivar_nonneg μ), heq, div_div, mul_comm]

theorem claim_C096 (hα : 0 < α) (t : ℝ) :
    theta t α μ = (∫ v, max (t - v) 0 ^ α ∂μ) ^ (1 / α) ∧
    (∀ g < t, Integrable (rLowerDev α t) μ →
      μ.real (Iic g) ≤ (theta t α μ / (t - g)) ^ α) ∧
    ((Measure.dirac t).real (Iic t) = 1 ∧ theta t α (Measure.dirac t) = 0 ∧ t - t = 0 ∧
      ¬ (Measure.dirac t).real (Iic t) ≤ (theta t α (Measure.dirac t) / (t - t)) ^ α) := by
  refine ⟨by simp only [theta, rlpm, rLowerDev_eq_max hα.ne'],
    fun _ hg hint => real_Iic_le_theta hα hg hint, ?_⟩
  simp [theta, rlpm_dirac_self hα, Real.zero_rpow hα.ne', Real.zero_rpow (inv_ne_zero hα.ne')]

theorem claim_C098 (μ : Measure ℝ) [IsProbabilityMeasure μ] (g : ℝ) :
    (MemLp id 2 μ → semivar μ + upm 2 (pmean μ) μ = pvar μ) ∧
    (g < pmean μ → Integrable (lowerDev 2 (pmean μ)) μ →
      μ.real (Iic g) ≤ semivar μ / (pmean μ - g) ^ 2) ∧
    (∀ (α t : ℝ), 0 < α → g < t → Integrable (rLowerDev α t) μ →
      μ.real (Iic g) ≤ (theta t α μ / (t - g)) ^ α) ∧
    (theta (pmean μ) 2 μ / (pmean μ - g)) ^ (2 : ℝ) = semivar μ / (pmean μ - g) ^ 2 :=
  ⟨fun h => semivar_add_upm_two h, fun hg => (claim_C093 hg).1,
    fun _ _ ha hg h => real_Iic_le_theta ha hg h, theta_mean_two_div μ g⟩

theorem claim_C099 :
    let ν := empirical wSample
    let gap := pmean ν - (-50)
    ν.real (Iic (-50)) < pvar ν / (pvar ν + gap ^ 2) ∧
      pvar ν / (pvar ν + gap ^ 2) < semivar ν / gap ^ 2 := by
  dsimp only
  rw [wSample_bounds.1, wSample_bounds.2.1, wSample_bounds.2.2.1]
  norm_num

theorem claim_C108 (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ) :
    (∀ r : ℝ, 0 < r →
      rlpm r t μ = ∫ v, max (t - v) 0 ^ r ∂μ ∧
      GapOperators.rupm r t μ = ∫ v, max (v - t) 0 ^ r ∂μ) ∧
    (∀ v, rLowerDev 0 t v = (if v ≤ t then 1 else 0) ∧
      rUpperDev 0 t v = (if t < v then 1 else 0)) ∧
    rlpm 0 t μ = μ.real (Iic t) ∧ GapOperators.rupm 0 t μ = μ.real (Ioi t) :=
  ⟨fun _ hr => ⟨GapOperators.rlpm_eq_integral_max_rpow hr t μ,
      GapOperators.rupm_eq_integral_max_rpow hr.ne' t μ⟩,
    fun v => ⟨rLowerDev_zero t v, rUpperDev_zero t v⟩,
    GapOperators.rlpm_zero t μ, GapOperators.rupm_zero t μ⟩

theorem claim_C115 (hμ : MemLp id 2 μ) :
    variance id μ = ∫ v, (v - pmean μ) ^ 2 ∂μ ∧
    variance id μ = upm 2 (pmean μ) μ + lpm 2 (pmean μ) μ := by
  rw [← pvar_eq_variance]
  exact ⟨rfl, (semivar_add_upm_two hμ).symm.trans (add_comm _ _)⟩

theorem claim_C354 (t : ℝ) :
    lpm 0 t μ = μ.real (Iic t) ∧ lpm 0 t μ = cdf μ t :=
  ⟨lpm_zero t μ, lpm_zero_eq_cdf t μ⟩

theorem claim_C361 (t : ℝ) :
    lpm 0 t μ = μ.real (Iic t) ∧ upm 0 t μ = μ.real (Ioi t) ∧
    lpm 0 t μ + upm 0 t μ = 1 ∧ mPMRatio 0 t μ = cdf μ t :=
  ⟨lpm_zero t μ, mUPM_zero_eq t, mLPM_zero_add_mUPM_zero t, mPMRatio_zero_eq_cdf t⟩

theorem claim_C364 (hX : MemLp X 2 P) (hg : g < ∫ ω, X ω ∂P) :
    P.real {ω | X ω ≤ g} ≤ variance X P / (variance X P + ((∫ ω, X ω ∂P) - g) ^ 2) :=
  claim_C089 hX hg

theorem claim_C365 (hμ : MemLp id 2 μ) :
    (μ.map (fun v => 2 * pmean μ - v) = μ → ∀ g < pmean μ,
      μ.real (Iic g) ≤ pvar μ / (2 * (pmean μ - g) ^ 2)) ∧
    (∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧ MemLp id 2 ν ∧
      ν.map (fun v => 2 * pmean ν - v) ≠ ν ∧
      ∀ g < pmean ν, ν.real (Iic g) ≤ pvar ν / (2 * (pmean ν - g) ^ 2)) ∧
    pvar (empirical zSample) / (2 * (pmean (empirical zSample) - (-1)) ^ 2) <
      (empirical zSample).real (Iic (-1)) :=
  claim_C090 hμ

theorem claim_C366 (hmean : Integrable id μ) (hg : g < pmean μ) :
    μ (Iic g) ≤ (∫⁻ v, ENNReal.ofReal (max (pmean μ - v) 0 ^ 2) ∂μ) /
      ENNReal.ofReal ((pmean μ - g) ^ 2) ∧
    (Integrable (lowerDev 2 (pmean μ)) μ →
      μ.real (Iic g) ≤ (Real.sqrt (∫ v, max (pmean μ - v) 0 ^ 2 ∂μ) / (pmean μ - g)) ^ 2) :=
  claim_C092 hmean hg

theorem claim_C367 (μ : Measure ℝ) (hα : 0 < α) (t : ℝ) :
    theta t α μ = (∫ v, max (t - v) 0 ^ α ∂μ) ^ (1 / α) := by
  simp only [theta, rlpm, rLowerDev_eq_max hα.ne']

theorem claim_C368 (hα : 0 < α) (t : ℝ) :
    theta t α μ = (∫ v, max (t - v) 0 ^ α ∂μ) ^ (1 / α) ∧
    (∀ g < t, Integrable (rLowerDev α t) μ →
      μ.real (Iic g) ≤ (theta t α μ / (t - g)) ^ α) ∧
    ((Measure.dirac t).real (Iic t) = 1 ∧ theta t α (Measure.dirac t) = 0 ∧ t - t = 0 ∧
      ¬ (Measure.dirac t).real (Iic t) ≤ (theta t α (Measure.dirac t) / (t - t)) ^ α) :=
  claim_C096 hα t

theorem claim_C369 (μ : Measure ℝ) [IsProbabilityMeasure μ] (g : ℝ) :
    (g < pmean μ → MemLp id 2 μ →
      μ.real (Iic g) ≤ pvar μ / (pvar μ + (pmean μ - g) ^ 2)) ∧
    (g < pmean μ → Integrable (lowerDev 2 (pmean μ)) μ →
      μ.real (Iic g) ≤ semivar μ / (pmean μ - g) ^ 2) ∧
    (∀ (α t : ℝ), 0 < α → g < t → Integrable (rLowerDev α t) μ →
      μ.real (Iic g) ≤ (theta t α μ / (t - g)) ^ α) ∧
    (theta (pmean μ) 2 μ / (pmean μ - g)) ^ (2 : ℝ) = semivar μ / (pmean μ - g) ^ 2 :=
  ⟨fun hg h => cantelli h hg, (claim_C098 μ g).2⟩

theorem claim_C375 :
    (∀ t, lpm 0 t μ = μ.real (Iic t) ∧ lpm 0 t μ = cdf μ t) ∧
    (∀ p ∈ Ioo (0 : ℝ) 1, sInf {t | p ≤ lpm 0 t μ} = GapQuantile.quantile μ p ∧
      μ.real (Iio (GapQuantile.quantile μ p)) ≤ p ∧ p ≤ μ.real (Iic (GapQuantile.quantile μ p))) ∧
    (∀ g < pmean μ, Integrable (lowerDev 2 (pmean μ)) μ →
      μ.real (Iic g) ≤ semivar μ / (pmean μ - g) ^ 2) ∧
    ∀ (α t g : ℝ), 0 < α → g < t →
      μ (Iic g) ≤ (∫⁻ v, ENNReal.ofReal (max (t - v) 0 ^ α) ∂μ) / ENNReal.ofReal ((t - g) ^ α) := by
  refine ⟨fun b => claim_C354 b, fun p hp => ⟨?_, ?_⟩, fun _ hg h => real_Iic_le_lpm_div 2 hg h,
    fun a b c ha hc => ?_⟩
  · simp only [lpm_zero_eq_cdf, GapQuantile.quantile_eq_sInf]
  · exact GapQuantile.measureReal_Iio_quantile_le_and_le_Iic hp.1 hp.2
  · simpa only [rLowerDev_eq_max ha.ne'] using measure_Iic_le_lintegral_div μ ha.le hc

theorem claim_C611 :
    (∀ g < pmean μ, MemLp id 2 μ →
      μ.real (Iic g) ≤ pvar μ / (pvar μ + (pmean μ - g) ^ 2)) ∧
    (∀ (α t g : ℝ), 0 < α → g < t → Integrable (rLowerDev α t) μ →
      μ.real (Iic g) ≤ rlpm α t μ / (t - g) ^ α) ∧
    (MemLp id 2 μ → μ.map (fun v => 2 * pmean μ - v) = μ → ∀ g < pmean μ,
      μ.real (Iic g) ≤ pvar μ / (2 * (pmean μ - g) ^ 2)) ∧
    pvar (empirical zSample) / (2 * (pmean (empirical zSample) - (-1)) ^ 2) <
      (empirical zSample).real (Iic (-1)) :=
  ⟨fun _ hg h => cantelli h hg, fun _ _ _ ha hg h => real_Iic_le_rlpm_div ha.le hg h,
    fun h hs _ hg => halved_chebyshev_of_symm h hg hs, by
      rw [zSample_bounds.2.1, zSample_bounds.1]
      norm_num⟩

end

end NNS.ClaimsA
