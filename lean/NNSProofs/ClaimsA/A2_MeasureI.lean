-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.GapDensity

namespace NNS.ClaimsA

open MeasureTheory ProbabilityTheory Set GapQuantile GapDensity Filter Topology
open scoped ENNReal

@[expose] public noncomputable section

variable {Ω : Type*} {μ : Measure ℝ}

private theorem parts_mul (x t : ℝ) : max (x - t) 0 * max (t - x) 0 = 0 := by
  rcases le_total x t with h | h
  · simp [max_eq_right (sub_nonpos.mpr h)]
  · simp [max_eq_right (sub_nonpos.mpr h)]

theorem claim_C003 (x t : ℝ) :
    |x - t| = max (x - t) 0 + max (t - x) 0 ∧ max (x - t) 0 * max (t - x) 0 = 0 :=
  ⟨by simpa only [neg_sub] using (max_add_max_neg (x - t)).symm, parts_mul x t⟩

theorem claim_C004 (x : ℝ) : x = max x 0 - max (-x) 0 := (max_sub_max_neg x).symm

theorem claim_C005 (x t : ℝ) : x - t = max (x - t) 0 - max (t - x) 0 := by
  simpa only [neg_sub] using claim_C004 (x - t)

theorem claim_C006 (x t : ℝ) : |x - t| = max (x - t) 0 + max (t - x) 0 :=
  (claim_C003 x t).1

theorem claim_C007 (x t : ℝ) :
    x - t = max (x - t) 0 - max (t - x) 0 ∧ |x - t| = max (x - t) 0 + max (t - x) 0 :=
  ⟨claim_C005 x t, claim_C006 x t⟩

theorem claim_C008 (X : Ω → ℝ) (t : ℝ) :
    {ω | X ω ≤ t} ∪ {ω | t < X ω} = univ ∧ Disjoint {ω | X ω ≤ t} {ω | t < X ω} ∧
    (∀ ω, X ω ≤ t → max (t - X ω) 0 = t - X ω ∧ max (X ω - t) 0 = 0) ∧
    (∀ ω, t < X ω → max (X ω - t) 0 = X ω - t ∧ max (t - X ω) 0 = 0) := by
  refine ⟨eq_univ_of_forall (fun ω => le_or_gt (X ω) t),
    ?_, ?_, ?_⟩
  · simp [Set.disjoint_left]
  · exact fun _ h => ⟨max_eq_left (sub_nonneg.mpr h), max_eq_right (sub_nonpos.mpr h)⟩
  · exact fun _ h => ⟨max_eq_left (sub_nonneg.mpr h.le), max_eq_right (sub_nonpos.mpr h.le)⟩

theorem claim_C009 (r : ℕ) (t : ℝ) (μ : Measure ℝ) :
    eLPM r t μ = ∫⁻ x, ENNReal.ofReal
      (if r = 0 then (if x ≤ t then 1 else 0) else max (t - x) 0 ^ r) ∂μ := by
  cases r with
  | zero => simp [eLPM, lowerDev]
  | succ k => simp [eLPM, lowerDev_eq_max_pow]

theorem claim_C010 (r : ℕ) (t : ℝ) (μ : Measure ℝ) :
    eUPM r t μ = ∫⁻ x, ENNReal.ofReal
      (if r = 0 then (if t < x then 1 else 0) else max (x - t) 0 ^ r) ∂μ := by
  cases r with
  | zero => simp [eUPM, upperDev]
  | succ k => simp [eUPM, upperDev_eq_max_pow]

theorem claim_C011 [IsProbabilityMeasure μ] (r : ℕ) (t : ℝ) :
    (Integrable (lowerDev r t) μ ↔ eLPM r t μ < ∞) ∧
    (Integrable (upperDev r t) μ ↔ eUPM r t μ < ∞) ∧
    (MemLp (fun x : ℝ => x) r μ → eLPM r t μ < ∞ ∧ eUPM r t μ < ∞) :=
  ⟨(eLPM_lt_top_iff r t).symm, (eUPM_lt_top_iff r t).symm,
    fun h => ⟨eLPM_lt_top h t, eUPM_lt_top h t⟩⟩

theorem claim_C012 [IsProbabilityMeasure μ] (t : ℝ) :
    (∀ x, lowerDev 0 t x = (if x ≤ t then 1 else 0) ∧
      upperDev 0 t x = (if t < x then 1 else 0)) ∧ mLPM 0 t μ + mUPM 0 t μ = 1 ∧
    (∀ x : ℝ, max (t - x) 0 ^ 0 = 1 ∧ max (x - t) 0 ^ 0 = 1) :=
  ⟨fun x => by simp [lowerDev, upperDev], mLPM_zero_add_mUPM_zero t, fun _ => by simp⟩

theorem claim_C015 [IsProbabilityMeasure μ] (t : ℝ) :
    (eLPM 0 t μ = μ (Iic t) ∧ eUPM 0 t μ = μ (Ioi t)) ∧
    (∀ r : ℕ, 0 < r →
      eLPM r t μ = (∫⁻ x, ENNReal.ofReal (max (t - x) 0 ^ r) ∂μ) ∧
      eUPM r t μ = (∫⁻ x, ENNReal.ofReal (max (x - t) 0 ^ r) ∂μ)) ∧
    (MemLp (fun x : ℝ => x) 2 μ →
      ENNReal.ofReal (variance (fun x : ℝ => x) μ) =
        eUPM 2 (∫ x, x ∂μ) μ + eLPM 2 (∫ x, x ∂μ) μ) := by
  refine ⟨⟨?_, ?_⟩, fun r hr => ⟨?_, ?_⟩, fun h => ?_⟩
  · simpa [eLPM, lowerDev, indicator, apply_ite ENNReal.ofReal] using
      (lintegral_indicator_one (μ := μ) (s := Iic t) measurableSet_Iic)
  · simpa [eUPM, upperDev, indicator, apply_ite ENNReal.ofReal] using
      (lintegral_indicator_one (μ := μ) (s := Ioi t) measurableSet_Ioi)
  · simpa [hr.ne'] using claim_C009 r t μ
  · simpa [hr.ne'] using claim_C010 r t μ
  · have hu := (eUPM_lt_top h (∫ x, x ∂μ)).ne
    have hl := (eLPM_lt_top h (∫ x, x ∂μ)).ne
    rw [variance_eq_mUPM_two_add_mLPM_two h, mUPM_eq_toReal, mLPM_eq_toReal,
      ← ENNReal.toReal_add hu hl, ENNReal.ofReal_toReal (ENNReal.add_ne_top.mpr ⟨hu, hl⟩)]

theorem claim_C017 {n : ℕ} (r : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    mLPM r t (empMeasure x) = LPM r t x ∧ mUPM r t (empMeasure x) = UPM r t x :=
  ⟨mLPM_empMeasure r t x, mUPM_empMeasure r t x⟩

theorem claim_C021 [IsProbabilityMeasure μ] (k : ℕ) (t : ℝ)
    (h : MemLp (fun x : ℝ => x) (k + 1 : ℕ) μ) :
    ∫ x, (x - t) ^ (k + 1) ∂μ = mUPM (k + 1) t μ + (-1) ^ (k + 1) * mLPM (k + 1) t μ :=
  (mUPM_add_neg_one_pow_mul_mLPM h t).symm

theorem claim_C022 [IsProbabilityMeasure μ] (h : Integrable (fun x : ℝ => x) μ) :
    ∫ x, x ∂μ = mUPM 1 0 μ - mLPM 1 0 μ := mean_eq_mUPM_one_sub_mLPM_one h

theorem claim_C023 [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ) :
    variance (fun x : ℝ => x) μ = mUPM 2 (∫ x, x ∂μ) μ + mLPM 2 (∫ x, x ∂μ) μ :=
  variance_eq_mUPM_two_add_mLPM_two h

theorem claim_C029 [IsProbabilityMeasure μ] (t : ℝ) : cdf μ t = mLPM 0 t μ :=
  (mLPM_zero_eq_cdf t).symm

theorem claim_C030 (r : ℕ) (t : ℝ) (μ : Measure ℝ) :
    eLPM r t μ = (∫⁻ x, ENNReal.ofReal
      (if r = 0 then (if x ≤ t then 1 else 0) else max (t - x) 0 ^ r) ∂μ) ∧
    eUPM r t μ = (∫⁻ x, ENNReal.ofReal
      (if r = 0 then (if t < x then 1 else 0) else max (x - t) 0 ^ r) ∂μ) :=
  ⟨claim_C009 r t μ, claim_C010 r t μ⟩

theorem claim_C031 (t x : ℝ) :
    lowerDev 0 t x = (if x ≤ t then 1 else 0) ∧ upperDev 0 t x = (if t < x then 1 else 0) := by
  simp [lowerDev, upperDev]

theorem claim_C032 (t : ℝ) : lowerDev 0 t = (Iic t).indicator (fun _ => (1 : ℝ)) := by
  simpa only [pow_zero] using lowerDev_eq_indicator 0 t

theorem claim_C033 [IsProbabilityMeasure μ] (t : ℝ) :
    mLPM 0 t μ = μ.real (Iic t) ∧ μ.real (Iic t) = cdf μ t :=
  ⟨(mLPM_zero_eq_cdf t).trans (cdf_eq_real _ _), (cdf_eq_real _ _).symm⟩

theorem claim_C040 [IsProbabilityMeasure μ] (t : ℝ) :
    mLPM 0 t μ = μ.real (Iic t) ∧ mUPM 0 t μ = μ.real (Ioi t) :=
  ⟨(claim_C033 t).1, mUPM_zero_eq t⟩

theorem claim_C051 [IsProbabilityMeasure μ] (p : ℝ) :
    quantile μ p = sInf {t | p ≤ cdf μ t} ∧
    sInf {t | p ≤ cdf μ t} = sInf {t | p ≤ mLPM 0 t μ} ∧
    (p ∈ Ioo 0 1 → {t | p ≤ mLPM 0 t μ}.Nonempty ∧ BddBelow {t | p ≤ mLPM 0 t μ}) := by
  refine ⟨rfl, ?_, ?_⟩
  · simp only [mLPM_zero_eq_cdf]
  · intro hp
    simpa only [mLPM_zero_eq_cdf] using
      And.intro (quantile_set_nonempty μ hp.2) (quantile_set_bddBelow μ hp.1)

theorem claim_C056 [IsProbabilityMeasure μ] {p : ℝ} (hp : p ∈ Ioo 0 1) :
    quantile μ p = sInf {t | p ≤ mLPM 0 t μ} ∧ IsLeast {t | p ≤ mLPM 0 t μ} (quantile μ p) := by
  refine ⟨(claim_C051 p).1.trans (claim_C051 p).2.1, ?_⟩
  simp only [mLPM_zero_eq_cdf, setOf_le_cdf_eq_Ici μ hp.1 hp.2]
  exact isLeast_Ici

theorem claim_C060 [IsProbabilityMeasure μ] (t : ℝ) {p : ℝ} (hp : p ∈ Ioo 0 1) :
    cdf μ t = mLPM 0 t μ ∧ 1 - cdf μ t = mUPM 0 t μ ∧
    (μ.real (Iio (quantile μ p)) ≤ p ∧ p ≤ μ.real (Iic (quantile μ p))) ∧
    (quantile μ p ≤ t ↔ p ≤ mLPM 0 t μ) ∧
    (ContinuousAt (cdf μ) (quantile μ p) → mLPM 0 (quantile μ p) μ = p) ∧
    cdf (empMeasure ![(0 : ℝ), 0, 1]) (quantile (empMeasure ![(0 : ℝ), 0, 1]) (1/2)) = 2/3 := by
  refine ⟨claim_C029 t, ?_, measureReal_Iio_quantile_le_and_le_Iic hp.1 hp.2, ?_, ?_,
    cdf_quantile_tied_sample⟩
  · rw [claim_C029]
    linarith [mLPM_zero_add_mUPM_zero (μ := μ) t]
  · rw [mLPM_zero_eq_cdf]
    exact quantile_le_iff μ hp.1 hp.2
  · intro hc
    rw [mLPM_zero_eq_cdf]
    refine le_antisymm ?_ (le_cdf_quantile μ hp.2)
    have hl : Tendsto (cdf μ) (𝓝[<] (quantile μ p)) (𝓝 (cdf μ (quantile μ p))) :=
      hc.continuousWithinAt.tendsto
    apply le_of_tendsto hl
    filter_upwards [self_mem_nhdsWithin] with s hs
    exact (cdf_lt_of_lt_quantile μ hp.1 hs).le

theorem claim_C061 [IsProbabilityMeasure μ] (t : ℝ) : cdf μ t = mLPM 0 t μ := claim_C029 t

theorem claim_C105 [IsProbabilityMeasure μ] (h : Integrable (fun x : ℝ => x) μ) (t : ℝ) :
    (∫ x, x - t ∂μ) = (∫ x, x ∂μ) - t ∧ (∫ x, x - (∫ y, y ∂μ) ∂μ) = 0 := by
  simp [integral_sub h (integrable_const _), integral_const]

theorem claim_C109 (x t : ℝ) : x - t = max (x - t) 0 - max (t - x) 0 := claim_C005 x t

theorem claim_C110 [IsProbabilityMeasure μ] (r : ℕ) (t : ℝ) (h : MemLp (fun x : ℝ => x) r μ) :
    ∫ x, (x - t) ^ r ∂μ = mUPM r t μ + (-1) ^ r * mLPM r t μ :=
  (mUPM_add_neg_one_pow_mul_mLPM h t).symm

theorem claim_C111 [IsProbabilityMeasure μ] (r : ℕ) (t : ℝ) (h : MemLp (fun x : ℝ => x) r μ) :
    ∫ x, (x - t) ^ r ∂μ = mUPM r t μ + (-1) ^ r * mLPM r t μ := claim_C110 r t h

theorem claim_C120 [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 4 μ) :
    ∫ x, (x - ∫ y, y ∂μ) ^ 4 ∂μ = mUPM 4 (∫ x, x ∂μ) μ + mLPM 4 (∫ x, x ∂μ) μ := by
  have h4 : MemLp (fun x : ℝ => x) ((4 : ℕ) : ℝ≥0∞) μ := by exact_mod_cast h
  convert! claim_C110 4 (∫ x, x ∂μ) h4 using 1
  norm_num

private theorem setIntegral_density {f : ℝ → ℝ} (hf : Measurable f) (hf0 : ∀ x, 0 ≤ f x)
    {s : Set ℝ} (hs : MeasurableSet s) (g : ℝ → ℝ) :
    ∫ x in s, g x ∂densityLaw f = ∫ x in s, g x * f x := by
  rw [← integral_indicator hs, ← integral_indicator hs, densityLaw,
    integral_withDensity_ofReal hf hf0]
  simp only [indicator_mul_left]

theorem claim_C130 {f : ℝ → ℝ} (hf : Measurable f) (hf0 : ∀ x, 0 ≤ f x)
    [IsProbabilityMeasure (densityLaw f)] (r : ℕ) (t : ℝ)
    (h : MemLp (fun x : ℝ => x) r (densityLaw f)) :
    (∫ x, (x - t) ^ r ∂densityLaw f) =
      (∫ x in Ioi t, (x - t) ^ r * f x) + (∫ x in Iic t, (x - t) ^ r * f x) ∧
    (∫ x in Ioi t, (x - t) ^ r * f x) = mUPM r t (densityLaw f) ∧
    (∫ x in Iic t, (x - t) ^ r * f x) = (-1) ^ r * mLPM r t (densityLaw f) := by
  have hu := (mUPM_eq_setIntegral r t (densityLaw f)).trans
    (setIntegral_density hf hf0 measurableSet_Ioi _)
  have hl : (-1) ^ r * mLPM r t (densityLaw f) = ∫ x in Iic t, (x - t) ^ r * f x := by
    rw [mLPM_eq_setIntegral, setIntegral_density hf hf0 measurableSet_Iic, ← integral_const_mul]
    apply integral_congr_ae
    filter_upwards with x
    rw [← mul_assoc, ← mul_pow]
    ring_nf
  exact ⟨by rw [claim_C110 r t h, hu, hl], hu.symm, hl.symm⟩

theorem claim_C131 {f : ℝ → ℝ} (hf : Measurable f) (hf0 : ∀ x, 0 ≤ f x)
    [IsProbabilityMeasure (densityLaw f)] :
    let m := ∫ x, x ∂densityLaw f
    mUPM 2 m (densityLaw f) = ∫ x in Ioi m, (x - m) ^ 2 * f x := by
  exact (mUPM_eq_setIntegral 2 _ _).trans (setIntegral_density hf hf0 measurableSet_Ioi _)

theorem claim_C136 [MeasurableSpace Ω] {X : Ω → ℝ} (hX : Measurable X) (t : ℝ) :
    MeasurableSet {ω | X ω ≤ t} ∧ MeasurableSet {ω | t < X ω} ∧
    Disjoint {ω | X ω ≤ t} {ω | t < X ω} ∧ {ω | X ω ≤ t} ∪ {ω | t < X ω} = univ :=
  ⟨hX measurableSet_Iic, hX measurableSet_Ioi, (claim_C008 X t).2.1, (claim_C008 X t).1⟩

theorem claim_C137 (f : Ω → ℝ) : f = fun ω => max (f ω) 0 - max (-f ω) 0 :=
  funext fun ω => claim_C004 (f ω)

theorem claim_C138 [MeasurableSpace Ω] {X : Ω → ℝ} (hX : Measurable X) (t : ℝ) :
    (∀ ω, 0 ≤ max (X ω - t) 0 ∧ 0 ≤ max (t - X ω) 0) ∧
    Measurable (fun ω => max (X ω - t) 0) ∧ Measurable (fun ω => max (t - X ω) 0) :=
  ⟨fun _ => ⟨le_max_right _ _, le_max_right _ _⟩,
    (hX.sub measurable_const).max measurable_const,
    (measurable_const.sub hX).max measurable_const⟩

theorem claim_C140 (X : Ω → ℝ) (t : ℝ) :
    (fun ω => X ω - t) = fun ω => max (X ω - t) 0 - max (t - X ω) 0 :=
  funext fun ω => claim_C005 (X ω) t

theorem claim_C141 (k : ℕ) (x t : ℝ) :
    max (x - t) 0 * max (t - x) 0 = 0 ∧
    (max (x - t) 0 - max (t - x) 0) ^ (k + 1) =
      max (x - t) 0 ^ (k + 1) + (-1) ^ (k + 1) * max (t - x) 0 ^ (k + 1) := by
  refine ⟨parts_mul x t, ?_⟩
  rw [← claim_C005]
  simpa only [upperDev_eq_max_pow (Nat.succ_ne_zero k), lowerDev_eq_max_pow (Nat.succ_ne_zero k)]
    using (upperDev_add_neg_one_pow_mul_lowerDev (k + 1) t x).symm

end

end NNS.ClaimsA
