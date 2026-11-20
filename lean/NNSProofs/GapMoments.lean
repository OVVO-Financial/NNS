-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.MeasurePM
public import NNSProofs.Moments

namespace NNS.GapMoments

open MeasureTheory ProbabilityTheory Finset
open scoped ENNReal

@[expose] public noncomputable section

variable {n : ℕ} {μ : Measure ℝ}

theorem mUPM_one_eq_mLPM_one_iff [IsProbabilityMeasure μ] (h : Integrable (fun x : ℝ => x) μ)
    (t : ℝ) : mUPM 1 t μ = mLPM 1 t μ ↔ t = ∫ x, x ∂μ := by
  have := mUPM_one_sub_mLPM_one h t
  constructor <;> intro h' <;> linarith

theorem mLPM_one_lt_mUPM_one_iff [IsProbabilityMeasure μ] (h : Integrable (fun x : ℝ => x) μ)
    (t : ℝ) : mLPM 1 t μ < mUPM 1 t μ ↔ t < ∫ x, x ∂μ := by
  have := mUPM_one_sub_mLPM_one h t
  constructor <;> intro h' <;> linarith

theorem mUPM_one_lt_mLPM_one_iff [IsProbabilityMeasure μ] (h : Integrable (fun x : ℝ => x) μ)
    (t : ℝ) : mUPM 1 t μ < mLPM 1 t μ ↔ ∫ x, x ∂μ < t := by
  have := mUPM_one_sub_mLPM_one h t
  constructor <;> intro h' <;> linarith

theorem UPM_one_eq_LPM_one_iff (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    UPM 1 t x = LPM 1 t x ↔ t = mean x := by
  have := UPM_one_sub_LPM_one hn t x
  constructor <;> intro h' <;> linarith

theorem mPMRatio_one_mean_of_not_ae_eq [IsProbabilityMeasure μ]
    (h : Integrable (fun x : ℝ => x) μ) (hμ : ¬ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ) :
    mPMRatio 1 (∫ y, y ∂μ) μ = 1 / 2 :=
  mPMRatio_one_mean h fun h0 => hμ (ae_iff.mpr h0)

theorem mPMRatio_one_mean_of_not_ae_const [IsProbabilityMeasure μ]
    (h : Integrable (fun x : ℝ => x) μ) (hμ : ¬ ∃ c : ℝ, ∀ᵐ x ∂μ, x = c) :
    mPMRatio 1 (∫ y, y ∂μ) μ = 1 / 2 :=
  mPMRatio_one_mean_of_not_ae_eq h fun h0 => hμ ⟨_, h0⟩

lemma mLPM_eq_zero_of_ae_eq {r : ℕ} (hr : r ≠ 0) {t : ℝ} (h : ∀ᵐ x ∂μ, x = t) :
    mLPM r t μ = 0 :=
  integral_eq_zero_of_ae (h.mono fun x hx => by simp [lowerDev, hx, zero_pow hr])

lemma mUPM_eq_zero_of_ae_eq (r : ℕ) {t : ℝ} (h : ∀ᵐ x ∂μ, x = t) : mUPM r t μ = 0 :=
  integral_eq_zero_of_ae (h.mono fun x hx => by simp [upperDev, hx])

lemma lowerDev_pos (r : ℕ) {t x : ℝ} (h : x < t) : 0 < lowerDev r t x := by
  simp only [lowerDev, h.le, ↓reduceIte]
  exact pow_pos (sub_pos.mpr h) r

lemma upperDev_pos (r : ℕ) {t x : ℝ} (h : t < x) : 0 < upperDev r t x := by
  simp only [upperDev, h, ↓reduceIte]
  exact pow_pos (sub_pos.mpr h) r

lemma variance_eq (μ : Measure ℝ) :
    variance (fun x : ℝ => x) μ = ∫ x, (x - ∫ y, y ∂μ) ^ 2 ∂μ :=
  variance_eq_integral aemeasurable_id'

lemma ae_le_of_mLPM_eq_zero {r : ℕ} {t : ℝ} (hi : Integrable (lowerDev r t) μ)
    (h0 : mLPM r t μ = 0) : ∀ᵐ x ∂μ, t ≤ x := by
  filter_upwards [(integral_eq_zero_iff_of_nonneg (lowerDev_nonneg r t) hi).mp h0] with x hx
  by_contra hlt
  have hlt' : x < t := not_le.mp hlt
  have hpos : 0 < lowerDev r t x := lowerDev_pos r hlt'
  exact hpos.ne' hx

lemma ae_le_of_mUPM_eq_zero {r : ℕ} {t : ℝ} (hi : Integrable (upperDev r t) μ)
    (h0 : mUPM r t μ = 0) : ∀ᵐ x ∂μ, x ≤ t := by
  filter_upwards [(integral_eq_zero_iff_of_nonneg (upperDev_nonneg r t) hi).mp h0] with x hx
  by_contra hlt
  have hlt' : t < x := not_le.mp hlt
  have hpos : 0 < upperDev r t x := upperDev_pos r hlt'
  exact hpos.ne' hx

lemma ae_le_of_eLPM_eq_zero {r : ℕ} {t : ℝ} (h0 : eLPM r t μ = 0) : ∀ᵐ x ∂μ, t ≤ x := by
  have hm : Measurable fun x => ENNReal.ofReal (lowerDev r t x) :=
    (lowerDev_measurable r t).ennreal_ofReal
  filter_upwards [(lintegral_eq_zero_iff hm).mp h0] with x hx
  by_contra hlt
  have hlt' : x < t := not_le.mp hlt
  have hpos : 0 < lowerDev r t x := lowerDev_pos r hlt'
  exact (ENNReal.ofReal_pos.mpr hpos).ne' hx

lemma ae_le_of_eUPM_eq_zero {r : ℕ} {t : ℝ} (h0 : eUPM r t μ = 0) : ∀ᵐ x ∂μ, x ≤ t := by
  have hm : Measurable fun x => ENNReal.ofReal (upperDev r t x) :=
    (upperDev_measurable r t).ennreal_ofReal
  filter_upwards [(lintegral_eq_zero_iff hm).mp h0] with x hx
  by_contra hlt
  have hlt' : t < x := not_le.mp hlt
  have hpos : 0 < upperDev r t x := upperDev_pos r hlt'
  exact (ENNReal.ofReal_pos.mpr hpos).ne' hx

theorem ae_eq_mean_of_ae_ge [IsProbabilityMeasure μ] (h : Integrable (fun x : ℝ => x) μ)
    (hge : ∀ᵐ x ∂μ, ∫ y, y ∂μ ≤ x) : ∀ᵐ x ∂μ, x = ∫ y, y ∂μ := by
  have hi : Integrable (fun x : ℝ => x - ∫ y, y ∂μ) μ := h.sub (integrable_const _)
  have h0 : ∫ x, (x - ∫ y, y ∂μ) ∂μ = 0 := by
    rw [integral_sub h (integrable_const _), integral_const, probReal_univ, one_smul, sub_self]
  have hnn : 0 ≤ᵐ[μ] fun x : ℝ => x - ∫ y, y ∂μ := hge.mono fun x hx => sub_nonneg.mpr hx
  filter_upwards [(integral_eq_zero_iff_of_nonneg_ae hnn hi).mp h0] with x hx
  exact sub_eq_zero.mp hx

theorem ae_eq_mean_of_ae_le [IsProbabilityMeasure μ] (h : Integrable (fun x : ℝ => x) μ)
    (hle : ∀ᵐ x ∂μ, x ≤ ∫ y, y ∂μ) : ∀ᵐ x ∂μ, x = ∫ y, y ∂μ := by
  have hi : Integrable (fun x : ℝ => (∫ y, y ∂μ) - x) μ := (integrable_const _).sub h
  have h0 : ∫ x, ((∫ y, y ∂μ) - x) ∂μ = 0 := by
    rw [integral_sub (integrable_const _) h, integral_const, probReal_univ, one_smul, sub_self]
  have hnn : 0 ≤ᵐ[μ] fun x : ℝ => (∫ y, y ∂μ) - x := hle.mono fun x hx => sub_nonneg.mpr hx
  filter_upwards [(integral_eq_zero_iff_of_nonneg_ae hnn hi).mp h0] with x hx
  exact (sub_eq_zero.mp hx).symm

theorem ae_eq_mean_of_mLPM_two_eq_zero [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 2 μ) (h0 : mLPM 2 (∫ y, y ∂μ) μ = 0) :
    (∀ᵐ x ∂μ, ∫ y, y ∂μ ≤ x) ∧ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ := by
  have h2 : MemLp (fun x : ℝ => x) ((2 : ℕ) : ℝ≥0∞) μ := by exact_mod_cast h
  have hge := ae_le_of_mLPM_eq_zero (integrable_lowerDev h2 _) h0
  exact ⟨hge, ae_eq_mean_of_ae_ge (h.integrable one_le_two) hge⟩

theorem ae_eq_mean_of_mUPM_two_eq_zero [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 2 μ) (h0 : mUPM 2 (∫ y, y ∂μ) μ = 0) :
    (∀ᵐ x ∂μ, x ≤ ∫ y, y ∂μ) ∧ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ := by
  have h2 : MemLp (fun x : ℝ => x) ((2 : ℕ) : ℝ≥0∞) μ := by exact_mod_cast h
  have hle := ae_le_of_mUPM_eq_zero (integrable_upperDev h2 _) h0
  exact ⟨hle, ae_eq_mean_of_ae_le (h.integrable one_le_two) hle⟩

theorem ae_eq_mean_of_eLPM_two_eq_zero [IsProbabilityMeasure μ]
    (h : Integrable (fun x : ℝ => x) μ) (h0 : eLPM 2 (∫ y, y ∂μ) μ = 0) :
    (∀ᵐ x ∂μ, ∫ y, y ∂μ ≤ x) ∧ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ :=
  ⟨ae_le_of_eLPM_eq_zero h0, ae_eq_mean_of_ae_ge h (ae_le_of_eLPM_eq_zero h0)⟩

theorem ae_eq_mean_of_eUPM_two_eq_zero [IsProbabilityMeasure μ]
    (h : Integrable (fun x : ℝ => x) μ) (h0 : eUPM 2 (∫ y, y ∂μ) μ = 0) :
    (∀ᵐ x ∂μ, x ≤ ∫ y, y ∂μ) ∧ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ :=
  ⟨ae_le_of_eUPM_eq_zero h0, ae_eq_mean_of_ae_le h (ae_le_of_eUPM_eq_zero h0)⟩

theorem mLPM_two_mean_eq_zero_iff [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ) :
    mLPM 2 (∫ y, y ∂μ) μ = 0 ↔ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ :=
  ⟨fun h0 => (ae_eq_mean_of_mLPM_two_eq_zero h h0).2, mLPM_eq_zero_of_ae_eq two_ne_zero⟩

theorem mUPM_two_mean_eq_zero_iff [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ) :
    mUPM 2 (∫ y, y ∂μ) μ = 0 ↔ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ :=
  ⟨fun h0 => (ae_eq_mean_of_mUPM_two_eq_zero h h0).2, mUPM_eq_zero_of_ae_eq 2⟩

theorem dirac_degenerate (c : ℝ) :
    ∫ x, x ∂Measure.dirac c = c ∧ mUPM 2 c (Measure.dirac c) = 0 ∧
      mLPM 2 c (Measure.dirac c) = 0 ∧ variance (fun x : ℝ => x) (Measure.dirac c) = 0 := by
  refine ⟨integral_dirac _ c, ?_, ?_, variance_dirac c⟩
  · simp [mUPM, upperDev]
  · simp [mLPM, lowerDev]

theorem degenerate_of_ae_eq [IsProbabilityMeasure μ] {c : ℝ} (h : ∀ᵐ x ∂μ, x = c) :
    ∫ x, x ∂μ = c ∧ mUPM 2 c μ = 0 ∧ mLPM 2 c μ = 0 ∧ variance (fun x : ℝ => x) μ = 0 := by
  have hm : ∫ y, y ∂μ = c := by
    rw [integral_congr_ae h, integral_const, probReal_univ, one_smul]
  refine ⟨hm, mUPM_eq_zero_of_ae_eq 2 h, mLPM_eq_zero_of_ae_eq two_ne_zero h, ?_⟩
  rw [variance_eq, hm]
  exact integral_eq_zero_of_ae (h.mono fun x hx => by simp [hx])

theorem mLPM_two_pos_and_mUPM_two_pos [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 2 μ) (hv : 0 < variance (fun x : ℝ => x) μ) :
    0 < mLPM 2 (∫ y, y ∂μ) μ ∧ 0 < mUPM 2 (∫ y, y ∂μ) μ := by
  have hvar := variance_eq_mUPM_two_add_mLPM_two h
  have hL := mLPM_nonneg 2 (∫ y, y ∂μ) μ
  have hU := mUPM_nonneg 2 (∫ y, y ∂μ) μ
  constructor
  · refine hL.lt_of_ne' fun h0 => ?_
    have := mUPM_eq_zero_of_ae_eq 2 ((mLPM_two_mean_eq_zero_iff h).mp h0)
    linarith
  · refine hU.lt_of_ne' fun h0 => ?_
    have := mLPM_eq_zero_of_ae_eq two_ne_zero ((mUPM_two_mean_eq_zero_iff h).mp h0)
    linarith

theorem variance_eq_zero_iff_ae_eq_mean [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 2 μ) :
    variance (fun x : ℝ => x) μ = 0 ↔ ∀ᵐ x ∂μ, x = ∫ y, y ∂μ := by
  refine ⟨fun h0 => ?_, fun h' => (degenerate_of_ae_eq h').2.2.2⟩
  have hvar := variance_eq_mUPM_two_add_mLPM_two h
  have hL := mLPM_nonneg 2 (∫ y, y ∂μ) μ
  have hU := mUPM_nonneg 2 (∫ y, y ∂μ) μ
  exact (mLPM_two_mean_eq_zero_iff h).mp (by linarith)

def mSkewness (μ : Measure ℝ) : ℝ :=
  (∫ x, (x - ∫ y, y ∂μ) ^ 3 ∂μ) / variance (fun x : ℝ => x) μ ^ ((3 : ℝ) / 2)

def mKurtosis (μ : Measure ℝ) : ℝ :=
  (∫ x, (x - ∫ y, y ∂μ) ^ 4 ∂μ) / variance (fun x : ℝ => x) μ ^ 2

theorem integral_id_empMeasure (x : Fin n → ℝ) : ∫ v, v ∂empMeasure x = mean x :=
  integral_empMeasure x _

theorem variance_empMeasure (x : Fin n → ℝ) :
    variance (fun v : ℝ => v) (empMeasure x) = popVariance x := by
  rw [variance_eq, integral_id_empMeasure, integral_empMeasure]
  rfl

theorem mSkewness_empMeasure (x : Fin n → ℝ) : mSkewness (empMeasure x) = popSkewness x := by
  rw [mSkewness, variance_empMeasure, integral_id_empMeasure, integral_empMeasure]
  rfl

theorem mKurtosis_empMeasure (x : Fin n → ℝ) :
    mKurtosis (empMeasure x) = popKurtosis x + 3 := by
  rw [mKurtosis, variance_empMeasure, integral_id_empMeasure, integral_empMeasure, popKurtosis,
    centralMoment]
  ring

theorem mSkewness_eq_PM [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 3 μ) :
    mSkewness μ = (mUPM 3 (∫ y, y ∂μ) μ - mLPM 3 (∫ y, y ∂μ) μ) /
      (mUPM 2 (∫ y, y ∂μ) μ + mLPM 2 (∫ y, y ∂μ) μ) ^ ((3 : ℝ) / 2) := by
  have h2 : MemLp (fun x : ℝ => x) 2 μ := h.mono_exponent (by norm_num)
  rw [mSkewness, thirdCentralMoment_eq_mUPM_three_sub_mLPM_three h,
    variance_eq_mUPM_two_add_mLPM_two h2]

lemma mUPM_three_eq_mLPM_three_of_variance_eq_zero [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 2 μ) (hv : variance (fun x : ℝ => x) μ = 0) :
    mUPM 3 (∫ y, y ∂μ) μ = 0 ∧ mLPM 3 (∫ y, y ∂μ) μ = 0 :=
  have hae := (variance_eq_zero_iff_ae_eq_mean h).mp hv
  ⟨mUPM_eq_zero_of_ae_eq 3 hae, mLPM_eq_zero_of_ae_eq three_ne_zero hae⟩

theorem mSkewness_pos_iff [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 3 μ) :
    0 < mSkewness μ ↔ mLPM 3 (∫ y, y ∂μ) μ < mUPM 3 (∫ y, y ∂μ) μ := by
  have h2 : MemLp (fun x : ℝ => x) 2 μ := h.mono_exponent (by norm_num)
  rcases (variance_nonneg (fun x : ℝ => x) μ).eq_or_lt with hv | hv
  · obtain ⟨hU, hL⟩ := mUPM_three_eq_mLPM_three_of_variance_eq_zero h2 hv.symm
    rw [mSkewness, ← hv, Real.zero_rpow (by norm_num), div_zero, hU, hL]
  · rw [mSkewness, thirdCentralMoment_eq_mUPM_three_sub_mLPM_three h,
      div_pos_iff_of_pos_right (Real.rpow_pos_of_pos hv _), sub_pos]

theorem mSkewness_neg_iff [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 3 μ) :
    mSkewness μ < 0 ↔ mUPM 3 (∫ y, y ∂μ) μ < mLPM 3 (∫ y, y ∂μ) μ := by
  have h2 : MemLp (fun x : ℝ => x) 2 μ := h.mono_exponent (by norm_num)
  rcases (variance_nonneg (fun x : ℝ => x) μ).eq_or_lt with hv | hv
  · obtain ⟨hU, hL⟩ := mUPM_three_eq_mLPM_three_of_variance_eq_zero h2 hv.symm
    rw [mSkewness, ← hv, Real.zero_rpow (by norm_num), div_zero, hU, hL]
  · rw [mSkewness, thirdCentralMoment_eq_mUPM_three_sub_mLPM_three h, ← neg_pos, ← neg_div,
      div_pos_iff_of_pos_right (Real.rpow_pos_of_pos hv _), neg_pos, sub_neg]

theorem fourthCentralMoment_eq_mUPM_four_add_mLPM_four [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 4 μ) :
    ∫ x, (x - ∫ y, y ∂μ) ^ 4 ∂μ = mUPM 4 (∫ y, y ∂μ) μ + mLPM 4 (∫ y, y ∂μ) μ := by
  have h4 : MemLp (fun x : ℝ => x) ((4 : ℕ) : ℝ≥0∞) μ := by exact_mod_cast h
  have := mUPM_add_neg_one_pow_mul_mLPM h4 (∫ y, y ∂μ)
  norm_num at this
  linarith

theorem mKurtosis_eq_PM [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 4 μ) :
    mKurtosis μ = (mUPM 4 (∫ y, y ∂μ) μ + mLPM 4 (∫ y, y ∂μ) μ) /
      (mUPM 2 (∫ y, y ∂μ) μ + mLPM 2 (∫ y, y ∂μ) μ) ^ 2 := by
  have h2 : MemLp (fun x : ℝ => x) 2 μ := h.mono_exponent (by norm_num)
  rw [mKurtosis, fourthCentralMoment_eq_mUPM_four_add_mLPM_four h,
    variance_eq_mUPM_two_add_mLPM_two h2]

theorem classical_moments_from_pairs_sample (hn : 0 < n) (x : Fin n → ℝ) :
    (∀ t : ℝ, UPM 1 t x = LPM 1 t x ↔ t = mean x) ∧
      mean x = UPM 1 0 x - LPM 1 0 x ∧
      popVariance x = LPM 2 (mean x) x + UPM 2 (mean x) x ∧
      popSkewness x = (UPM 3 (mean x) x - LPM 3 (mean x) x) /
        (LPM 2 (mean x) x + UPM 2 (mean x) x) ^ ((3 : ℝ) / 2) ∧
      popKurtosis x + 3 = (UPM 4 (mean x) x + LPM 4 (mean x) x) /
        (LPM 2 (mean x) x + UPM 2 (mean x) x) ^ 2 := by
  refine ⟨fun t => UPM_one_eq_LPM_one_iff hn t x, ?_, popVariance_eq_PM x, popSkewness_eq_PM x,
    ?_⟩
  · rw [UPM_one_sub_LPM_one hn, sub_zero]
  · rw [popKurtosis_eq_PM]
    ring

def condMean (μ : Measure ℝ) (A : Set ℝ) : ℝ := (∫ x in A, x ∂μ) / μ.real A

def condVar (μ : Measure ℝ) (A : Set ℝ) : ℝ :=
  (∫ x in A, (x - condMean μ A) ^ 2 ∂μ) / μ.real A

lemma integral_sub_sq_eq_add {ν : Measure ℝ} [IsFiniteMeasure ν]
    (h : MemLp (fun x : ℝ => x) 2 ν) {a : ℝ} (ha : ∫ x, x ∂ν = ν.real Set.univ * a) (t : ℝ) :
    ∫ x, (x - t) ^ 2 ∂ν = ∫ x, (x - a) ^ 2 ∂ν + ν.real Set.univ * (a - t) ^ 2 := by
  have h2 : MemLp (fun x : ℝ => x) ((2 : ℕ) : ℝ≥0∞) ν := by exact_mod_cast h
  have hi : Integrable (fun x : ℝ => x) ν := h.integrable one_le_two
  have hia : Integrable (fun x : ℝ => (x - a) ^ 2) ν := integrable_sub_pow h2 a
  have hib : Integrable (fun x : ℝ => 2 * (a - t) * (x - a)) ν :=
    (hi.sub (integrable_const a)).const_mul _
  have hiab : Integrable (fun x : ℝ => (x - a) ^ 2 + 2 * (a - t) * (x - a)) ν := hia.add hib
  have key : ∫ x, (x - t) ^ 2 ∂ν =
      ∫ x, ((x - a) ^ 2 + 2 * (a - t) * (x - a) + (a - t) ^ 2) ∂ν :=
    integral_congr_ae (Filter.Eventually.of_forall fun x => by ring)
  rw [key, integral_add hiab (integrable_const _), integral_add hia hib, integral_const_mul,
    integral_sub hi (integrable_const a), integral_const, integral_const, ha]
  simp only [smul_eq_mul]
  ring

theorem setIntegral_sub_sq_eq [IsFiniteMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ)
    {A : Set ℝ} (hA : μ A ≠ 0) (t : ℝ) :
    ∫ x in A, (x - t) ^ 2 ∂μ = μ.real A * (condVar μ A + (condMean μ A - t) ^ 2) := by
  have hp : μ.real A ≠ 0 := (measureReal_ne_zero_iff (measure_ne_top μ A)).mpr hA
  have ha : ∫ x in A, x ∂μ = (μ.restrict A).real Set.univ * condMean μ A := by
    rw [measureReal_restrict_apply_univ, condMean]
    field_simp
  rw [integral_sub_sq_eq_add (h.restrict A) ha t, measureReal_restrict_apply_univ, condVar]
  field_simp

theorem mul_condVar_le_setIntegral_sub_sq [IsFiniteMeasure μ]
    (h : MemLp (fun x : ℝ => x) 2 μ) (A : Set ℝ) (t : ℝ) :
    μ.real A * condVar μ A ≤ ∫ x in A, (x - t) ^ 2 ∂μ := by
  by_cases hA : μ A = 0
  · rw [(measureReal_eq_zero_iff (measure_ne_top μ A)).mpr hA, zero_mul]
    exact integral_nonneg fun _ => sq_nonneg _
  · rw [setIntegral_sub_sq_eq h hA t]
    exact mul_le_mul_of_nonneg_left (le_add_of_nonneg_right (sq_nonneg _)) measureReal_nonneg

lemma mLPM_eq_setIntegral_Iio {r : ℕ} (hr : r ≠ 0) (t : ℝ) (μ : Measure ℝ) :
    mLPM r t μ = ∫ x in Set.Iio t, (t - x) ^ r ∂μ := by
  have h : lowerDev r t = (Set.Iio t).indicator fun v => (t - v) ^ r := by
    ext v
    rcases lt_trichotomy v t with h | rfl | h
    · simp [lowerDev, Set.indicator, h, h.le]
    · simp [lowerDev, Set.indicator, zero_pow hr]
    · simp [lowerDev, Set.indicator, not_le.mpr h, not_lt.mpr h.le]
  rw [mLPM, h, integral_indicator measurableSet_Iio]

theorem mLPM_two_eq_cond [IsFiniteMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ) (t : ℝ)
    (hA : μ (Set.Iic t) ≠ 0) :
    mLPM 2 t μ = μ.real (Set.Iic t) *
      (condVar μ (Set.Iic t) + (condMean μ (Set.Iic t) - t) ^ 2) := by
  rw [mLPM_eq_setIntegral, ← setIntegral_sub_sq_eq h hA t]
  exact integral_congr_ae (Filter.Eventually.of_forall fun x => sub_sq_comm t x)

theorem mUPM_two_eq_cond [IsFiniteMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ) (t : ℝ)
    (hA : μ (Set.Ioi t) ≠ 0) :
    mUPM 2 t μ = μ.real (Set.Ioi t) *
      (condVar μ (Set.Ioi t) + (condMean μ (Set.Ioi t) - t) ^ 2) := by
  rw [mUPM_eq_setIntegral, ← setIntegral_sub_sq_eq h hA t]

theorem variance_eq_setIntegral_add_compl [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 2 μ) {A : Set ℝ} (hA : MeasurableSet A) :
    variance (fun x : ℝ => x) μ =
      ∫ x in A, (x - ∫ y, y ∂μ) ^ 2 ∂μ + ∫ x in Aᶜ, (x - ∫ y, y ∂μ) ^ 2 ∂μ := by
  have h2 : MemLp (fun x : ℝ => x) ((2 : ℕ) : ℝ≥0∞) μ := by exact_mod_cast h
  rw [variance_eq, integral_add_compl hA (integrable_sub_pow h2 _)]

theorem variance_eq_within_add_between [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 2 μ) {A : Set ℝ} (hA : MeasurableSet A) (h1 : μ A ≠ 0)
    (h2 : μ Aᶜ ≠ 0) :
    variance (fun x : ℝ => x) μ = (μ.real A * condVar μ A + μ.real Aᶜ * condVar μ Aᶜ) +
      (μ.real A * (condMean μ A - ∫ y, y ∂μ) ^ 2 +
        μ.real Aᶜ * (condMean μ Aᶜ - ∫ y, y ∂μ) ^ 2) := by
  rw [variance_eq_setIntegral_add_compl h hA, setIntegral_sub_sq_eq h h1,
    setIntegral_sub_sq_eq h h2]
  ring

theorem between_nonneg (μ : Measure ℝ) (A : Set ℝ) (m : ℝ) :
    0 ≤ μ.real A * (condMean μ A - m) ^ 2 + μ.real Aᶜ * (condMean μ Aᶜ - m) ^ 2 :=
  add_nonneg (mul_nonneg measureReal_nonneg (sq_nonneg _))
    (mul_nonneg measureReal_nonneg (sq_nonneg _))

theorem within_le_variance [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ)
    {A : Set ℝ} (hA : MeasurableSet A) :
    μ.real A * condVar μ A + μ.real Aᶜ * condVar μ Aᶜ ≤ variance (fun x : ℝ => x) μ := by
  rw [variance_eq_setIntegral_add_compl h hA]
  exact add_le_add (mul_condVar_le_setIntegral_sub_sq h A _)
    (mul_condVar_le_setIntegral_sub_sq h Aᶜ _)

theorem variance_eq_total_variance [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ)
    {A : Set ℝ} (hA : MeasurableSet A) (h1 : μ A ≠ 0) (h2 : μ Aᶜ ≠ 0) :
    variance (fun x : ℝ => x) μ = μ.real A * condVar μ A + (1 - μ.real A) * condVar μ Aᶜ +
      μ.real A * (1 - μ.real A) * (condMean μ A - condMean μ Aᶜ) ^ 2 := by
  have hq : μ.real Aᶜ = 1 - μ.real A := probReal_compl_eq_one_sub hA
  have hp1 : μ.real A ≠ 0 := (measureReal_ne_zero_iff (measure_ne_top μ A)).mpr h1
  have hp2 : μ.real Aᶜ ≠ 0 := (measureReal_ne_zero_iff (measure_ne_top μ Aᶜ)).mpr h2
  have hm : ∫ y, y ∂μ = μ.real A * condMean μ A + μ.real Aᶜ * condMean μ Aᶜ := by
    rw [← integral_add_compl hA (h.integrable one_le_two), condMean, condMean]
    field_simp
  rw [variance_eq_within_add_between h hA h1 h2, hm, hq]
  ring

theorem variance_eq_total_variance_Ici [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 2 μ) (t : ℝ) (h1 : μ (Set.Ici t) ≠ 0)
    (h2 : μ (Set.Iio t) ≠ 0) :
    variance (fun x : ℝ => x) μ =
      μ.real (Set.Ici t) * condVar μ (Set.Ici t) +
        (1 - μ.real (Set.Ici t)) * condVar μ (Set.Iio t) +
        μ.real (Set.Ici t) * (1 - μ.real (Set.Ici t)) *
          (condMean μ (Set.Ici t) - condMean μ (Set.Iio t)) ^ 2 := by
  have := variance_eq_total_variance h (measurableSet_Ici (a := t)) h1 (by rwa [Set.compl_Ici])
  rwa [Set.compl_Ici] at this

theorem within_le_variance_Ici [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ)
    (t : ℝ) :
    μ.real (Set.Ici t) * condVar μ (Set.Ici t) +
      (1 - μ.real (Set.Ici t)) * condVar μ (Set.Iio t) ≤ variance (fun x : ℝ => x) μ := by
  have := within_le_variance h (measurableSet_Ici (a := t))
  rwa [probReal_compl_eq_one_sub measurableSet_Ici, Set.compl_Ici] at this

lemma measure_Iic_ne_zero_of_mLPM_ne_zero {r : ℕ} {t : ℝ} (h : mLPM r t μ ≠ 0) :
    μ (Set.Iic t) ≠ 0 := fun h0 =>
  h (by rw [mLPM_eq_setIntegral, Measure.restrict_eq_zero.mpr h0, integral_zero_measure])

lemma measure_Ioi_ne_zero_of_mUPM_ne_zero {r : ℕ} {t : ℝ} (h : mUPM r t μ ≠ 0) :
    μ (Set.Ioi t) ≠ 0 := fun h0 =>
  h (by rw [mUPM_eq_setIntegral, Measure.restrict_eq_zero.mpr h0, integral_zero_measure])

theorem variance_split_at_mean [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ)
    (hv : 0 < variance (fun x : ℝ => x) μ) :
    mLPM 2 (∫ y, y ∂μ) μ = μ.real (Set.Iic (∫ y, y ∂μ)) *
        (condVar μ (Set.Iic (∫ y, y ∂μ)) + (condMean μ (Set.Iic (∫ y, y ∂μ)) - ∫ y, y ∂μ) ^ 2) ∧
      mUPM 2 (∫ y, y ∂μ) μ = μ.real (Set.Ioi (∫ y, y ∂μ)) *
        (condVar μ (Set.Ioi (∫ y, y ∂μ)) + (condMean μ (Set.Ioi (∫ y, y ∂μ)) - ∫ y, y ∂μ) ^ 2) ∧
      variance (fun x : ℝ => x) μ =
        (μ.real (Set.Iic (∫ y, y ∂μ)) * condVar μ (Set.Iic (∫ y, y ∂μ)) +
          μ.real (Set.Ioi (∫ y, y ∂μ)) * condVar μ (Set.Ioi (∫ y, y ∂μ))) +
        (μ.real (Set.Iic (∫ y, y ∂μ)) * (condMean μ (Set.Iic (∫ y, y ∂μ)) - ∫ y, y ∂μ) ^ 2 +
          μ.real (Set.Ioi (∫ y, y ∂μ)) * (condMean μ (Set.Ioi (∫ y, y ∂μ)) - ∫ y, y ∂μ) ^ 2) := by
  obtain ⟨hL, hU⟩ := mLPM_two_pos_and_mUPM_two_pos h hv
  have eL := mLPM_two_eq_cond h _ (measure_Iic_ne_zero_of_mLPM_ne_zero hL.ne')
  have eU := mUPM_two_eq_cond h _ (measure_Ioi_ne_zero_of_mUPM_ne_zero hU.ne')
  refine ⟨eL, eU, ?_⟩
  rw [variance_eq_mUPM_two_add_mLPM_two h, eL, eU]
  ring

theorem variance_eq_total_variance_mean [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 2 μ) (hv : 0 < variance (fun x : ℝ => x) μ) :
    variance (fun x : ℝ => x) μ =
      μ.real (Set.Ici (∫ y, y ∂μ)) * condVar μ (Set.Ici (∫ y, y ∂μ)) +
        (1 - μ.real (Set.Ici (∫ y, y ∂μ))) * condVar μ (Set.Iio (∫ y, y ∂μ)) +
        μ.real (Set.Ici (∫ y, y ∂μ)) * (1 - μ.real (Set.Ici (∫ y, y ∂μ))) *
          (condMean μ (Set.Ici (∫ y, y ∂μ)) - condMean μ (Set.Iio (∫ y, y ∂μ))) ^ 2 := by
  obtain ⟨hL, hU⟩ := mLPM_two_pos_and_mUPM_two_pos h hv
  have h1 : μ (Set.Ici (∫ y, y ∂μ)) ≠ 0 := fun h0 =>
    measure_Ioi_ne_zero_of_mUPM_ne_zero hU.ne' (measure_mono_null Set.Ioi_subset_Ici_self h0)
  have h2 : μ (Set.Iio (∫ y, y ∂μ)) ≠ 0 := fun h0 => hL.ne' (by
    rw [mLPM_eq_setIntegral_Iio two_ne_zero, Measure.restrict_eq_zero.mpr h0,
      integral_zero_measure])
  exact variance_eq_total_variance_Ici h _ h1 h2

def threePoint (a : ℝ) : Fin 18 → ℝ := Fin.cons (-a) (Fin.cons a fun _ => 0)

lemma sum_threePoint (a : ℝ) (g : ℝ → ℝ) :
    ∑ i, g (threePoint a i) = g (-a) + g a + 16 * g 0 := by
  rw [Fin.sum_univ_succ, Fin.sum_univ_succ]
  simp only [threePoint, Fin.cons_zero, Fin.cons_succ, sum_const, card_univ, Fintype.card_fin,
    nsmul_eq_mul, Nat.cast_ofNat]
  ring

theorem threePoint_bounded (a : ℝ) (i : Fin 18) : |threePoint a i| ≤ |a| := by
  refine Fin.cases ?_ (fun j => Fin.cases ?_ (fun k => ?_) j) i
  · simp [threePoint]
  · simp [threePoint]
  · simp [threePoint]

theorem threePoint_mean (a : ℝ) : mean (threePoint a) = 0 := by
  have h := sum_threePoint a fun v => v
  rw [mean, h]
  ring

theorem threePoint_popVariance (a : ℝ) : popVariance (threePoint a) = a ^ 2 / 9 := by
  have h := sum_threePoint a fun v => (v - 0) ^ 2
  rw [popVariance, centralMoment, threePoint_mean, h]
  push_cast
  ring

theorem threePoint_centralMoment_four (a : ℝ) :
    centralMoment 4 (mean (threePoint a)) (threePoint a) = a ^ 4 / 9 := by
  have h := sum_threePoint a fun v => (v - 0) ^ 4
  rw [centralMoment, threePoint_mean, h]
  push_cast
  ring

theorem threePoint_kurtosis {a : ℝ} (ha : a ≠ 0) :
    0 < popVariance (threePoint a) ∧ popKurtosis (threePoint a) = 6 ∧
      mKurtosis (empMeasure (threePoint a)) = 9 := by
  have hk : popKurtosis (threePoint a) = 6 := by
    rw [popKurtosis, threePoint_centralMoment_four, threePoint_popVariance]
    field_simp
    norm_num
  refine ⟨?_, hk, ?_⟩
  · rw [threePoint_popVariance]
    positivity
  · rw [mKurtosis_empMeasure, hk]
    norm_num

theorem mean_neg (x : Fin n → ℝ) : mean (-x) = -mean x := by
  simp [mean, neg_div]

lemma upperDev_neg {d : ℕ} (hd : d ≠ 0) (t v : ℝ) : upperDev d (-t) (-v) = lowerDev d t v := by
  rcases lt_trichotomy v t with h | rfl | h
  · simp [upperDev, lowerDev, h, h.le, neg_add_eq_sub]
  · simp [upperDev, lowerDev, zero_pow hd]
  · simp [upperDev, lowerDev, not_le.mpr h, not_lt.mpr h.le]

lemma lowerDev_neg {d : ℕ} (hd : d ≠ 0) (t v : ℝ) : lowerDev d (-t) (-v) = upperDev d t v := by
  rcases lt_trichotomy t v with h | rfl | h
  · simp [upperDev, lowerDev, h, h.le, neg_add_eq_sub]
  · simp [upperDev, lowerDev, zero_pow hd]
  · simp [upperDev, lowerDev, not_le.mpr h, not_lt.mpr h.le]

theorem UPM_neg {d : ℕ} (hd : d ≠ 0) (t : ℝ) (x : Fin n → ℝ) :
    UPM d (-t) (-x) = LPM d t x := by
  simp only [UPM, LPM, Pi.neg_apply, upperDev_neg hd]

theorem LPM_neg {d : ℕ} (hd : d ≠ 0) (t : ℝ) (x : Fin n → ℝ) :
    LPM d (-t) (-x) = UPM d t x := by
  simp only [UPM, LPM, Pi.neg_apply, lowerDev_neg hd]

theorem tails_four_neg (x : Fin n → ℝ) :
    UPM 4 (mean (-x)) (-x) = LPM 4 (mean x) x ∧ LPM 4 (mean (-x)) (-x) = UPM 4 (mean x) x := by
  rw [mean_neg]
  exact ⟨UPM_neg (by norm_num) _ x, LPM_neg (by norm_num) _ x⟩

theorem popKurtosis_neg (x : Fin n → ℝ) :
    popVariance (-x) = popVariance x ∧ popKurtosis (-x) = popKurtosis x := by
  have h : -x = fun i => -1 * x i + 0 := by
    ext i
    simp
  rw [h, popVariance_affine, popKurtosis_affine _ _ _ (by norm_num)]
  norm_num

def tailSample : Fin 3 → ℝ := ![-1, -1, 2]

theorem tailSample_witness :
    0 < popVariance tailSample ∧ popVariance (-tailSample) = popVariance tailSample ∧
      popKurtosis (-tailSample) = popKurtosis tailSample ∧
      UPM 4 (mean tailSample) tailSample = 16 / 3 ∧
      LPM 4 (mean tailSample) tailSample = 2 / 3 ∧
      UPM 4 (mean (-tailSample)) (-tailSample) = 2 / 3 ∧
      LPM 4 (mean (-tailSample)) (-tailSample) = 16 / 3 := by
  have hm : mean tailSample = 0 := by
    norm_num [mean, tailSample, Fin.sum_univ_three]
  have hv : popVariance tailSample = 2 := by
    rw [popVariance, centralMoment, hm]
    norm_num [tailSample, Fin.sum_univ_three]
  have hU : UPM 4 (mean tailSample) tailSample = 16 / 3 := by
    rw [hm]
    norm_num [UPM, upperDev, tailSample, Fin.sum_univ_three]
  have hL : LPM 4 (mean tailSample) tailSample = 2 / 3 := by
    rw [hm]
    norm_num [LPM, lowerDev, tailSample, Fin.sum_univ_three]
  obtain ⟨h1, h2⟩ := tails_four_neg tailSample
  obtain ⟨h3, h4⟩ := popKurtosis_neg tailSample
  exact ⟨by rw [hv]; norm_num, h3, h4, hU, hL, h1.trans hL, h2.trans hU⟩

section StudentT

open Filter Topology Set

def t5Kernel (y : ℝ) : ℝ := ((1 + y ^ 2) ^ 3)⁻¹

lemma t5Kernel_pos (y : ℝ) : 0 < t5Kernel y := by unfold t5Kernel; positivity

lemma continuous_t5Kernel : Continuous t5Kernel :=
  (by fun_prop : Continuous fun y : ℝ => (1 + y ^ 2) ^ 3).inv₀ fun y => by positivity

lemma abs_pow_le_one_add_sq_sq {k : ℕ} (hk : k ≤ 4) (y : ℝ) : |y| ^ k ≤ (1 + y ^ 2) ^ 2 := by
  rcases le_total |y| 1 with h | h
  · calc |y| ^ k ≤ 1 := pow_le_one₀ (abs_nonneg y) h
      _ ≤ (1 + y ^ 2) ^ 2 := one_le_pow₀ (by nlinarith [sq_nonneg y])
  · calc |y| ^ k ≤ |y| ^ 4 := pow_le_pow_right₀ h hk
      _ = (y ^ 2) ^ 2 := by rw [← sq_abs y]; ring
      _ ≤ (1 + y ^ 2) ^ 2 := by gcongr; linarith

lemma integrable_pow_mul_t5Kernel {k : ℕ} (hk : k ≤ 4) :
    Integrable fun y => y ^ k * t5Kernel y := by
  refine integrable_inv_one_add_sq.mono'
    ((continuous_pow k).mul continuous_t5Kernel).aestronglyMeasurable (ae_of_all _ fun y => ?_)
  rw [Real.norm_eq_abs, abs_mul, abs_pow, abs_of_pos (t5Kernel_pos y), t5Kernel]
  have hu : 0 < 1 + y ^ 2 := by positivity
  calc |y| ^ k * ((1 + y ^ 2) ^ 3)⁻¹ ≤ (1 + y ^ 2) ^ 2 * ((1 + y ^ 2) ^ 3)⁻¹ := by
        gcongr; exact abs_pow_le_one_add_sq_sq hk y
    _ = (1 + y ^ 2)⁻¹ := by field_simp

lemma tendsto_zero_of_le_inv_abs {φ : ℝ → ℝ} (hb : ∀ y, 1 ≤ |y| → |φ y| ≤ |y|⁻¹) :
    Tendsto φ atBot (𝓝 0) ∧ Tendsto φ atTop (𝓝 0) := by
  have key : ∀ l : Filter ℝ, Tendsto (fun y : ℝ => |y|) l atTop →
      Tendsto φ l (𝓝 0) := fun l hl => by
    refine squeeze_zero_norm' ?_ (tendsto_inv_atTop_zero.comp hl)
    filter_upwards [hl.eventually (eventually_ge_atTop 1)] with y hy
    exact hb y hy
  exact ⟨key _ tendsto_abs_atBot_atTop, key _ tendsto_abs_atTop_atTop⟩

lemma integral_eq_zero_of_deriv {φ φ' : ℝ → ℝ} (hd : ∀ y, HasDerivAt φ (φ' y) y)
    (hi : Integrable φ') (hb : ∀ y, 1 ≤ |y| → |φ y| ≤ |y|⁻¹) : ∫ y, φ' y = 0 := by
  obtain ⟨h1, h2⟩ := tendsto_zero_of_le_inv_abs hb
  simpa using integral_of_hasDerivAt_of_tendsto hd hi h1 h2

lemma hasDerivAt_pow_div (k : ℕ) (y : ℝ) :
    HasDerivAt (fun y : ℝ => y ^ (k + 1) / (1 + y ^ 2) ^ 2)
      ((k + 1) * (y ^ k * t5Kernel y) + (k - 3) * (y ^ (k + 2) * t5Kernel y)) y := by
  have hu : (1 + y ^ 2) ≠ 0 := by positivity
  have h := (hasDerivAt_pow (k + 1) y).div (((hasDerivAt_pow 2 y).const_add 1).pow 2)
    (pow_ne_zero 2 hu)
  convert h using 1
  simp only [t5Kernel, Pi.pow_apply, Nat.cast_ofNat, Nat.cast_add, Nat.cast_one, Nat.add_sub_cancel]
  field_simp
  ring

lemma abs_pow_div_le {k : ℕ} (hk : k ≤ 2) (y : ℝ) (hy : 1 ≤ |y|) :
    |y ^ (k + 1) / (1 + y ^ 2) ^ 2| ≤ |y|⁻¹ := by
  have hy0 : 0 < |y| := by linarith
  rw [abs_div, abs_pow, abs_of_pos (by positivity : (0 : ℝ) < (1 + y ^ 2) ^ 2),
    div_le_iff₀ (by positivity), inv_mul_eq_div, le_div_iff₀ hy0, ← pow_succ]
  exact abs_pow_le_one_add_sq_sq (by omega) y

lemma t5_recursion {k : ℕ} (hk : k ≤ 2) :
    (k + 1) * ∫ y, y ^ k * t5Kernel y = (3 - k) * ∫ y, y ^ (k + 2) * t5Kernel y := by
  have ik := integrable_pow_mul_t5Kernel (k := k) (by omega)
  have ik2 := integrable_pow_mul_t5Kernel (k := k + 2) (by omega)
  have h0 := integral_eq_zero_of_deriv (hasDerivAt_pow_div k)
    ((ik.const_mul _).add (ik2.const_mul _)) fun y hy => abs_pow_div_le hk y hy
  rw [integral_add (ik.const_mul _) (ik2.const_mul _), integral_const_mul,
    integral_const_mul] at h0
  linarith

lemma integrable_t5Kernel : Integrable t5Kernel := by
  simpa using integrable_pow_mul_t5Kernel (k := 0) (by norm_num)

lemma integral_t5Kernel_pos : 0 < ∫ y, t5Kernel y := by
  refine (integral_pos_iff_support_of_nonneg (fun y => (t5Kernel_pos y).le)
    integrable_t5Kernel).mpr ?_
  have : Function.support t5Kernel = univ := by
    ext y; simp [(t5Kernel_pos y).ne']
  simp [this]

lemma t5_moments :
    ∫ y, y ^ 2 * t5Kernel y = (∫ y, t5Kernel y) / 3 ∧
      ∫ y, y ^ 4 * t5Kernel y = ∫ y, t5Kernel y := by
  have h0 := t5_recursion (k := 0) (by norm_num)
  have h2 := t5_recursion (k := 2) (by norm_num)
  simp only [pow_zero, one_mul, CharP.cast_eq_zero, zero_add, Nat.cast_ofNat] at h0 h2
  constructor <;> linarith

lemma integral_id_mul_t5Kernel : ∫ y, y * t5Kernel y = 0 := by
  have h := integral_neg_eq_self (fun y => y * t5Kernel y) volume
  simp only [t5Kernel, even_two.neg_pow, neg_mul, integral_neg] at h
  simp only [t5Kernel]
  linarith

def scaledT5 : Measure ℝ :=
  volume.withDensity fun y => ENNReal.ofReal (t5Kernel y / ∫ z, t5Kernel z)

lemma integral_scaledT5 (g : ℝ → ℝ) :
    ∫ y, g y ∂scaledT5 = (∫ y, g y * t5Kernel y) / ∫ z, t5Kernel z := by
  rw [scaledT5, integral_withDensity_eq_integral_toReal_smul
    (continuous_t5Kernel.div_const _).measurable.ennreal_ofReal
    (Eventually.of_forall fun _ => ENNReal.ofReal_lt_top), ← integral_div]
  refine integral_congr_ae (Eventually.of_forall fun y => ?_)
  simp only [smul_eq_mul, ENNReal.toReal_ofReal
    (div_nonneg (t5Kernel_pos y).le integral_t5Kernel_pos.le)]
  ring

lemma integral_id_scaledT5 : ∫ y, y ∂scaledT5 = 0 := by
  rw [integral_scaledT5, integral_id_mul_t5Kernel, zero_div]

theorem mKurtosis_scaledT5 : mKurtosis scaledT5 = 9 := by
  obtain ⟨h2, h4⟩ := t5_moments
  have hZ := integral_t5Kernel_pos.ne'
  rw [mKurtosis, variance_eq_integral (X := fun x : ℝ => x) aemeasurable_id',
    integral_id_scaledT5]
  simp only [sub_zero]
  rw [integral_scaledT5, integral_scaledT5, h2, h4]
  field_simp
  norm_num

theorem scaledT5_unbounded (M : ℝ) : scaledT5 {y | M < |y|} ≠ 0 := by
  rw [scaledT5, Ne, withDensity_apply_eq_zero'
    (continuous_t5Kernel.div_const _).measurable.ennreal_ofReal.aemeasurable]
  have hpos : ∀ y, ENNReal.ofReal (t5Kernel y / ∫ z, t5Kernel z) ≠ 0 := fun y =>
    (ENNReal.ofReal_pos.mpr (div_pos (t5Kernel_pos y) integral_t5Kernel_pos)).ne'
  have hset : {y | ENNReal.ofReal (t5Kernel y / ∫ z, t5Kernel z) ≠ 0} ∩ {y : ℝ | M < |y|} ⊇
      Ioi |M| := fun y hy => ⟨hpos y, (le_abs_self M).trans_lt (hy.trans_le (le_abs_self y))⟩
  intro h0
  have := measure_mono_null hset h0
  simp at this

end StudentT

end

end NNS.GapMoments
