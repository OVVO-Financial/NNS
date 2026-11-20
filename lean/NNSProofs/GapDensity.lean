-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.CoConditional
public import NNSProofs.GapQuantile
public import NNSProofs.GapSampling
public import Mathlib.Probability.Distributions.Exponential
public import Mathlib.Probability.Distributions.Uniform
public import Mathlib.Analysis.Complex.ExponentialBounds
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.AbsolutelyContinuousFun
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.LebesgueDifferentiationThm

namespace NNS.GapDensity

open MeasureTheory ProbabilityTheory Set Filter Topology
open scoped ENNReal

@[expose] public noncomputable section

theorem integral_withDensity_ofReal {f : ℝ → ℝ} (hf : Measurable f) (hf0 : ∀ x, 0 ≤ f x)
    (g : ℝ → ℝ) :
    ∫ x, g x ∂(volume.withDensity fun x => ENNReal.ofReal (f x)) = ∫ x, g x * f x := by
  rw [integral_withDensity_eq_integral_toReal_smul hf.ennreal_ofReal
    (Eventually.of_forall fun _ => ENNReal.ofReal_lt_top)]
  refine integral_congr_ae (Eventually.of_forall fun x => ?_)
  simp only [ENNReal.toReal_ofReal (hf0 x), smul_eq_mul, mul_comm]

theorem expectation_eq_integral_mul_density {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    {X : Ω → ℝ} {f g : ℝ → ℝ} (hX : AEMeasurable X P) (hf : Measurable f) (hf0 : ∀ x, 0 ≤ f x)
    (hlaw : P.map X = volume.withDensity fun x => ENNReal.ofReal (f x))
    (hg : AEStronglyMeasurable g (P.map X)) : ∫ ω, g (X ω) ∂P = ∫ x, g x * f x := by
  rw [← integral_map hX hg, hlaw, integral_withDensity_ofReal hf hf0]

theorem expectation_eq_integral_cdf_measure {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    [IsProbabilityMeasure P] {X : Ω → ℝ} {g : ℝ → ℝ} (hX : AEMeasurable X P)
    (hg : AEStronglyMeasurable g (P.map X)) :
    ∫ ω, g (X ω) ∂P = ∫ x, g x ∂(cdf (P.map X)).measure := by
  rw [measure_cdf, integral_map hX hg]

section Uniform

variable {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {X : Ω → ℝ} {a b : ℝ} {f : ℝ → ℝ}

theorem IsUniform.expectation_eq (hab : a ≤ b) (hX : pdf.IsUniform X (Icc a b) P)
    (hf : AEStronglyMeasurable f (volume.restrict (Icc a b))) :
    ∫ ω, f (X ω) ∂P = (b - a)⁻¹ * ∫ x in a..b, f x := by
  have hmap : P.map X = (volume (Icc a b))⁻¹ • volume.restrict (Icc a b) := hX.map_eq
  rw [← integral_map hX.aemeasurable (by rw [hmap]; exact hf.smul_measure _), hmap,
    integral_smul_measure, Real.volume_Icc, ENNReal.toReal_inv,
    ENNReal.toReal_ofReal (sub_nonneg.mpr hab), intervalIntegral.integral_of_le hab,
    integral_Icc_eq_integral_Ioc, smul_eq_mul]

theorem IsUniform.intervalIntegral_eq (hab : a < b) (hX : pdf.IsUniform X (Icc a b) P)
    (hf : AEStronglyMeasurable f (volume.restrict (Icc a b))) :
    ∫ x in a..b, f x = (b - a) * ∫ ω, f (X ω) ∂P := by
  rw [IsUniform.expectation_eq hab.le hX hf, mul_inv_cancel_left₀ (sub_ne_zero.mpr hab.ne')]

lemma IsUniform.aux (hab : a < b) (hX : pdf.IsUniform X (Icc a b) P)
    (hf : IntervalIntegrable f volume a b) :
    IsProbabilityMeasure P ∧ Integrable (fun ω => f (X ω)) P ∧
      AEStronglyMeasurable f (volume.restrict (Icc a b)) := by
  have hIcc : IntegrableOn f (Icc a b) := by
    rw [integrableOn_Icc_iff_integrableOn_Ioc]; exact hf.1
  have hne : volume (Icc a b) ≠ 0 := by
    rw [Real.volume_Icc]; exact (ENNReal.ofReal_pos.mpr (sub_pos.mpr hab)).ne'
  have hmap : P.map X = (volume (Icc a b))⁻¹ • volume.restrict (Icc a b) := hX.map_eq
  refine ⟨hX.isProbabilityMeasure hne (by simp), ?_, hIcc.aestronglyMeasurable⟩
  have h1 : Integrable f (P.map X) := by
    rw [hmap]; exact hIcc.smul_measure (ENNReal.inv_ne_top.mpr hne)
  exact (integrable_map_measure h1.aestronglyMeasurable hX.aemeasurable).mp h1

theorem IsUniform.intervalIntegral_eq_partialMoments (hab : a < b)
    (hX : pdf.IsUniform X (Icc a b) P) (hf : IntervalIntegrable f volume a b) :
    ∫ x in a..b, f x = (b - a) *
      (mUPM 1 0 (P.map fun ω => f (X ω)) - mLPM 1 0 (P.map fun ω => f (X ω))) := by
  obtain ⟨hP, hint, hfm⟩ := IsUniform.aux hab hX hf
  have hY : AEMeasurable (fun ω => f (X ω)) P := hint.aemeasurable
  have hid : Integrable (fun y : ℝ => y) (P.map fun ω => f (X ω)) :=
    (integrable_map_measure aestronglyMeasurable_id hY).mpr hint
  have hm := mean_eq_mUPM_one_sub_mLPM_one hid
  rw [integral_map (f := fun y : ℝ => y) hY aestronglyMeasurable_id] at hm
  rw [← hm, IsUniform.intervalIntegral_eq hab hX hfm]

theorem IsUniform.intervalIntegral_abs_eq (hab : a < b)
    (hX : pdf.IsUniform X (Icc a b) P) (hf : IntervalIntegrable f volume a b) :
    ∫ x in a..b, |f x| = (b - a) *
      (mUPM 1 0 (P.map fun ω => f (X ω)) + mLPM 1 0 (P.map fun ω => f (X ω))) := by
  obtain ⟨hP, hint, hfm⟩ := IsUniform.aux hab hX hf
  have hY : AEMeasurable (fun ω => f (X ω)) P := hint.aemeasurable
  have hid : Integrable (fun y : ℝ => y) (P.map fun ω => f (X ω)) :=
    (integrable_map_measure aestronglyMeasurable_id hY).mpr hint
  have habs : AEStronglyMeasurable (fun x => |f x|) (volume.restrict (Icc a b)) :=
    (continuous_abs.comp_aestronglyMeasurable hfm)
  rw [mUPM_one_add_mLPM_one hid, IsUniform.intervalIntegral_eq hab hX habs]
  simp_rw [sub_zero]
  rw [integral_map hY (by fun_prop)]

lemma UPM_one_add_LPM_one_zero {n : ℕ} (y : Fin n → ℝ) :
    UPM 1 0 y + LPM 1 0 y = (∑ i, |y i|) / n := by
  simp only [UPM, LPM, ← add_div, ← Finset.sum_add_distrib, upperDev_add_lowerDev, sub_zero,
    pow_one]

theorem ae_tendsto_uniform_estimate [IsProbabilityMeasure P] (hab : a < b) {X : ℕ → Ω → ℝ}
    (hX : ∀ i, pdf.IsUniform (X i) (Icc a b) P) (hindep : Pairwise fun i j => X i ⟂ᵢ[P] X j)
    (hfm : Measurable f) (hf : IntervalIntegrable f volume a b) :
    ∀ᵐ ω ∂P,
      Tendsto (fun n : ℕ => (b - a) * (UPM 1 0 (fun i : Fin n => f (X i ω)) -
        LPM 1 0 (fun i : Fin n => f (X i ω)))) atTop (𝓝 (∫ x in a..b, f x)) ∧
      Tendsto (fun n : ℕ => (b - a) * (UPM 1 0 (fun i : Fin n => f (X i ω)) +
        LPM 1 0 (fun i : Fin n => f (X i ω)))) atTop (𝓝 (∫ x in a..b, |f x|)) := by
  have hident : ∀ i, IdentDistrib (X i) (X 0) P P := fun i =>
    ⟨(hX i).aemeasurable, (hX 0).aemeasurable, by rw [(hX i).map_eq, (hX 0).map_eq]⟩
  have hint : Integrable (fun ω => f (X 0 ω)) P := (IsUniform.aux hab (hX 0) hf).2.1
  have habs : Measurable fun x => |f x| := continuous_abs.measurable.comp hfm
  have h1 := GapSampling.ae_tendsto_average X hfm hint hindep hident
  have h2 := GapSampling.ae_tendsto_average X habs hint.abs hindep hident
  rw [integral_map (hX 0).aemeasurable hfm.aestronglyMeasurable] at h1
  rw [integral_map (hX 0).aemeasurable habs.aestronglyMeasurable] at h2
  have e1 : ∫ x in a..b, f x = (b - a) * ∫ ω, f (X 0 ω) ∂P :=
    IsUniform.intervalIntegral_eq hab (hX 0) hfm.aestronglyMeasurable
  have e2 : ∫ x in a..b, |f x| = (b - a) * ∫ ω, |f (X 0 ω)| ∂P :=
    IsUniform.intervalIntegral_eq hab (hX 0) habs.aestronglyMeasurable
  filter_upwards [h1, h2] with ω hω1 hω2
  rw [e1, e2]
  refine ⟨(hω1.const_mul (b - a)).congr' ?_, (hω2.const_mul (b - a)).congr' ?_⟩ <;>
    filter_upwards [eventually_gt_atTop 0] with n hn
  · rw [UPM_one_sub_LPM_one hn, sub_zero, mean]
  · rw [UPM_one_add_LPM_one_zero]

theorem ae_tendsto_uniform_estimate_of_intervalIntegrable (hab : a < b) {X : ℕ → Ω → ℝ}
    (hX : ∀ i, pdf.IsUniform (X i) (Icc a b) P) (hindep : Pairwise fun i j => X i ⟂ᵢ[P] X j)
    (hf : IntervalIntegrable f volume a b) :
    ∀ᵐ ω ∂P,
      Tendsto (fun n : ℕ => (b - a) * (UPM 1 0 (fun i : Fin n => f (X i ω)) -
        LPM 1 0 (fun i : Fin n => f (X i ω)))) atTop (𝓝 (∫ x in a..b, f x)) ∧
      Tendsto (fun n : ℕ => (b - a) * (UPM 1 0 (fun i : Fin n => f (X i ω)) +
        LPM 1 0 (fun i : Fin n => f (X i ω)))) atTop (𝓝 (∫ x in a..b, |f x|)) := by
  obtain ⟨hP, -, hm⟩ := IsUniform.aux hab (hX 0) hf
  set g := hm.mk f
  have hg : Measurable g := hm.stronglyMeasurable_mk.measurable
  have he : f =ᵐ[volume.restrict (Icc a b)] g := hm.ae_eq_mk
  have hi : f =ᵐ[volume.restrict (uIoc a b)] g := by
    rw [uIoc_of_le hab.le]
    exact ae_mono (Measure.restrict_mono Ioc_subset_Icc_self le_rfl) he
  have heq (i : ℕ) : (fun ω => f (X i ω)) =ᵐ[P] fun ω => g (X i ω) := by
    apply ae_eq_comp (hX i).aemeasurable
    rw [(hX i).map_eq]
    exact Measure.ae_smul_measure he _
  have hlim := ae_tendsto_uniform_estimate hab hX hindep hg (hf.congr_ae hi)
  have habs := intervalIntegral.integral_congr_ae_restrict (hi.fun_comp abs)
  simp only [Function.comp_apply] at habs
  rw [intervalIntegral.integral_congr_ae_restrict hi, habs]
  filter_upwards [ae_all_iff.mpr heq, hlim] with ω hω h
  simpa only [hω] using h

end Uniform

def kde {n : ℕ} (K : ℝ → ℝ) (h : ℝ) (x : Fin n → ℝ) (t : ℝ) : ℝ :=
  ((n : ℝ) * h)⁻¹ * ∑ i, K ((t - x i) / h)

theorem kde_nonneg {n : ℕ} {K : ℝ → ℝ} {h : ℝ} (hh : 0 ≤ h) (hK : ∀ u, 0 ≤ K u)
    (x : Fin n → ℝ) (t : ℝ) : 0 ≤ kde K h x t :=
  mul_nonneg (inv_nonneg.mpr (mul_nonneg n.cast_nonneg hh)) (Finset.sum_nonneg fun _ _ => hK _)

theorem integral_kde {n : ℕ} (hn : 0 < n) {K : ℝ → ℝ} {h : ℝ} (hh : 0 < h)
    (hK : ∫ u, K u = 1) (x : Fin n → ℝ) : ∫ t, kde K h x t = 1 := by
  have hKi : Integrable K := by
    by_contra hni
    rw [integral_undef hni] at hK
    exact zero_ne_one hK
  have hterm : ∀ i, Integrable fun t => K ((t - x i) / h) := fun i =>
    (hKi.comp_div hh.ne').comp_sub_right (x i)
  have hval : ∀ i, ∫ t, K ((t - x i) / h) = h := fun i => by
    rw [integral_sub_right_eq_self (fun t => K (t / h)) (x i), Measure.integral_comp_div, hK,
      abs_of_pos hh, smul_eq_mul, mul_one]
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simp_rw [kde]
  rw [integral_const_mul, integral_finsetSum _ fun i _ => hterm i]
  simp_rw [hval]
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul,
    inv_mul_cancel₀ (mul_ne_zero hn' hh.ne')]

def densityLaw (f : ℝ → ℝ) : Measure ℝ := volume.withDensity fun x => ENNReal.ofReal (f x)

lemma mUPM_zero_eq_one_sub_cdf {μ : Measure ℝ} [IsProbabilityMeasure μ] (t : ℝ) :
    mUPM 0 t μ = 1 - cdf μ t := by
  rw [← mLPM_zero_eq_cdf, ← mLPM_zero_add_mUPM_zero (μ := μ) t, add_sub_cancel_left]

theorem hasDerivAt_neg_log_mUPM_zero {μ : Measure ℝ} [IsProbabilityMeasure μ] {t f' : ℝ}
    (hF : HasDerivAt (fun s => cdf μ s) f' t) (hS : 0 < mUPM 0 t μ) :
    HasDerivAt (fun s => -Real.log (mUPM 0 s μ)) (f' / mUPM 0 t μ) t := by
  have h1 : HasDerivAt (fun s => mUPM 0 s μ) (-f') t := by
    simp_rw [mUPM_zero_eq_one_sub_cdf]
    exact hF.const_sub 1
  convert (h1.log hS.ne').neg using 1
  ring

section Density

variable {f : ℝ → ℝ}

lemma integrable_density [IsFiniteMeasure (densityLaw f)] (hf : Measurable f)
    (hf0 : ∀ x, 0 ≤ f x) : Integrable f := by
  refine ⟨hf.aestronglyMeasurable, (hasFiniteIntegral_iff_ofReal (ae_of_all _ hf0)).mpr ?_⟩
  have h := measure_lt_top (densityLaw f) univ
  rwa [densityLaw, withDensity_apply _ MeasurableSet.univ, Measure.restrict_univ] at h

lemma cdf_densityLaw [IsProbabilityMeasure (densityLaw f)] (hf : Measurable f)
    (hf0 : ∀ x, 0 ≤ f x) (t : ℝ) : cdf (densityLaw f) t = ∫ x in Iic t, f x := by
  rw [cdf_eq_real, measureReal_def, densityLaw, withDensity_apply _ measurableSet_Iic,
    ← ofReal_integral_eq_lintegral_ofReal (integrable_density hf hf0).integrableOn
      (ae_of_all _ hf0),
    ENNReal.toReal_ofReal (integral_nonneg hf0)]

lemma cdf_densityLaw_eq_add [IsProbabilityMeasure (densityLaw f)] (hf : Measurable f)
    (hf0 : ∀ x, 0 ≤ f x) (a s : ℝ) :
    cdf (densityLaw f) s = (∫ x in Iic a, f x) + ∫ x in a..s, f x := by
  have hint := integrable_density hf hf0
  rw [cdf_densityLaw hf hf0, ← intervalIntegral.integral_Iic_sub_Iic hint.integrableOn
    hint.integrableOn, add_sub_cancel]

lemma hasDerivAt_cdf_densityLaw [IsProbabilityMeasure (densityLaw f)] (hf : Measurable f)
    (hf0 : ∀ x, 0 ≤ f x) {t : ℝ} (hc : ContinuousAt f t) :
    HasDerivAt (fun s => cdf (densityLaw f) s) (f t) t := by
  simp_rw [cdf_densityLaw_eq_add hf hf0 t]
  exact (intervalIntegral.integral_hasDerivAt_right IntervalIntegrable.refl
    hf.stronglyMeasurable.stronglyMeasurableAtFilter hc).const_add _

theorem hasDerivAt_cumHazard [IsProbabilityMeasure (densityLaw f)] (hf : Measurable f)
    (hf0 : ∀ x, 0 ≤ f x) {t : ℝ} (hc : ContinuousAt f t) (hS : 0 < mUPM 0 t (densityLaw f)) :
    HasDerivAt (fun s => -Real.log (mUPM 0 s (densityLaw f)))
      (f t / mUPM 0 t (densityLaw f)) t :=
  hasDerivAt_neg_log_mUPM_zero (hasDerivAt_cdf_densityLaw hf hf0 hc) hS

lemma continuous_mUPM_zero_densityLaw [IsProbabilityMeasure (densityLaw f)] (hf : Measurable f)
    (hf0 : ∀ x, 0 ≤ f x) : Continuous fun s => mUPM 0 s (densityLaw f) := by
  simp_rw [mUPM_zero_eq_one_sub_cdf, cdf_densityLaw_eq_add hf hf0 0]
  exact continuous_const.sub (continuous_const.add (intervalIntegral.continuous_primitive
    (fun _ _ => (integrable_density hf hf0).intervalIntegrable) 0))

lemma mUPM_zero_antitone {μ : Measure ℝ} [IsProbabilityMeasure μ] :
    Antitone fun s => mUPM 0 s μ := fun s t hst => by
  simp_rw [mUPM_zero_eq_one_sub_cdf]
  exact sub_le_sub_left (monotone_cdf μ hst) 1

lemma integrableOn_hazard [IsProbabilityMeasure (densityLaw f)] (hf : Measurable f)
    (hf0 : ∀ x, 0 ≤ f x) {t : ℝ} (hS : 0 < mUPM 0 t (densityLaw f)) :
    IntegrableOn (fun s => f s / mUPM 0 s (densityLaw f)) (Iic t) := by
  have hint := integrable_density hf hf0
  refine Integrable.mono' (hint.div_const (mUPM 0 t (densityLaw f))).integrableOn
    (hf.div (continuous_mUPM_zero_densityLaw hf hf0).measurable).aestronglyMeasurable ?_
  filter_upwards [ae_restrict_mem measurableSet_Iic] with s hs
  have hSs : mUPM 0 t (densityLaw f) ≤ mUPM 0 s (densityLaw f) := mUPM_zero_antitone hs
  rw [Real.norm_of_nonneg (div_nonneg (hf0 s) (hS.trans_le hSs).le)]
  exact div_le_div_of_nonneg_left (hf0 s) hS hSs

lemma mUPM_zero_zero_of_nonneg [IsProbabilityMeasure (densityLaw f)]
    (h0 : densityLaw f (Iio 0) = 0) : mUPM 0 0 (densityLaw f) = 1 := by
  have hsing : densityLaw f {0} = 0 :=
    withDensity_absolutelyContinuous _ _ Real.volume_singleton
  have hIic : densityLaw f (Iic 0) = 0 := by
    rw [← Iio_insert, insert_eq]
    exact measure_union_null hsing h0
  rw [mUPM_zero_eq_one_sub_cdf, cdf_eq_real, measureReal_def, hIic]
  simp

lemma lipschitzOnWith_neg_log {m : ℝ} (hm : 0 < m) :
    LipschitzOnWith (Real.toNNReal m⁻¹) (fun x => -Real.log x) (Ici m) := by
  refine (convex_Ici m).lipschitzOnWith_of_nnnorm_hasDerivWithin_le (f' := fun x => -x⁻¹)
    (fun x hx => ((Real.hasDerivAt_log (hm.trans_le hx).ne').neg).hasDerivWithinAt)
    fun x hx => ?_
  have hx0 : 0 < x := hm.trans_le hx
  rw [← NNReal.coe_le_coe, coe_nnnorm, Real.norm_eq_abs, abs_neg, abs_of_pos (inv_pos.mpr hx0),
    Real.coe_toNNReal _ (inv_pos.mpr hm).le]
  exact inv_anti₀ hm hx

lemma ae_hasDerivAt_cdf_densityLaw [IsProbabilityMeasure (densityLaw f)] (hf : Measurable f)
    (hf0 : ∀ x, 0 ≤ f x) : ∀ᵐ t, HasDerivAt (fun s => cdf (densityLaw f) s) (f t) t := by
  filter_upwards [LocallyIntegrable.ae_hasDerivAt_integral
    (integrable_density hf hf0).locallyIntegrable] with t ht
  simp_rw [cdf_densityLaw_eq_add hf hf0 0]
  exact (ht 0).const_add _

theorem ae_hasDerivAt_cumHazard [IsProbabilityMeasure (densityLaw f)] (hf : Measurable f)
    (hf0 : ∀ x, 0 ≤ f x) :
    ∀ᵐ t, 0 < mUPM 0 t (densityLaw f) → HasDerivAt (fun s => -Real.log (mUPM 0 s (densityLaw f)))
      (f t / mUPM 0 t (densityLaw f)) t := by
  filter_upwards [ae_hasDerivAt_cdf_densityLaw hf hf0] with t ht hS
  exact hasDerivAt_neg_log_mUPM_zero ht hS

lemma absolutelyContinuousOnInterval_cumHazard [IsProbabilityMeasure (densityLaw f)]
    (hf : Measurable f) (hf0 : ∀ x, 0 ≤ f x) {a t : ℝ} (hat : a ≤ t)
    (hS : 0 < mUPM 0 t (densityLaw f)) :
    AbsolutelyContinuousOnInterval (fun s => -Real.log (mUPM 0 s (densityLaw f))) a t := by
  have hI : AbsolutelyContinuousOnInterval (fun s => ∫ v in a..s, f v) a t :=
    (integrable_density hf hf0).intervalIntegrable.absolutelyContinuousOnInterval_intervalIntegral
      left_mem_uIcc
  have hC : AbsolutelyContinuousOnInterval (fun _ : ℝ => 1 - ∫ x in Iic a, f x) a t :=
    (LipschitzWith.const (1 - ∫ x in Iic a, f x)).lipschitzOnWith.absolutelyContinuousOnInterval
  have hSac : AbsolutelyContinuousOnInterval (fun s => mUPM 0 s (densityLaw f)) a t :=
    (hC.sub hI).congr fun s _ => by
      simp only [Pi.sub_apply, mUPM_zero_eq_one_sub_cdf, cdf_densityLaw_eq_add hf hf0 a s]
      ring
  refine (lipschitzOnWith_neg_log hS).comp_absolutelyContinuousOnInterval (fun s hs => ?_) hSac
  rw [uIcc_of_le hat] at hs
  exact mUPM_zero_antitone hs.2

theorem cumHazard_sub_eq_intervalIntegral [IsProbabilityMeasure (densityLaw f)]
    (hf : Measurable f) (hf0 : ∀ x, 0 ≤ f x) {a t : ℝ} (hat : a ≤ t)
    (hS : 0 < mUPM 0 t (densityLaw f)) :
    -Real.log (mUPM 0 t (densityLaw f)) - -Real.log (mUPM 0 a (densityLaw f)) =
      ∫ s in a..t, f s / mUPM 0 s (densityLaw f) := by
  rw [← (absolutelyContinuousOnInterval_cumHazard hf hf0 hat hS).integral_deriv_eq_sub]
  refine intervalIntegral.integral_congr_ae ?_
  filter_upwards [ae_hasDerivAt_cumHazard hf hf0] with s hs hmem
  rw [uIoc_of_le hat] at hmem
  exact (hs (hS.trans_le (mUPM_zero_antitone hmem.2))).deriv

theorem cumHazard_eq_integral_Iic_of_density [IsProbabilityMeasure (densityLaw f)]
    (hf : Measurable f) (hf0 : ∀ x, 0 ≤ f x) {t : ℝ} (hS : 0 < mUPM 0 t (densityLaw f)) :
    -Real.log (mUPM 0 t (densityLaw f)) = ∫ s in Iic t, f s / mUPM 0 s (densityLaw f) := by
  have hlim : Tendsto (fun s => -Real.log (mUPM 0 s (densityLaw f))) atBot (𝓝 0) := by
    have h1 : Tendsto (fun s => mUPM 0 s (densityLaw f)) atBot (𝓝 1) := by
      simp_rw [mUPM_zero_eq_one_sub_cdf]
      simpa using (tendsto_cdf_atBot (densityLaw f)).const_sub 1
    simpa using ((Real.continuousAt_log one_ne_zero).tendsto.comp h1).neg
  have hL : Tendsto (fun a => -Real.log (mUPM 0 t (densityLaw f)) -
      -Real.log (mUPM 0 a (densityLaw f))) atBot
      (𝓝 (-Real.log (mUPM 0 t (densityLaw f)) - 0)) := tendsto_const_nhds.sub hlim
  have hR := intervalIntegral_tendsto_integral_Iic t (integrableOn_hazard hf hf0 hS) tendsto_id
  rw [sub_zero] at hL
  refine tendsto_nhds_unique (hL.congr' ?_) hR
  filter_upwards [eventually_le_atBot t] with a ha
  exact cumHazard_sub_eq_intervalIntegral hf hf0 ha hS

theorem cumHazard_eq_intervalIntegral_of_density [IsProbabilityMeasure (densityLaw f)]
    (hf : Measurable f) (hf0 : ∀ x, 0 ≤ f x) (h0 : densityLaw f (Iio 0) = 0) {t : ℝ}
    (ht : 0 ≤ t) (hS : 0 < mUPM 0 t (densityLaw f)) :
    -Real.log (mUPM 0 t (densityLaw f)) = ∫ s in (0)..t, f s / mUPM 0 s (densityLaw f) := by
  rw [← cumHazard_sub_eq_intervalIntegral hf hf0 ht hS, mUPM_zero_zero_of_nonneg h0]
  simp

end Density

section Exponential

variable {r : ℝ}

lemma expMeasure_Iio_zero (r : ℝ) : expMeasure r (Iio 0) = 0 := by
  rw [expMeasure, gammaMeasure, withDensity_apply _ measurableSet_Iio]
  exact lintegral_exponentialPDF_of_nonpos le_rfl

lemma ae_nonneg_expMeasure (r : ℝ) : ∀ᵐ x ∂expMeasure r, 0 ≤ x := by
  have hset : {x : ℝ | ¬0 ≤ x} = Iio 0 := by ext x; simp
  rw [ae_iff, hset]
  exact expMeasure_Iio_zero r

lemma mUPM_zero_expMeasure (hr : 0 < r) {t : ℝ} (ht : 0 ≤ t) :
    mUPM 0 t (expMeasure r) = Real.exp (-(r * t)) := by
  have := isProbabilityMeasure_expMeasure hr
  rw [mUPM_zero_eq_one_sub_cdf, cdf_expMeasure_eq hr, ite_eq_left ht, sub_sub_cancel]

lemma eUPM_expMeasure (hr : 0 < r) (d : ℕ) {t : ℝ} (ht : 0 ≤ t) :
    eUPM d t (expMeasure r) =
      ENNReal.ofReal (d.factorial * Real.exp (-(r * t)) / r ^ d) := by
  have := isProbabilityMeasure_expMeasure hr
  induction d generalizing t with
  | zero =>
    have hfin : eUPM 0 t (expMeasure r) ≠ ∞ :=
      (eUPM_lt_top (r := 0) (by simpa using continuous_id'.aestronglyMeasurable) t).ne
    rw [← ENNReal.ofReal_toReal hfin, ← mUPM_eq_toReal, mUPM_zero_expMeasure hr ht]
    simp
  | succ d ih =>
    have hI : IntegrableOn (fun s => (d.factorial : ℝ) * Real.exp (-(r * s)) / r ^ d) (Ici t) := by
      rw [integrableOn_Ici_iff_integrableOn_Ioi]
      refine IntegrableOn.congr_fun (((integrableOn_exp_mul_Ioi (neg_neg_of_pos hr) t).const_mul
        (d.factorial : ℝ)).div_const (r ^ d)) (fun s _ => ?_) measurableSet_Ioi
      simp only [neg_mul]
    have hval : ∫ s in Ici t, (d.factorial : ℝ) * Real.exp (-(r * s)) / r ^ d =
        d.factorial * Real.exp (-(r * t)) / r ^ (d + 1) := by
      simp_rw [mul_div_assoc, ← neg_mul]
      rw [integral_Ici_eq_integral_Ioi, integral_const_mul, integral_div,
        integral_exp_mul_Ioi (neg_neg_of_pos hr)]
      field_simp
      ring
    rw [eUPM_succ, setLIntegral_congr_fun measurableSet_Ici fun s hs => ih (ht.trans hs),
      ← ofReal_integral_eq_lintegral_ofReal hI (ae_of_all _ fun s => by positivity), hval,
      ← ENNReal.ofReal_natCast, ← ENNReal.ofReal_one, ← ENNReal.ofReal_add (by positivity)
        zero_le_one, ← ENNReal.ofReal_mul (by positivity)]
    congr 1
    push_cast [Nat.factorial_succ]
    ring

theorem mUPM_expMeasure (hr : 0 < r) (d : ℕ) {t : ℝ} (ht : 0 ≤ t) :
    mUPM d t (expMeasure r) = d.factorial * Real.exp (-(r * t)) / r ^ d := by
  rw [mUPM_eq_toReal, eUPM_expMeasure hr d ht, ENNReal.toReal_ofReal (by positivity)]

lemma integrable_id_expMeasure (hr : 0 < r) : Integrable (fun x : ℝ => x) (expMeasure r) := by
  have h : Integrable (upperDev 1 0) (expMeasure r) :=
    (eUPM_lt_top_iff 1 0).mp (by rw [eUPM_expMeasure hr 1 le_rfl]; exact ENNReal.ofReal_lt_top)
  refine h.congr ?_
  filter_upwards [ae_nonneg_expMeasure r] with x hx
  rcases hx.eq_or_lt with rfl | hx <;> simp [upperDev, *]

theorem integral_id_expMeasure (hr : 0 < r) : ∫ x, x ∂expMeasure r = 1 / r := by
  have := isProbabilityMeasure_expMeasure hr
  have hL : mLPM 1 0 (expMeasure r) = 0 := by
    refine integral_eq_zero_of_ae ?_
    filter_upwards [ae_nonneg_expMeasure r] with x hx
    rcases hx.eq_or_lt with rfl | hx <;> simp [lowerDev, not_le.mpr, *]
  rw [mean_eq_mUPM_one_sub_mLPM_one (integrable_id_expMeasure hr), hL,
    mUPM_expMeasure hr 1 le_rfl]
  simp

theorem cdf_expMeasure_mean (hr : 0 < r) :
    cdf (expMeasure r) (∫ x, x ∂expMeasure r) = 1 - Real.exp (-1) := by
  rw [integral_id_expMeasure hr, cdf_expMeasure_eq hr, ite_eq_left (by positivity),
    mul_one_div_cancel hr.ne']

theorem one_sub_exp_neg_one_bounds :
    0.632 < 1 - Real.exp (-1) ∧ 1 - Real.exp (-1) < 0.633 := by
  have h1 := Real.exp_one_gt_d9
  have h2 := Real.exp_one_lt_d9
  have he : Real.exp (-1) * Real.exp 1 = 1 := by rw [← Real.exp_add]; simp
  have hpos := Real.exp_pos (-1)
  constructor <;> nlinarith

theorem mLPM_one_expMeasure (hr : 0 < r) {t : ℝ} (ht : 0 ≤ t) :
    mLPM 1 t (expMeasure r) = t - 1 / r + Real.exp (-(r * t)) / r := by
  have := isProbabilityMeasure_expMeasure hr
  have h := mUPM_one_sub_mLPM_one (integrable_id_expMeasure hr) t
  rw [integral_id_expMeasure hr, mUPM_expMeasure hr 1 ht] at h
  simp only [Nat.factorial_one, Nat.cast_one, one_mul, pow_one] at h
  linarith

instance : IsProbabilityMeasure (expMeasure 1) := isProbabilityMeasure_expMeasure one_pos

theorem mLPM_one_expMeasure_one {t : ℝ} (ht : 0 ≤ t) :
    mLPM 1 t (expMeasure 1) = t - 1 + Real.exp (-t) := by
  rw [mLPM_one_expMeasure one_pos ht]; simp

theorem mPMRatio_one_expMeasure (hr : 0 < r) {t : ℝ} (ht : 0 ≤ t) :
    mPMRatio 1 t (expMeasure r) =
      (r * t - 1 + Real.exp (-(r * t))) / (r * t - 1 + 2 * Real.exp (-(r * t))) := by
  have hL : t - 1 / r + Real.exp (-(r * t)) / r = (r * t - 1 + Real.exp (-(r * t))) / r := by
    field_simp
  have hS : (r * t - 1 + Real.exp (-(r * t))) / r + Real.exp (-(r * t)) / r =
      (r * t - 1 + 2 * Real.exp (-(r * t))) / r := by ring
  rw [mPMRatio, mLPM_one_expMeasure hr ht, mUPM_expMeasure hr 1 ht]
  simp only [Nat.factorial_one, Nat.cast_one, one_mul, pow_one]
  rw [hL, hS, div_div_div_cancel_right₀ hr.ne']

lemma exp_neg_bounds : Real.exp (-0.29) < 0.749 ∧ 0.74 < Real.exp (-0.30) := by
  have h1 := Real.exp_bound (x := -0.29) (by norm_num [abs_le]) (n := 5) (by norm_num)
  have h2 := Real.exp_bound (x := -0.30) (by norm_num [abs_le]) (n := 5) (by norm_num)
  rw [abs_le] at h1 h2
  norm_num [Finset.sum_range_succ, Nat.factorial] at h1 h2
  constructor <;> linarith [h1.1, h1.2, h2.1, h2.2]

lemma mPMRatio_one_expMeasure_lt (hr : 0 < r) : mPMRatio 1 (0.29 / r) (expMeasure r) < 1 / 20 := by
  have he := exp_neg_bounds.1
  have hlow := Real.add_one_le_exp (-0.29 : ℝ)
  have hu : r * (0.29 / r) = 0.29 := by field_simp
  rw [mPMRatio_one_expMeasure hr (by positivity), hu,
    div_lt_div_iff₀ (by linarith) (by norm_num)]
  linarith

lemma lt_mPMRatio_one_expMeasure (hr : 0 < r) : 1 / 20 < mPMRatio 1 (0.30 / r) (expMeasure r) := by
  have he := exp_neg_bounds.2
  have hu : r * (0.30 / r) = 0.30 := by field_simp
  rw [mPMRatio_one_expMeasure hr (by positivity), hu,
    div_lt_div_iff₀ (by norm_num) (by linarith)]
  linarith

lemma memLp_one_expMeasure (hr : 0 < r) :
    MemLp (fun x : ℝ => x) ((1 : ℕ) : ℝ≥0∞) (expMeasure r) :=
  memLp_one_of_integrable (integrable_id_expMeasure hr)

theorem mPMRatio_one_expMeasure_eq_one_div_twenty (hr : 0 < r) {t : ℝ}
    (h : mPMRatio 1 t (expMeasure r) = 1 / 20) :
    0.29 < r * t ∧ r * t < 0.30 ∧ 0.25 < cdf (expMeasure r) t ∧ cdf (expMeasure r) t < 0.26 := by
  have := isProbabilityMeasure_expMeasure hr
  have hmono := mPMRatio_monotone (memLp_one_expMeasure hr)
  have h1 : 0.29 < r * t := by
    by_contra hle
    have := hmono ((le_div_iff₀ hr).mpr (by linarith : t * r ≤ 0.29))
    linarith [mPMRatio_one_expMeasure_lt hr]
  have h2 : r * t < 0.30 := by
    by_contra hle
    have := hmono ((div_le_iff₀ hr).mpr (by linarith : 0.30 ≤ t * r))
    linarith [lt_mPMRatio_one_expMeasure hr]
  have ht : 0 ≤ t := by
    by_contra hneg
    nlinarith [not_le.mp hneg]
  have he1 : Real.exp (-(r * t)) < Real.exp (-0.29) := Real.exp_lt_exp.mpr (by linarith)
  have he2 : Real.exp (-0.30) < Real.exp (-(r * t)) := Real.exp_lt_exp.mpr (by linarith)
  have hb := exp_neg_bounds
  rw [cdf_expMeasure_eq hr, ite_eq_left ht]
  refine ⟨h1, h2, ?_, ?_⟩ <;> linarith

lemma expMeasure_compl_singleton_ne_zero (hr : 0 < r) (t : ℝ) : expMeasure r {t}ᶜ ≠ 0 := by
  have := isProbabilityMeasure_expMeasure hr
  intro h
  have hsing : expMeasure r {t} = 0 :=
    withDensity_absolutelyContinuous _ _ Real.volume_singleton
  have huniv := measure_union_null hsing h
  rw [union_compl_self, measure_univ] at huniv
  exact one_ne_zero huniv

theorem threshold_one_expMeasure (hr : 0 < r) :
    0.29 < r * GapQuantile.threshold 1 (expMeasure r) (1 / 20) ∧
      r * GapQuantile.threshold 1 (expMeasure r) (1 / 20) < 0.30 ∧
      0.25 < cdf (expMeasure r) (GapQuantile.threshold 1 (expMeasure r) (1 / 20)) ∧
      cdf (expMeasure r) (GapQuantile.threshold 1 (expMeasure r) (1 / 20)) < 0.26 := by
  have := isProbabilityMeasure_expMeasure hr
  exact mPMRatio_one_expMeasure_eq_one_div_twenty hr (GapQuantile.mPMRatio_threshold one_ne_zero
    (memLp_one_expMeasure hr) (expMeasure_compl_singleton_ne_zero hr) (by norm_num) (by norm_num))

theorem threshold_one_two_points {a b : ℝ} (hab : a < b) :
    GapQuantile.threshold 1 (empMeasure ![a, b]) (1 / 20) = (19 * a + b) / 20 ∧
      LPM 0 (GapQuantile.threshold 1 (empMeasure ![a, b]) (1 / 20)) ![a, b] = 1 / 2 := by
  set t := (19 * a + b) / 20 with ht
  have ha : a < t := by rw [ht]; linarith
  have hb : t < b := by rw [ht]; linarith
  have hx : ∃ i j, (![a, b] : Fin 2 → ℝ) i ≠ ![a, b] j := ⟨0, 1, by simpa using hab.ne⟩
  have hnc : ¬ ∀ i : Fin 2, (![a, b] : Fin 2 → ℝ) i = t := fun h => ha.ne (h 0)
  have hr : LPMratio 1 t ![a, b] = 1 / 20 := by
    simp only [LPMratio, LPM, UPM, Fin.sum_univ_two]
    simp [lowerDev, upperDev, hnc, ha.le, hb.not_ge, ha.not_gt, hb]
    field_simp
    rw [ht]
    ring
  have := isProbabilityMeasure_empMeasure (n := 2) (by norm_num) ![a, b]
  have heq : t = GapQuantile.threshold 1 (empMeasure ![a, b]) (1 / 20) :=
    GapQuantile.eq_threshold_of_mPMRatio_eq one_ne_zero
      (memLp_empMeasure (by norm_num) ![a, b] _)
      (GapQuantile.empMeasure_compl_singleton_ne_zero hx) (by norm_num) (by norm_num)
      ((GapQuantile.mPMRatio_empMeasure one_ne_zero hx t).trans hr)
  rw [← heq]
  refine ⟨rfl, ?_⟩
  simp only [LPM, Fin.sum_univ_two]
  simp [lowerDev, ha.le, hb.not_ge]

lemma integrable_lowerDev_expMeasure (hr : 0 < r) (d : ℕ) (t : ℝ) :
    Integrable (lowerDev d t) (expMeasure r) := by
  have := isProbabilityMeasure_expMeasure hr
  refine Integrable.mono' (integrable_const (|t| ^ d))
    (lowerDev_measurable d t).aestronglyMeasurable ?_
  filter_upwards [ae_nonneg_expMeasure r] with x hx
  rw [Real.norm_of_nonneg (lowerDev_nonneg d t x), lowerDev]
  split_ifs with hxt
  · exact pow_le_pow_left₀ (sub_nonneg.mpr hxt) (by rw [abs_of_nonneg (hx.trans hxt)]; linarith) d
  · positivity

lemma mLPM_expMeasure_of_neg (r : ℝ) (d : ℕ) {t : ℝ} (ht : t < 0) :
    mLPM d t (expMeasure r) = 0 := by
  refine integral_eq_zero_of_ae ?_
  filter_upwards [ae_nonneg_expMeasure r] with x hx
  simp [lowerDev, not_le.mpr (ht.trans_le hx)]

lemma setIntegral_Iic_eq_intervalIntegral {φ g : ℝ → ℝ} (h0 : ∀ s < 0, φ s = 0)
    (hg : ∀ s, 0 ≤ s → φ s = g s) {t : ℝ} (ht : 0 ≤ t) :
    ∫ s in Iic t, φ s = ∫ s in (0)..t, g s := by
  have hcongr : ∀ s ∈ Iic t, φ s = (Icc 0 t).indicator g s := fun s hs => by
    by_cases h : 0 ≤ s
    · have hmem : s ∈ Icc 0 t := ⟨h, hs⟩
      rw [indicator_of_mem hmem, hg s h]
    · rw [indicator_of_notMem (fun h' => h h'.1), h0 s (not_le.mp h)]
  rw [setIntegral_congr_fun measurableSet_Iic hcongr, setIntegral_indicator measurableSet_Icc,
    inter_eq_self_of_subset_right Icc_subset_Iic_self, integral_Icc_eq_integral_Ioc,
    intervalIntegral.integral_of_le ht]

theorem mLPM_two_expMeasure_one {t : ℝ} (ht : 0 ≤ t) :
    mLPM 2 t (expMeasure 1) = t ^ 2 - 2 * t + 2 - 2 * Real.exp (-t) := by
  have hderiv : ∀ s ∈ uIcc 0 t, HasDerivAt (fun s : ℝ => s ^ 2 / 2 - s - Real.exp (-s))
      (s - 1 + Real.exp (-s)) s := fun s _ => by
    exact ((((hasDerivAt_pow 2 s).div_const 2).sub (hasDerivAt_id s)).sub
      (hasDerivAt_neg s).exp).congr_deriv (by simp)
  change mLPM (1 + 1) t (expMeasure 1) = _
  rw [mLPM_succ 1 t (integrable_lowerDev_expMeasure one_pos _ t),
    setIntegral_Iic_eq_intervalIntegral (fun s hs => mLPM_expMeasure_of_neg 1 1 hs)
      (fun s hs => mLPM_one_expMeasure_one hs) ht,
    intervalIntegral.integral_eq_sub_of_hasDerivAt hderiv
      (Continuous.intervalIntegrable (by fun_prop) _ _)]
  simp; ring

theorem mLPM_three_expMeasure_one {t : ℝ} (ht : 0 ≤ t) :
    mLPM 3 t (expMeasure 1) = t ^ 3 - 3 * t ^ 2 + 6 * t - 6 + 6 * Real.exp (-t) := by
  have hderiv : ∀ s ∈ uIcc 0 t,
      HasDerivAt (fun s : ℝ => s ^ 3 / 3 - s ^ 2 + 2 * s + 2 * Real.exp (-s))
        (s ^ 2 - 2 * s + 2 - 2 * Real.exp (-s)) s := fun s _ => by
    exact (((((hasDerivAt_pow 3 s).div_const 3).sub (hasDerivAt_pow 2 s)).add
      ((hasDerivAt_id s).const_mul 2)).add ((hasDerivAt_neg s).exp.const_mul 2)).congr_deriv
      (by simp; ring)
  change mLPM (2 + 1) t (expMeasure 1) = _
  rw [mLPM_succ 2 t (integrable_lowerDev_expMeasure one_pos _ t),
    setIntegral_Iic_eq_intervalIntegral (fun s hs => mLPM_expMeasure_of_neg 1 2 hs)
      (fun s hs => mLPM_two_expMeasure_one hs) ht,
    intervalIntegral.integral_eq_sub_of_hasDerivAt hderiv
      (Continuous.intervalIntegrable (by fun_prop) _ _)]
  simp; ring

end Exponential

section NegExp

lemma mLPM_map_neg (μ : Measure ℝ) (d : ℕ) (t : ℝ) :
    mLPM (d + 1) (-t) (μ.map Neg.neg) = mUPM (d + 1) t μ := by
  rw [mLPM, mUPM, integral_map measurable_neg.aemeasurable
    (lowerDev_measurable _ _).aestronglyMeasurable]
  refine integral_congr_ae (Eventually.of_forall fun x => ?_)
  simp only [lowerDev, upperDev, neg_le_neg_iff]
  rcases lt_trichotomy t x with h | rfl | h
  · simp only [h, h.le, ↓reduceIte]; ring
  · simp
  · simp only [not_le.mpr h, not_lt.mpr h.le, ↓reduceIte]

lemma mUPM_map_neg (μ : Measure ℝ) (d : ℕ) (t : ℝ) :
    mUPM (d + 1) (-t) (μ.map Neg.neg) = mLPM (d + 1) t μ := by
  rw [mLPM, mUPM, integral_map measurable_neg.aemeasurable
    (upperDev_measurable _ _).aestronglyMeasurable]
  refine integral_congr_ae (Eventually.of_forall fun x => ?_)
  simp only [lowerDev, upperDev, neg_lt_neg_iff]
  rcases lt_trichotomy x t with h | rfl | h
  · simp only [h, h.le, ↓reduceIte]; ring
  · simp
  · simp only [not_le.mpr h, not_lt.mpr h.le, ↓reduceIte]

def negExp : Measure ℝ := (expMeasure 1).map Neg.neg

instance : IsProbabilityMeasure negExp := by unfold negExp; infer_instance

theorem integral_id_negExp : ∫ x, x ∂negExp = -1 := by
  rw [negExp, integral_map (f := fun y : ℝ => y) measurable_neg.aemeasurable
    aestronglyMeasurable_id, integral_neg, integral_id_expMeasure one_pos, div_one]

theorem mLPM_negExp (d : ℕ) :
    mLPM (d + 1) (-1) negExp = (d + 1).factorial / Real.exp 1 := by
  rw [negExp, mLPM_map_neg, mUPM_expMeasure one_pos _ zero_le_one]
  simp [Real.exp_neg, div_eq_mul_inv]

theorem mUPM_negExp_one : mUPM 1 (-1) negExp = 1 / Real.exp 1 := by
  rw [negExp, mUPM_map_neg _ 0, mLPM_one_expMeasure_one zero_le_one]
  simp [Real.exp_neg]

theorem mUPM_negExp_two : mUPM 2 (-1) negExp = 1 - 2 / Real.exp 1 := by
  rw [negExp, mUPM_map_neg _ 1, mLPM_two_expMeasure_one zero_le_one, Real.exp_neg]
  ring

theorem mUPM_negExp_three : mUPM 3 (-1) negExp = 6 / Real.exp 1 - 2 := by
  rw [negExp, mUPM_map_neg _ 2, mLPM_three_expMeasure_one zero_le_one, Real.exp_neg]
  ring

def negExpProfile (d : ℕ) : ℝ :=
  mCoLPM d d (∫ x, x ∂negExp) (∫ x, x ∂negExp) (negExp.prod negExp) /
    mCoUPM d d (∫ x, x ∂negExp) (∫ x, x ∂negExp) (negExp.prod negExp)

lemma negExpProfile_eq (d : ℕ) :
    negExpProfile d = (mLPM d (-1) negExp / mUPM d (-1) negExp) ^ 2 := by
  rw [negExpProfile, integral_id_negExp, mCoLPM_prod, mCoUPM_prod, div_pow, sq, sq]
  rfl

theorem negExpProfile_one : negExpProfile 1 = 1 := by
  rw [negExpProfile_eq, mLPM_negExp 0, mUPM_negExp_one]
  simp [(Real.exp_pos 1).ne']

theorem negExpProfile_two : negExpProfile 2 = (2 / (Real.exp 1 - 2)) ^ 2 := by
  have he := Real.exp_one_gt_d9
  rw [negExpProfile_eq, mLPM_negExp 1, mUPM_negExp_two]
  congr 1
  have h2 : Real.exp 1 - 2 ≠ 0 := by linarith
  field_simp
  norm_num [Nat.factorial]

theorem negExpProfile_three : negExpProfile 3 = (3 / (3 - Real.exp 1)) ^ 2 := by
  have he := Real.exp_one_lt_d9
  rw [negExpProfile_eq, mLPM_negExp 2, mUPM_negExp_three]
  congr 1
  have h3 : 3 - Real.exp 1 ≠ 0 := by linarith
  have key : 6 / Real.exp 1 - 2 = 2 * (3 - Real.exp 1) / Real.exp 1 := by
    field_simp
    ring
  have hfac : (((2 + 1).factorial : ℕ) : ℝ) = 6 := by norm_num [Nat.factorial]
  rw [key, hfac]
  field_simp
  norm_num

theorem negExpProfile_bounds :
    negExpProfile 1 = 1 ∧ 7.75 < negExpProfile 2 ∧ negExpProfile 2 < 7.76 ∧
      113.3 < negExpProfile 3 ∧ negExpProfile 3 < 113.5 := by
  have h1 := Real.exp_one_gt_d9
  have h2 := Real.exp_one_lt_d9
  have hpos2 : 0 < (Real.exp 1 - 2) ^ 2 := by nlinarith
  have hpos3 : 0 < (3 - Real.exp 1) ^ 2 := by nlinarith
  rw [negExpProfile_two, negExpProfile_three, div_pow, div_pow, lt_div_iff₀ hpos2,
    div_lt_iff₀ hpos2, lt_div_iff₀ hpos3, div_lt_iff₀ hpos3]
  refine ⟨negExpProfile_one, ?_, ?_, ?_, ?_⟩ <;> nlinarith

end NegExp

end

end NNS.GapDensity
