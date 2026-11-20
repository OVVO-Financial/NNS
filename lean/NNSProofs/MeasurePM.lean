-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.PartialMoments
public import Mathlib.Probability.CDF
public import Mathlib.Probability.Moments.Variance
public import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals

namespace NNS

open MeasureTheory ProbabilityTheory Set Filter Topology
open scoped ENNReal

@[expose] public noncomputable section

def mLPM (r : ℕ) (t : ℝ) (μ : Measure ℝ) : ℝ := ∫ x, lowerDev r t x ∂μ

def mUPM (r : ℕ) (t : ℝ) (μ : Measure ℝ) : ℝ := ∫ x, upperDev r t x ∂μ

def eLPM (r : ℕ) (t : ℝ) (μ : Measure ℝ) : ℝ≥0∞ := ∫⁻ x, ENNReal.ofReal (lowerDev r t x) ∂μ

def eUPM (r : ℕ) (t : ℝ) (μ : Measure ℝ) : ℝ≥0∞ := ∫⁻ x, ENNReal.ofReal (upperDev r t x) ∂μ

variable {μ : Measure ℝ}

lemma lowerDev_eq_indicator (r : ℕ) (t : ℝ) :
    lowerDev r t = (Iic t).indicator fun v => (t - v) ^ r := by
  ext v; simp [lowerDev, indicator]

lemma upperDev_eq_indicator (r : ℕ) (t : ℝ) :
    upperDev r t = (Ioi t).indicator fun v => (v - t) ^ r := by
  ext v; simp [upperDev, indicator]

@[fun_prop]
lemma lowerDev_measurable (r : ℕ) (t : ℝ) : Measurable (lowerDev r t) := by
  rw [lowerDev_eq_indicator]; exact (by fun_prop : Measurable _).indicator measurableSet_Iic

@[fun_prop]
lemma upperDev_measurable (r : ℕ) (t : ℝ) : Measurable (upperDev r t) := by
  rw [upperDev_eq_indicator]; exact (by fun_prop : Measurable _).indicator measurableSet_Ioi

lemma lowerDev_le_abs_pow (r : ℕ) (t x : ℝ) : lowerDev r t x ≤ |x - t| ^ r := by
  unfold lowerDev
  split_ifs with h
  · rw [abs_sub_comm, abs_of_nonneg (sub_nonneg.mpr h)]
  · positivity

lemma upperDev_le_abs_pow (r : ℕ) (t x : ℝ) : upperDev r t x ≤ |x - t| ^ r := by
  unfold upperDev
  split_ifs with h
  · rw [abs_of_nonneg (sub_nonneg.mpr h.le)]
  · positivity

lemma lowerDev_mono (r : ℕ) (x : ℝ) : Monotone fun t => lowerDev r t x := by
  intro s t hst
  simp only [lowerDev]
  split_ifs with h1 h2
  · exact pow_le_pow_left₀ (sub_nonneg.mpr h1) (by linarith) r
  · exact absurd (h1.trans hst) h2
  · exact pow_nonneg (sub_nonneg.mpr ‹x ≤ t›) r
  · exact le_rfl

lemma upperDev_anti (r : ℕ) (x : ℝ) : Antitone fun t => upperDev r t x := by
  intro s t hst
  simp only [upperDev]
  by_cases h1 : t < x
  · simp only [h1, hst.trans_lt h1, ↓reduceIte]
    exact pow_le_pow_left₀ (sub_nonneg.mpr h1.le) (by linarith) r
  · simp only [h1, ↓reduceIte]
    split_ifs with h2
    · exact pow_nonneg (sub_nonneg.mpr h2.le) r
    · exact le_rfl

lemma integrable_abs_sub_pow [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ)
    (t : ℝ) : Integrable (fun x : ℝ => |x - t| ^ r) μ :=
  (h.sub (memLp_const t)).integrable_norm_pow'

lemma integrable_lowerDev [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ)
    (t : ℝ) : Integrable (lowerDev r t) μ :=
  (integrable_abs_sub_pow h t).mono' (lowerDev_measurable r t).aestronglyMeasurable
    (Eventually.of_forall fun x => by
      rw [Real.norm_of_nonneg (lowerDev_nonneg r t x)]; exact lowerDev_le_abs_pow r t x)

lemma integrable_upperDev [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ)
    (t : ℝ) : Integrable (upperDev r t) μ :=
  (integrable_abs_sub_pow h t).mono' (upperDev_measurable r t).aestronglyMeasurable
    (Eventually.of_forall fun x => by
      rw [Real.norm_of_nonneg (upperDev_nonneg r t x)]; exact upperDev_le_abs_pow r t x)

lemma integrable_sub_pow [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ)
    (t : ℝ) : Integrable (fun x : ℝ => (x - t) ^ r) μ :=
  (integrable_abs_sub_pow h t).mono' (by fun_prop)
    (Eventually.of_forall fun x => by rw [Real.norm_eq_abs, abs_pow])

theorem mLPM_zero_eq_cdf [IsProbabilityMeasure μ] (t : ℝ) : mLPM 0 t μ = cdf μ t := by
  rw [cdf_eq_real, mLPM, lowerDev_eq_indicator]
  simp only [pow_zero]
  exact integral_indicator_one measurableSet_Iic

theorem mUPM_zero_eq [IsProbabilityMeasure μ] (t : ℝ) : mUPM 0 t μ = μ.real (Ioi t) := by
  rw [mUPM, upperDev_eq_indicator]
  simp only [pow_zero]
  exact integral_indicator_one measurableSet_Ioi

theorem mLPM_eq_setIntegral (r : ℕ) (t : ℝ) (μ : Measure ℝ) :
    mLPM r t μ = ∫ x in Iic t, (t - x) ^ r ∂μ := by
  rw [mLPM, lowerDev_eq_indicator, integral_indicator measurableSet_Iic]

theorem mUPM_eq_setIntegral (r : ℕ) (t : ℝ) (μ : Measure ℝ) :
    mUPM r t μ = ∫ x in Ioi t, (x - t) ^ r ∂μ := by
  rw [mUPM, upperDev_eq_indicator, integral_indicator measurableSet_Ioi]

theorem mLPM_nonneg (r : ℕ) (t : ℝ) (μ : Measure ℝ) : 0 ≤ mLPM r t μ :=
  integral_nonneg (lowerDev_nonneg r t)

theorem mUPM_nonneg (r : ℕ) (t : ℝ) (μ : Measure ℝ) : 0 ≤ mUPM r t μ :=
  integral_nonneg (upperDev_nonneg r t)

theorem mUPM_add_neg_one_pow_mul_mLPM [IsFiniteMeasure μ] {r : ℕ}
    (h : MemLp (fun x : ℝ => x) r μ) (t : ℝ) :
    mUPM r t μ + (-1) ^ r * mLPM r t μ = ∫ x, (x - t) ^ r ∂μ := by
  rw [mUPM, mLPM, ← integral_const_mul, ← integral_add (integrable_upperDev h t)
    ((integrable_lowerDev h t).const_mul _)]
  exact integral_congr_ae (Eventually.of_forall fun x =>
    upperDev_add_neg_one_pow_mul_lowerDev r t x)

theorem mLPM_zero_add_mUPM_zero [IsProbabilityMeasure μ] (t : ℝ) :
    mLPM 0 t μ + mUPM 0 t μ = 1 := by
  have := mUPM_add_neg_one_pow_mul_mLPM (μ := μ) (r := 0)
    (by rw [Nat.cast_zero, memLp_zero_iff_aestronglyMeasurable]; exact aestronglyMeasurable_id) t
  simp only [pow_zero, one_mul, integral_const, probReal_univ, smul_eq_mul, mul_one] at this
  linarith

lemma upperDev_add_lowerDev (r : ℕ) (t x : ℝ) :
    upperDev r t x + lowerDev r t x = |x - t| ^ r := by
  unfold upperDev lowerDev
  by_cases h : t < x
  · simp [h, not_le.mpr h, abs_of_pos (sub_pos.mpr h)]
  · have h' : x ≤ t := not_lt.mp h
    simp [h, h', abs_sub_comm x t, abs_of_nonneg (sub_nonneg.mpr h')]

lemma max_sub_max_neg (x : ℝ) : max x 0 - max (-x) 0 = x := by
  rcases le_total x 0 with h | h
  · rw [max_eq_right h, max_eq_left (neg_nonneg.mpr h)]; ring
  · rw [max_eq_left h, max_eq_right (neg_nonpos.mpr h)]; ring

lemma max_add_max_neg (x : ℝ) : max x 0 + max (-x) 0 = |x| := by
  rcases le_total x 0 with h | h
  · rw [max_eq_right h, max_eq_left (neg_nonneg.mpr h), abs_of_nonpos h]; ring
  · rw [max_eq_left h, max_eq_right (neg_nonpos.mpr h), abs_of_nonneg h]; ring

lemma upperDev_one (t x : ℝ) : upperDev 1 t x = max (x - t) 0 := by
  unfold upperDev
  split_ifs with h
  · rw [pow_one, max_eq_left (sub_nonneg.mpr h.le)]
  · rw [max_eq_right (sub_nonpos.mpr (not_lt.mp h))]

lemma lowerDev_one (t x : ℝ) : lowerDev 1 t x = max (t - x) 0 := by
  unfold lowerDev
  split_ifs with h
  · rw [pow_one, max_eq_left (sub_nonneg.mpr h)]
  · rw [max_eq_right (sub_nonpos.mpr (not_le.mp h).le)]

lemma lowerDev_eq_max_pow {r : ℕ} (hr : r ≠ 0) (t x : ℝ) :
    lowerDev r t x = max (t - x) 0 ^ r := by
  rw [← lowerDev_one]; unfold lowerDev
  split_ifs <;> simp [zero_pow hr]

lemma upperDev_eq_max_pow {r : ℕ} (hr : r ≠ 0) (t x : ℝ) :
    upperDev r t x = max (x - t) 0 ^ r := by
  rw [← upperDev_one]; unfold upperDev
  split_ifs <;> simp [zero_pow hr]

theorem mUPM_add_mLPM [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ) (t : ℝ) :
    mUPM r t μ + mLPM r t μ = ∫ x, |x - t| ^ r ∂μ := by
  rw [mUPM, mLPM, ← integral_add (integrable_upperDev h t) (integrable_lowerDev h t)]
  exact integral_congr_ae (Eventually.of_forall fun x => upperDev_add_lowerDev r t x)

lemma memLp_one_of_integrable (h : Integrable (fun x : ℝ => x) μ) :
    MemLp (fun x : ℝ => x) ((1 : ℕ) : ℝ≥0∞) μ := by
  rw [Nat.cast_one, memLp_one_iff_integrable]; exact h

theorem mUPM_one_sub_mLPM_one [IsProbabilityMeasure μ] (h : Integrable (fun x : ℝ => x) μ)
    (t : ℝ) : mUPM 1 t μ - mLPM 1 t μ = (∫ x, x ∂μ) - t := by
  have := mUPM_add_neg_one_pow_mul_mLPM (memLp_one_of_integrable h) t
  simp_rw [pow_one] at this
  rw [integral_sub h (integrable_const t), integral_const, probReal_univ,
    one_smul] at this
  linarith

theorem mean_eq_mUPM_one_sub_mLPM_one [IsProbabilityMeasure μ]
    (h : Integrable (fun x : ℝ => x) μ) : ∫ x, x ∂μ = mUPM 1 0 μ - mLPM 1 0 μ := by
  rw [mUPM_one_sub_mLPM_one h, sub_zero]

theorem mUPM_one_add_mLPM_one [IsFiniteMeasure μ] (h : Integrable (fun x : ℝ => x) μ)
    (t : ℝ) : mUPM 1 t μ + mLPM 1 t μ = ∫ x, |x - t| ∂μ := by
  rw [mUPM_add_mLPM (memLp_one_of_integrable h)]; simp

theorem mUPM_two_add_mLPM_two [IsProbabilityMeasure μ] (h : MemLp (fun x : ℝ => x) 2 μ)
    (t : ℝ) :
    mUPM 2 t μ + mLPM 2 t μ = variance (fun x : ℝ => x) μ + ((∫ x, x ∂μ) - t) ^ 2 := by
  have h2 : MemLp (fun x : ℝ => x) ((2 : ℕ) : ℝ≥0∞) μ := by exact_mod_cast h
  have := mUPM_add_neg_one_pow_mul_mLPM h2 t
  rw [neg_one_sq, one_mul] at this
  have hY : MemLp (fun x : ℝ => x - t) 2 μ := h.sub (memLp_const t)
  rw [this, ← variance_sub_const (X := fun x : ℝ => x) aestronglyMeasurable_id t,
    variance_eq_sub hY]
  have hi : Integrable (fun x : ℝ => x) μ := h.integrable one_le_two
  simp only [Pi.pow_apply]
  rw [integral_sub hi (integrable_const t), integral_const, probReal_univ, one_smul]
  ring

theorem variance_eq_mUPM_two_add_mLPM_two [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 2 μ) :
    variance (fun x : ℝ => x) μ = mUPM 2 (∫ x, x ∂μ) μ + mLPM 2 (∫ x, x ∂μ) μ := by
  rw [mUPM_two_add_mLPM_two h, sub_self]; ring

theorem thirdCentralMoment_eq_mUPM_three_sub_mLPM_three [IsProbabilityMeasure μ]
    (h : MemLp (fun x : ℝ => x) 3 μ) :
    ∫ x, (x - ∫ y, y ∂μ) ^ 3 ∂μ = mUPM 3 (∫ y, y ∂μ) μ - mLPM 3 (∫ y, y ∂μ) μ := by
  have h3 : MemLp (fun x : ℝ => x) ((3 : ℕ) : ℝ≥0∞) μ := by exact_mod_cast h
  have := mUPM_add_neg_one_pow_mul_mLPM h3 (∫ y, y ∂μ)
  norm_num at this
  linarith

theorem mLPM_eq_toReal (r : ℕ) (t : ℝ) (μ : Measure ℝ) : mLPM r t μ = (eLPM r t μ).toReal :=
  integral_eq_lintegral_of_nonneg_ae (Eventually.of_forall (lowerDev_nonneg r t))
    (lowerDev_measurable r t).aestronglyMeasurable

theorem mUPM_eq_toReal (r : ℕ) (t : ℝ) (μ : Measure ℝ) : mUPM r t μ = (eUPM r t μ).toReal :=
  integral_eq_lintegral_of_nonneg_ae (Eventually.of_forall (upperDev_nonneg r t))
    (upperDev_measurable r t).aestronglyMeasurable

lemma eLPM_lt_top_iff (r : ℕ) (t : ℝ) : eLPM r t μ < ∞ ↔ Integrable (lowerDev r t) μ := by
  rw [Integrable, and_iff_right (lowerDev_measurable r t).aestronglyMeasurable,
    hasFiniteIntegral_iff_ofReal (Eventually.of_forall (lowerDev_nonneg r t)), eLPM]

lemma eUPM_lt_top_iff (r : ℕ) (t : ℝ) : eUPM r t μ < ∞ ↔ Integrable (upperDev r t) μ := by
  rw [Integrable, and_iff_right (upperDev_measurable r t).aestronglyMeasurable,
    hasFiniteIntegral_iff_ofReal (Eventually.of_forall (upperDev_nonneg r t)), eUPM]

theorem eLPM_lt_top [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ) (t : ℝ) :
    eLPM r t μ < ∞ := (eLPM_lt_top_iff r t).mpr (integrable_lowerDev h t)

theorem eUPM_lt_top [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ) (t : ℝ) :
    eUPM r t μ < ∞ := (eUPM_lt_top_iff r t).mpr (integrable_upperDev h t)

def empMeasure {n : ℕ} (x : Fin n → ℝ) : Measure ℝ :=
  (n : ℝ≥0∞)⁻¹ • ∑ i, Measure.dirac (x i)

theorem integral_empMeasure {n : ℕ} (x : Fin n → ℝ) (f : ℝ → ℝ) :
    ∫ v, f v ∂empMeasure x = (∑ i, f (x i)) / n := by
  rw [empMeasure, integral_smul_measure,
    integral_finsetSum_measure fun i _ => integrable_dirac enorm_lt_top]
  simp [integral_dirac, div_eq_inv_mul]

theorem mLPM_empMeasure {n : ℕ} (d : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    mLPM d t (empMeasure x) = LPM d t x := integral_empMeasure x _

theorem mUPM_empMeasure {n : ℕ} (d : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    mUPM d t (empMeasure x) = UPM d t x := integral_empMeasure x _

theorem isProbabilityMeasure_empMeasure {n : ℕ} (hn : 0 < n) (x : Fin n → ℝ) :
    IsProbabilityMeasure (empMeasure x) := by
  constructor
  simp only [empMeasure, Measure.smul_apply, Measure.coe_finsetSum, Finset.sum_apply,
    measure_univ, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one,
    smul_eq_mul]
  exact ENNReal.inv_mul_cancel (by exact_mod_cast hn.ne') (ENNReal.natCast_ne_top n)

theorem cdf_empMeasure {n : ℕ} (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) :
    haveI := isProbabilityMeasure_empMeasure hn x
    cdf (empMeasure x) t = LPM 0 t x := by
  have := isProbabilityMeasure_empMeasure hn x
  rw [← mLPM_zero_eq_cdf, mLPM_empMeasure]

theorem memLp_empMeasure {n : ℕ} (hn : 0 < n) (x : Fin n → ℝ) (p : ℝ≥0∞) :
    MemLp (fun v : ℝ => v) p (empMeasure x) := by
  refine MemLp.smul_measure ?_ (ENNReal.inv_ne_top.mpr (by exact_mod_cast hn.ne'))
  exact Finset.sum_induction _ (fun ν => MemLp (fun v : ℝ => v) p ν)
    (fun _ _ h1 h2 => h1.add_measure h2) memLp_measure_zero
    fun i _ => memLp_dirac enorm_lt_top

lemma lowerDev_target_eq_indicator (r : ℕ) (x : ℝ) :
    (fun s => lowerDev r s x) = (Ici x).indicator fun s => (s - x) ^ r := by
  ext s; simp [lowerDev, indicator]

lemma measurable_lowerDev_uncurry (r : ℕ) : Measurable fun p : ℝ × ℝ => lowerDev r p.1 p.2 := by
  unfold lowerDev
  exact Measurable.ite (measurableSet_le measurable_snd measurable_fst) (by fun_prop)
    measurable_const

lemma integrableOn_lowerDev_target (r : ℕ) (t x : ℝ) :
    IntegrableOn (fun s => lowerDev r s x) (Iic t) := by
  rw [lowerDev_target_eq_indicator, IntegrableOn,
    integrable_indicator_iff measurableSet_Ici, IntegrableOn,
    Measure.restrict_restrict measurableSet_Ici, Ici_inter_Iic]
  exact (by fun_prop : Continuous fun s : ℝ => (s - x) ^ r).integrableOn_Icc

lemma lowerDev_succ_eq_setIntegral_Iic (r : ℕ) (t x : ℝ) :
    lowerDev (r + 1) t x = ((r : ℝ) + 1) * ∫ s in Iic t, lowerDev r s x := by
  rw [lowerDev_target_eq_indicator, setIntegral_indicator measurableSet_Ici, inter_comm,
    Ici_inter_Iic]
  by_cases hxt : x ≤ t
  · rw [integral_Icc_eq_integral_Ioc, ← intervalIntegral.integral_of_le hxt,
      intervalIntegral.integral_comp_sub_right (fun s => s ^ r), integral_pow, sub_self,
      zero_pow (Nat.succ_ne_zero r), sub_zero, lowerDev]
    simp only [hxt, ↓reduceIte]
    field_simp
  · rw [Icc_eq_empty hxt, Measure.restrict_empty, integral_zero_measure, mul_zero, lowerDev]
    simp only [hxt, ↓reduceIte]

lemma ofReal_lowerDev_succ (r : ℕ) (t x : ℝ) :
    ENNReal.ofReal (lowerDev (r + 1) t x) =
      ((r : ℝ≥0∞) + 1) * ∫⁻ s in Iic t, ENNReal.ofReal (lowerDev r s x) := by
  rw [lowerDev_succ_eq_setIntegral_Iic, ENNReal.ofReal_mul (by positivity),
    ofReal_integral_eq_lintegral_ofReal (integrableOn_lowerDev_target r t x)
      (Eventually.of_forall fun s => lowerDev_nonneg r s x),
    ENNReal.ofReal_add (by positivity) zero_le_one, ENNReal.ofReal_natCast, ENNReal.ofReal_one]

theorem eLPM_succ (r : ℕ) (t : ℝ) (μ : Measure ℝ) [SFinite μ] :
    eLPM (r + 1) t μ = ((r : ℝ≥0∞) + 1) * ∫⁻ s in Iic t, eLPM r s μ := by
  simp_rw [eLPM, ofReal_lowerDev_succ]
  rw [lintegral_const_mul' _ _ (by simp)]
  congr 1
  refine lintegral_lintegral_swap ?_
  exact ((measurable_lowerDev_uncurry r).comp measurable_swap).ennreal_ofReal.aemeasurable

lemma measurable_eLPM_target (r : ℕ) (μ : Measure ℝ) [SFinite μ] :
    Measurable fun s => eLPM r s μ :=
  (measurable_lowerDev_uncurry r).ennreal_ofReal.lintegral_prod_right'

theorem mLPM_succ [SFinite μ] (r : ℕ) (t : ℝ) (h : Integrable (lowerDev (r + 1) t) μ) :
    mLPM (r + 1) t μ = ((r : ℝ) + 1) * ∫ s in Iic t, mLPM r s μ := by
  have hfin := (eLPM_lt_top_iff (μ := μ) (r + 1) t).mpr h
  rw [eLPM_succ] at hfin
  have hI : ∫⁻ s in Iic t, eLPM r s μ ≠ ∞ := by
    intro hI
    rw [hI, ENNReal.mul_top (by simp)] at hfin
    exact lt_irrefl _ hfin
  simp_rw [mLPM_eq_toReal]
  rw [eLPM_succ, ENNReal.toReal_mul, integral_toReal (measurable_eLPM_target r μ).aemeasurable
    (ae_lt_top (measurable_eLPM_target r μ) hI)]
  congr 1

theorem mLPM_one_eq_integral_cdf [IsProbabilityMeasure μ] (t : ℝ)
    (h : Integrable (lowerDev 1 t) μ) : mLPM 1 t μ = ∫ s in Iic t, cdf μ s := by
  rw [mLPM_succ 0 t h]
  simp [mLPM_zero_eq_cdf]

lemma lowerDev_le_lowerDev_add_one {k r : ℕ} (hk : k ≤ r) {s t : ℝ} (hst : s ≤ t) (x : ℝ) :
    lowerDev k s x ≤ lowerDev r t x + 1 := by
  refine (lowerDev_mono k x hst).trans ?_
  change (if x ≤ t then (t - x) ^ k else 0) ≤ (if x ≤ t then (t - x) ^ r else 0) + 1
  split_ifs with h
  · rcases le_total (t - x) 1 with h1 | h1
    · have := pow_le_one₀ (n := k) (sub_nonneg.mpr h) h1
      have := pow_nonneg (sub_nonneg.mpr h) r
      linarith
    · have := pow_le_pow_right₀ h1 hk
      linarith
  · norm_num

lemma integrable_lowerDev_of_le [IsFiniteMeasure μ] {k r : ℕ} (hk : k ≤ r) {s t : ℝ}
    (hst : s ≤ t) (h : Integrable (lowerDev r t) μ) : Integrable (lowerDev k s) μ :=
  (h.add (integrable_const 1)).mono' (lowerDev_measurable k s).aestronglyMeasurable
    (Eventually.of_forall fun x => by
      rw [Real.norm_of_nonneg (lowerDev_nonneg k s x)]
      exact lowerDev_le_lowerDev_add_one hk hst x)

theorem mLPM_two_eq_iterated_integral_cdf [IsProbabilityMeasure μ] (t : ℝ)
    (h : Integrable (lowerDev 2 t) μ) :
    mLPM 2 t μ = 2 * ∫ s in Iic t, ∫ u in Iic s, cdf μ u := by
  rw [mLPM_succ 1 t h]
  norm_num
  exact setIntegral_congr_fun measurableSet_Iic fun s hs =>
    mLPM_one_eq_integral_cdf s (integrable_lowerDev_of_le one_le_two hs h)

theorem forall_mLPM_one_le_iff [IsProbabilityMeasure μ] {ν : Measure ℝ} [IsProbabilityMeasure ν]
    (hμ : ∀ t, Integrable (lowerDev 1 t) μ) (hν : ∀ t, Integrable (lowerDev 1 t) ν) :
    (∀ t, mLPM 1 t μ ≤ mLPM 1 t ν) ↔
      ∀ t, ∫ s in Iic t, cdf μ s ≤ ∫ s in Iic t, cdf ν s := by
  have eμ := fun t => mLPM_one_eq_integral_cdf t (hμ t)
  have eν := fun t => mLPM_one_eq_integral_cdf (μ := ν) t (hν t)
  simp_rw [eμ, eν]

theorem forall_mLPM_two_le_iff [IsProbabilityMeasure μ] {ν : Measure ℝ} [IsProbabilityMeasure ν]
    (hμ : ∀ t, Integrable (lowerDev 2 t) μ) (hν : ∀ t, Integrable (lowerDev 2 t) ν) :
    (∀ t, mLPM 2 t μ ≤ mLPM 2 t ν) ↔
      ∀ t, ∫ s in Iic t, ∫ u in Iic s, cdf μ u ≤ ∫ s in Iic t, ∫ u in Iic s, cdf ν u := by
  have eμ := fun t => mLPM_two_eq_iterated_integral_cdf t (hμ t)
  have eν := fun t => mLPM_two_eq_iterated_integral_cdf (μ := ν) t (hν t)
  simp_rw [eμ, eν]
  constructor <;> intro H t <;> linarith [H t]

lemma upperDev_target_eq_indicator (r : ℕ) (x : ℝ) :
    (fun s => upperDev r s x) = (Iio x).indicator fun s => (x - s) ^ r := by
  ext s; simp [upperDev, indicator]

lemma measurable_upperDev_uncurry (r : ℕ) : Measurable fun p : ℝ × ℝ => upperDev r p.1 p.2 := by
  unfold upperDev
  exact Measurable.ite (measurableSet_lt measurable_fst measurable_snd) (by fun_prop)
    measurable_const

lemma integrableOn_upperDev_target (r : ℕ) (t x : ℝ) :
    IntegrableOn (fun s => upperDev r s x) (Ici t) := by
  rw [upperDev_target_eq_indicator, IntegrableOn,
    integrable_indicator_iff measurableSet_Iio, IntegrableOn,
    Measure.restrict_restrict measurableSet_Iio, inter_comm, Ici_inter_Iio]
  exact (by fun_prop : Continuous fun s : ℝ => (x - s) ^ r).integrableOn_Icc.mono_set
    Ico_subset_Icc_self

lemma upperDev_succ_eq_setIntegral_Ici (r : ℕ) (t x : ℝ) :
    upperDev (r + 1) t x = ((r : ℝ) + 1) * ∫ s in Ici t, upperDev r s x := by
  rw [upperDev_target_eq_indicator, setIntegral_indicator measurableSet_Iio, Ici_inter_Iio]
  by_cases htx : t < x
  · rw [integral_Ico_eq_integral_Ioo, ← integral_Ioc_eq_integral_Ioo,
      ← intervalIntegral.integral_of_le htx.le,
      intervalIntegral.integral_comp_sub_left (fun s => s ^ r), integral_pow, sub_self,
      zero_pow (Nat.succ_ne_zero r), sub_zero, upperDev]
    simp only [htx, ↓reduceIte]
    field_simp
  · rw [Ico_eq_empty htx, Measure.restrict_empty, integral_zero_measure, mul_zero, upperDev]
    simp only [htx, ↓reduceIte]

lemma ofReal_upperDev_succ (r : ℕ) (t x : ℝ) :
    ENNReal.ofReal (upperDev (r + 1) t x) =
      ((r : ℝ≥0∞) + 1) * ∫⁻ s in Ici t, ENNReal.ofReal (upperDev r s x) := by
  rw [upperDev_succ_eq_setIntegral_Ici, ENNReal.ofReal_mul (by positivity),
    ofReal_integral_eq_lintegral_ofReal (integrableOn_upperDev_target r t x)
      (Eventually.of_forall fun s => upperDev_nonneg r s x),
    ENNReal.ofReal_add (by positivity) zero_le_one, ENNReal.ofReal_natCast, ENNReal.ofReal_one]

theorem eUPM_succ (r : ℕ) (t : ℝ) (μ : Measure ℝ) [SFinite μ] :
    eUPM (r + 1) t μ = ((r : ℝ≥0∞) + 1) * ∫⁻ s in Ici t, eUPM r s μ := by
  simp_rw [eUPM, ofReal_upperDev_succ]
  rw [lintegral_const_mul' _ _ (by simp)]
  congr 1
  refine lintegral_lintegral_swap ?_
  exact ((measurable_upperDev_uncurry r).comp measurable_swap).ennreal_ofReal.aemeasurable

lemma measurable_eUPM_target (r : ℕ) (μ : Measure ℝ) [SFinite μ] :
    Measurable fun s => eUPM r s μ :=
  (measurable_upperDev_uncurry r).ennreal_ofReal.lintegral_prod_right'

theorem mUPM_succ [SFinite μ] (r : ℕ) (t : ℝ) (h : Integrable (upperDev (r + 1) t) μ) :
    mUPM (r + 1) t μ = ((r : ℝ) + 1) * ∫ s in Ici t, mUPM r s μ := by
  have hfin := (eUPM_lt_top_iff (μ := μ) (r + 1) t).mpr h
  rw [eUPM_succ] at hfin
  have hI : ∫⁻ s in Ici t, eUPM r s μ ≠ ∞ := by
    intro hI
    rw [hI, ENNReal.mul_top (by simp)] at hfin
    exact lt_irrefl _ hfin
  simp_rw [mUPM_eq_toReal]
  rw [eUPM_succ, ENNReal.toReal_mul, integral_toReal (measurable_eUPM_target r μ).aemeasurable
    (ae_lt_top (measurable_eUPM_target r μ) hI)]
  congr 1

theorem mLPM_monotone [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ) :
    Monotone fun t => mLPM r t μ := fun s t hst =>
  integral_mono (integrable_lowerDev h s) (integrable_lowerDev h t) fun x =>
    lowerDev_mono r x hst

theorem mUPM_antitone [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ) :
    Antitone fun t => mUPM r t μ := fun s t hst =>
  integral_mono (integrable_upperDev h t) (integrable_upperDev h s) fun x =>
    upperDev_anti r x hst

lemma integrable_abs_add_pow [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ)
    (c : ℝ) : Integrable (fun x : ℝ => (|x| + |c| + 1) ^ r) μ := by
  have h' : MemLp (fun x : ℝ => ‖x‖ + (|c| + 1)) r μ := h.norm.add (memLp_const _)
  refine h'.integrable_norm_pow'.congr (Eventually.of_forall fun x => ?_)
  simp only [Real.norm_eq_abs]
  rw [abs_of_nonneg (by positivity), add_assoc]

lemma abs_sub_le_of_dist_lt {x t t₀ : ℝ} (ht : dist t t₀ < 1) : |x - t| ≤ |x| + |t₀| + 1 := by
  rw [Real.dist_eq] at ht
  have h1 : |x - t| ≤ |x| + |t| := abs_sub _ _
  have h2 : |t| ≤ |t₀| + |t - t₀| := by
    calc |t| = |t₀ + (t - t₀)| := by ring_nf
      _ ≤ |t₀| + |t - t₀| := abs_add_le _ _
  linarith

theorem continuous_mLPM [IsFiniteMeasure μ] {r : ℕ} (hr : r ≠ 0)
    (h : MemLp (fun x : ℝ => x) r μ) : Continuous fun t => mLPM r t μ := by
  refine continuous_iff_continuousAt.mpr fun t₀ => ?_
  refine continuousAt_of_dominated (bound := fun x => (|x| + |t₀| + 1) ^ r)
    (Eventually.of_forall fun t => (lowerDev_measurable r t).aestronglyMeasurable) ?_
    (integrable_abs_add_pow h t₀) (Eventually.of_forall fun x => ?_)
  · filter_upwards [Metric.ball_mem_nhds t₀ one_pos] with t ht
    refine Eventually.of_forall fun x => ?_
    rw [Real.norm_of_nonneg (lowerDev_nonneg r t x)]
    exact (lowerDev_le_abs_pow r t x).trans
      (pow_le_pow_left₀ (abs_nonneg _) (abs_sub_le_of_dist_lt ht) r)
  · simp_rw [lowerDev_eq_max_pow hr]
    fun_prop

theorem continuous_mUPM [IsFiniteMeasure μ] {r : ℕ} (hr : r ≠ 0)
    (h : MemLp (fun x : ℝ => x) r μ) : Continuous fun t => mUPM r t μ := by
  refine continuous_iff_continuousAt.mpr fun t₀ => ?_
  refine continuousAt_of_dominated (bound := fun x => (|x| + |t₀| + 1) ^ r)
    (Eventually.of_forall fun t => (upperDev_measurable r t).aestronglyMeasurable) ?_
    (integrable_abs_add_pow h t₀) (Eventually.of_forall fun x => ?_)
  · filter_upwards [Metric.ball_mem_nhds t₀ one_pos] with t ht
    refine Eventually.of_forall fun x => ?_
    rw [Real.norm_of_nonneg (upperDev_nonneg r t x)]
    exact (upperDev_le_abs_pow r t x).trans
      (pow_le_pow_left₀ (abs_nonneg _) (abs_sub_le_of_dist_lt ht) r)
  · simp_rw [upperDev_eq_max_pow hr]
    fun_prop

theorem tendsto_mLPM_atBot [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ) :
    Tendsto (fun t => mLPM r t μ) atBot (𝓝 0) := by
  have := tendsto_integral_filter_of_dominated_convergence (μ := μ) (l := atBot)
    (F := fun t x => lowerDev r t x) (f := fun _ => (0 : ℝ)) (lowerDev r 0)
    (Eventually.of_forall fun t => (lowerDev_measurable r t).aestronglyMeasurable) ?_
    (integrable_lowerDev h 0) (Eventually.of_forall fun x => ?_)
  · rw [integral_zero] at this; exact this
  · filter_upwards [eventually_le_atBot 0] with t ht
    refine Eventually.of_forall fun x => ?_
    rw [Real.norm_of_nonneg (lowerDev_nonneg r t x)]
    exact lowerDev_mono r x ht
  · refine tendsto_const_nhds.congr' ?_
    filter_upwards [eventually_lt_atBot x] with t ht
    simp [lowerDev, not_le.mpr ht]

theorem tendsto_mUPM_atTop [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ) :
    Tendsto (fun t => mUPM r t μ) atTop (𝓝 0) := by
  have := tendsto_integral_filter_of_dominated_convergence (μ := μ) (l := atTop)
    (F := fun t x => upperDev r t x) (f := fun _ => (0 : ℝ)) (upperDev r 0)
    (Eventually.of_forall fun t => (upperDev_measurable r t).aestronglyMeasurable) ?_
    (integrable_upperDev h 0) (Eventually.of_forall fun x => ?_)
  · rw [integral_zero] at this; exact this
  · filter_upwards [eventually_ge_atTop 0] with t ht
    refine Eventually.of_forall fun x => ?_
    rw [Real.norm_of_nonneg (upperDev_nonneg r t x)]
    exact upperDev_anti r x ht
  · refine tendsto_const_nhds.congr' ?_
    filter_upwards [eventually_gt_atTop x] with t ht
    simp [upperDev, not_lt.mpr ht.le]

lemma exists_mUPM_pos [IsProbabilityMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ) :
    ∃ t, 0 < mUPM r t μ := by
  obtain ⟨t, ht⟩ := ((tendsto_cdf_atBot μ).eventually (gt_mem_nhds one_pos)).exists
  refine ⟨t, ?_⟩
  have hI : μ (Ioi t) ≠ 0 := by
    intro h0
    rw [← compl_Iic, prob_compl_eq_zero_iff measurableSet_Iic] at h0
    rw [cdf_eq_real, measureReal_def, h0, ENNReal.toReal_one] at ht
    exact lt_irrefl _ ht
  rw [mUPM, integral_pos_iff_support_of_nonneg (upperDev_nonneg r t) (integrable_upperDev h t)]
  refine (pos_iff_ne_zero.mpr hI).trans_le (measure_mono fun x hx => ?_)
  have hx' : t < x := hx
  simp only [Function.mem_support, upperDev, hx', ↓reduceIte]
  exact pow_ne_zero _ (sub_pos.mpr hx').ne'

lemma exists_mLPM_pos [IsProbabilityMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ) :
    ∃ t, 0 < mLPM r t μ := by
  obtain ⟨t, ht⟩ := ((tendsto_cdf_atTop μ).eventually (lt_mem_nhds one_pos)).exists
  refine ⟨t + 1, ?_⟩
  have hI : μ (Iic t) ≠ 0 := by
    intro h0
    rw [cdf_eq_real, measureReal_def, h0, ENNReal.toReal_zero] at ht
    exact lt_irrefl _ ht
  rw [mLPM, integral_pos_iff_support_of_nonneg (lowerDev_nonneg r _)
    (integrable_lowerDev h _)]
  refine (pos_iff_ne_zero.mpr hI).trans_le (measure_mono fun x hx => ?_)
  have hx' : x ≤ t := hx
  have hx1 : x ≤ t + 1 := by linarith
  simp only [Function.mem_support, lowerDev, hx1, ↓reduceIte]
  exact pow_ne_zero _ (by linarith)

def mPMRatio (r : ℕ) (t : ℝ) (μ : Measure ℝ) : ℝ := mLPM r t μ / (mLPM r t μ + mUPM r t μ)

theorem mPMRatio_nonneg (r : ℕ) (t : ℝ) (μ : Measure ℝ) : 0 ≤ mPMRatio r t μ :=
  div_nonneg (mLPM_nonneg r t μ) (add_nonneg (mLPM_nonneg r t μ) (mUPM_nonneg r t μ))

theorem mPMRatio_le_one (r : ℕ) (t : ℝ) (μ : Measure ℝ) : mPMRatio r t μ ≤ 1 :=
  div_le_one_of_le₀ (le_add_of_nonneg_right (mUPM_nonneg r t μ))
    (add_nonneg (mLPM_nonneg r t μ) (mUPM_nonneg r t μ))

theorem mPMRatio_zero_eq_cdf [IsProbabilityMeasure μ] (t : ℝ) : mPMRatio 0 t μ = cdf μ t := by
  rw [mPMRatio, mLPM_zero_add_mUPM_zero, div_one, mLPM_zero_eq_cdf]

theorem mPMRatio_monotone [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ) :
    Monotone fun t => mPMRatio r t μ := by
  intro s t hst
  have hL : mLPM r s μ ≤ mLPM r t μ := mLPM_monotone h hst
  have hU : mUPM r t μ ≤ mUPM r s μ := mUPM_antitone h hst
  have hLs := mLPM_nonneg r s μ
  have hLt := mLPM_nonneg r t μ
  have hUs := mUPM_nonneg r s μ
  have hUt := mUPM_nonneg r t μ
  change mPMRatio r s μ ≤ mPMRatio r t μ
  rcases (add_nonneg hLs hUs).eq_or_lt with h0 | h0
  · rw [mPMRatio, ← h0, div_zero]; exact mPMRatio_nonneg r t μ
  rcases (add_nonneg hLt hUt).eq_or_lt with h1 | h1
  · have hs0 : mLPM r s μ = 0 := by linarith
    rw [mPMRatio, hs0, zero_div]; exact mPMRatio_nonneg r t μ
  rw [mPMRatio, mPMRatio, div_le_div_iff₀ h0 h1]
  nlinarith [mul_le_mul hL hU hUt hLt]

theorem mLPM_add_mUPM_pos [IsFiniteMeasure μ] {r : ℕ} (h : MemLp (fun x : ℝ => x) r μ)
    {t : ℝ} (ht : μ {t}ᶜ ≠ 0) : 0 < mLPM r t μ + mUPM r t μ := by
  rw [add_comm, mUPM_add_mLPM h, integral_pos_iff_support_of_nonneg
    (fun x => pow_nonneg (abs_nonneg _) r) (integrable_abs_sub_pow h t)]
  refine (pos_iff_ne_zero.mpr ht).trans_le (measure_mono fun x hx => ?_)
  exact pow_ne_zero _ (abs_ne_zero.mpr (sub_ne_zero.mpr hx))

theorem continuous_mPMRatio [IsFiniteMeasure μ] {r : ℕ} (hr : r ≠ 0)
    (h : MemLp (fun x : ℝ => x) r μ) (hμ : ∀ t, μ {t}ᶜ ≠ 0) :
    Continuous fun t => mPMRatio r t μ :=
  (continuous_mLPM hr h).div ((continuous_mLPM hr h).add (continuous_mUPM hr h))
    fun t => (mLPM_add_mUPM_pos h (hμ t)).ne'

theorem tendsto_mPMRatio_atBot [IsProbabilityMeasure μ] {r : ℕ}
    (h : MemLp (fun x : ℝ => x) r μ) : Tendsto (fun t => mPMRatio r t μ) atBot (𝓝 0) := by
  obtain ⟨t₀, ht₀⟩ := exists_mUPM_pos h
  have hlim : Tendsto (fun t => mLPM r t μ / mUPM r t₀ μ) atBot (𝓝 0) := by
    simpa using (tendsto_mLPM_atBot h).div_const (mUPM r t₀ μ)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hlim
    (Eventually.of_forall fun t => mPMRatio_nonneg r t μ) ?_
  filter_upwards [eventually_le_atBot t₀] with t ht
  exact div_le_div_of_nonneg_left (mLPM_nonneg r t μ) ht₀
    ((mUPM_antitone h ht).trans (le_add_of_nonneg_left (mLPM_nonneg r t μ)))

theorem tendsto_mPMRatio_atTop [IsProbabilityMeasure μ] {r : ℕ}
    (h : MemLp (fun x : ℝ => x) r μ) : Tendsto (fun t => mPMRatio r t μ) atTop (𝓝 1) := by
  obtain ⟨t₀, ht₀⟩ := exists_mLPM_pos h
  have hlim : Tendsto (fun t => 1 - mUPM r t μ / mLPM r t₀ μ) atTop (𝓝 1) := by
    simpa using ((tendsto_mUPM_atTop h).div_const (mLPM r t₀ μ)).const_sub 1
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hlim tendsto_const_nhds ?_
    (Eventually.of_forall fun t => mPMRatio_le_one r t μ)
  filter_upwards [eventually_ge_atTop t₀] with t ht
  have hL : mLPM r t₀ μ ≤ mLPM r t μ := mLPM_monotone h ht
  have hU := mUPM_nonneg r t μ
  have hS : 0 < mLPM r t μ + mUPM r t μ := by linarith
  have h1 : mUPM r t μ / (mLPM r t μ + mUPM r t μ) ≤ mUPM r t μ / mLPM r t₀ μ :=
    div_le_div_of_nonneg_left hU ht₀ (by linarith)
  have h2 : mPMRatio r t μ = 1 - mUPM r t μ / (mLPM r t μ + mUPM r t μ) := by
    rw [mPMRatio, eq_sub_iff_add_eq, ← add_div, div_self hS.ne']
  linarith

theorem mPMRatio_one_mean [IsProbabilityMeasure μ] (h : Integrable (fun x : ℝ => x) μ)
    (hμ : μ {∫ y, y ∂μ}ᶜ ≠ 0) : mPMRatio 1 (∫ y, y ∂μ) μ = 1 / 2 := by
  have hpos := mLPM_add_mUPM_pos (memLp_one_of_integrable h) hμ
  have hUL := mUPM_one_sub_mLPM_one h (∫ y, y ∂μ)
  rw [sub_self] at hUL
  rw [mPMRatio, div_eq_iff hpos.ne']
  linarith

lemma integrable_upperDev_of_integrable_pow {r : ℕ} {t : ℝ}
    (h : Integrable (fun x : ℝ => (x - t) ^ r) μ) : Integrable (upperDev r t) μ :=
  h.norm.mono' (upperDev_measurable r t).aestronglyMeasurable (Eventually.of_forall fun x => by
    rw [Real.norm_of_nonneg (upperDev_nonneg r t x), Real.norm_eq_abs, abs_pow]
    exact upperDev_le_abs_pow r t x)

lemma eLPM_add_measure (r : ℕ) (t : ℝ) (μ ν : Measure ℝ) :
    eLPM r t (μ + ν) = eLPM r t μ + eLPM r t ν := lintegral_add_measure _ _ _

lemma eUPM_add_measure (r : ℕ) (t : ℝ) (μ ν : Measure ℝ) :
    eUPM r t (μ + ν) = eUPM r t μ + eUPM r t ν := lintegral_add_measure _ _ _

lemma eLPM_smul_measure (r : ℕ) (t : ℝ) (c : ℝ≥0∞) (μ : Measure ℝ) :
    eLPM r t (c • μ) = c * eLPM r t μ := lintegral_smul_measure _ _

lemma eUPM_smul_measure (r : ℕ) (t : ℝ) (c : ℝ≥0∞) (μ : Measure ℝ) :
    eUPM r t (c • μ) = c * eUPM r t μ := lintegral_smul_measure _ _

def heavyDensity (x : ℝ) : ℝ≥0∞ := ENNReal.ofReal (2 * x ^ (-3 : ℝ))

lemma measurable_heavyDensity : Measurable heavyDensity := by
  unfold heavyDensity; fun_prop

def heavyTail : Measure ℝ := (volume.restrict (Ioi 1)).withDensity heavyDensity

lemma ae_heavyTail : ∀ᵐ x ∂heavyTail, 1 < x :=
  (withDensity_absolutelyContinuous _ _).ae_le (ae_restrict_mem measurableSet_Ioi)

lemma heavyTail_univ : heavyTail univ = 1 := by
  rw [heavyTail, withDensity_apply _ MeasurableSet.univ, Measure.restrict_univ]
  unfold heavyDensity
  rw [← ofReal_integral_eq_lintegral_ofReal
      ((integrableOn_Ioi_rpow_of_lt (by norm_num) one_pos).const_mul 2)]
  · rw [integral_const_mul, integral_Ioi_rpow_of_lt (by norm_num) one_pos]
    norm_num
  · filter_upwards [ae_restrict_mem measurableSet_Ioi] with x hx
    exact mul_nonneg zero_le_two (Real.rpow_nonneg (zero_le_one.trans (mem_Ioi.mp hx).le) _)

lemma pow_mul_heavyDensity (r : ℕ) {x : ℝ} (hx : 1 < x) :
    x ^ r * (heavyDensity x).toReal = 2 * x ^ ((r : ℝ) - 3) := by
  have hx0 : 0 < x := zero_lt_one.trans hx
  unfold heavyDensity
  rw [ENNReal.toReal_ofReal (by positivity), sub_eq_add_neg,
    Real.rpow_add hx0, Real.rpow_natCast]
  ring

lemma integrable_heavyTail_iff {g : ℝ → ℝ} :
    Integrable g heavyTail ↔
      IntegrableOn (fun x => g x * (heavyDensity x).toReal) (Ioi 1) :=
  integrable_withDensity_iff measurable_heavyDensity
    (Eventually.of_forall fun _ => ENNReal.ofReal_lt_top)

lemma integrable_id_heavyTail : Integrable (fun x : ℝ => x) heavyTail := by
  rw [integrable_heavyTail_iff]
  have h0 : IntegrableOn (fun x : ℝ => 2 * x ^ (-2 : ℝ)) (Ioi 1) :=
    (integrableOn_Ioi_rpow_of_lt (by norm_num : (-2 : ℝ) < -1) one_pos).const_mul 2
  refine h0.congr_fun (fun x hx => ?_) measurableSet_Ioi
  have := pow_mul_heavyDensity 1 (mem_Ioi.mp hx)
  rw [pow_one] at this
  rw [this]; norm_num

lemma not_integrable_upperDev_heavyTail {r : ℕ} (hr : 2 ≤ r) :
    ¬ Integrable (upperDev r 0) heavyTail := by
  rw [integrable_heavyTail_iff]
  intro h
  have h' : IntegrableOn (fun x : ℝ => x ^ ((r : ℝ) - 3)) (Ioi 1) := by
    have h0 : IntegrableOn (fun x => 2⁻¹ * (upperDev r 0 x * (heavyDensity x).toReal)) (Ioi 1) :=
      h.const_mul 2⁻¹
    refine h0.congr_fun (fun x hx => ?_) measurableSet_Ioi
    have hx' : 1 < x := mem_Ioi.mp hx
    simp only [upperDev, zero_lt_one.trans hx', ↓reduceIte, sub_zero]
    rw [pow_mul_heavyDensity r hx']
    ring
  rw [integrableOn_Ioi_rpow_iff one_pos] at h'
  have : (2 : ℝ) ≤ r := by exact_mod_cast hr
  linarith

theorem eUPM_heavyTail {r : ℕ} (hr : 2 ≤ r) : eUPM r 0 heavyTail = ∞ := by
  by_contra h
  exact not_integrable_upperDev_heavyTail hr
    ((eUPM_lt_top_iff r 0).mp (lt_top_iff_ne_top.mpr h))

theorem eLPM_map_neg_heavyTail (r : ℕ) :
    eLPM r 0 (heavyTail.map Neg.neg) = eUPM r 0 heavyTail := by
  rw [eLPM, lintegral_map (lowerDev_measurable r 0).ennreal_ofReal measurable_neg, eUPM]
  refine lintegral_congr_ae ?_
  filter_upwards [ae_heavyTail] with x hx
  have hx0 : 0 < x := zero_lt_one.trans hx
  simp [lowerDev, upperDev, hx0, hx0.le]

def heavyTailSym : Measure ℝ :=
  (2⁻¹ : ℝ≥0∞) • heavyTail + (2⁻¹ : ℝ≥0∞) • heavyTail.map Neg.neg

instance : IsProbabilityMeasure heavyTailSym := by
  constructor
  simp only [heavyTailSym, Measure.add_apply, Measure.smul_apply,
    Measure.map_apply measurable_neg MeasurableSet.univ, preimage_univ, heavyTail_univ,
    smul_eq_mul, mul_one]
  exact ENNReal.inv_two_add_inv_two

lemma integrable_id_heavyTailSym : Integrable (fun x : ℝ => x) heavyTailSym := by
  have h2 : Integrable (fun x : ℝ => x) (heavyTail.map Neg.neg) := by
    rw [integrable_map_measure (g := fun x : ℝ => x) measurable_id.aestronglyMeasurable
      measurable_neg.aemeasurable]
    exact integrable_id_heavyTail.neg
  exact (integrable_id_heavyTail.smul_measure (by simp)).add_measure
    (h2.smul_measure (by simp))

theorem eUPM_heavyTailSym {r : ℕ} (hr : 2 ≤ r) : eUPM r 0 heavyTailSym = ∞ := by
  rw [heavyTailSym, eUPM_add_measure, eUPM_smul_measure, eUPM_heavyTail hr]
  simp [ENNReal.mul_top]

theorem eLPM_heavyTailSym {r : ℕ} (hr : 2 ≤ r) : eLPM r 0 heavyTailSym = ∞ := by
  rw [heavyTailSym, eLPM_add_measure, eLPM_smul_measure, eLPM_smul_measure,
    eLPM_map_neg_heavyTail, eUPM_heavyTail hr]
  simp [ENNReal.mul_top]

theorem exists_integrable_eUPM_eLPM_eq_top :
    ∃ μ : Measure ℝ, IsProbabilityMeasure μ ∧ Integrable (fun x : ℝ => x) μ ∧
      ∀ r : ℕ, 2 ≤ r → eUPM r 0 μ = ∞ ∧ eLPM r 0 μ = ∞ ∧
        ¬ Integrable (fun x : ℝ => x ^ r) μ ∧ ¬ MemLp (fun x : ℝ => x) r μ := by
  refine ⟨heavyTailSym, inferInstance, integrable_id_heavyTailSym, fun r hr => ?_⟩
  have hU := eUPM_heavyTailSym hr
  refine ⟨hU, eLPM_heavyTailSym hr, fun h => ?_, fun h => ?_⟩
  · have h' : Integrable (fun x : ℝ => (x - 0) ^ r) heavyTailSym := by simpa using h
    exact (eUPM_lt_top_iff r 0).mpr (integrable_upperDev_of_integrable_pow h') |>.ne hU
  · exact (eUPM_lt_top h 0).ne hU

end

end NNS
