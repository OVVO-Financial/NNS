-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.RealDegree
public import NNSProofs.Moments
public import Mathlib.Probability.Moments.Variance
public import Mathlib.Probability.CDF

namespace NNS

open MeasureTheory ProbabilityTheory Set Finset

@[expose] public noncomputable section

def lpm (d : ℕ) (t : ℝ) (μ : Measure ℝ) : ℝ := ∫ v, lowerDev d t v ∂μ

def upm (d : ℕ) (t : ℝ) (μ : Measure ℝ) : ℝ := ∫ v, upperDev d t v ∂μ

def rlpm (α t : ℝ) (μ : Measure ℝ) : ℝ := ∫ v, rLowerDev α t v ∂μ

def theta (t α : ℝ) (μ : Measure ℝ) : ℝ := rlpm α t μ ^ (1 / α)

def pmean (μ : Measure ℝ) : ℝ := ∫ v, v ∂μ

def pvar (μ : Measure ℝ) : ℝ := ∫ v, (v - pmean μ) ^ 2 ∂μ

def semivar (μ : Measure ℝ) : ℝ := lpm 2 (pmean μ) μ

def empirical {n : ℕ} (x : Fin n → ℝ) : Measure ℝ :=
  (n : ENNReal)⁻¹ • ∑ i, Measure.dirac (x i)

variable {n : ℕ}

theorem measurable_lowerDev (d : ℕ) (t : ℝ) : Measurable (lowerDev d t) := by
  change Measurable fun v => if v ≤ t then (t - v) ^ d else 0
  exact Measurable.ite measurableSet_Iic (by fun_prop) measurable_const

theorem measurable_upperDev (d : ℕ) (t : ℝ) : Measurable (upperDev d t) := by
  change Measurable fun v => if t < v then (v - t) ^ d else 0
  exact Measurable.ite measurableSet_Ioi (by fun_prop) measurable_const

theorem measurable_rLowerDev (α t : ℝ) : Measurable (rLowerDev α t) := by
  change Measurable fun v => if v ≤ t then (t - v) ^ α else 0
  exact Measurable.ite measurableSet_Iic (by fun_prop) measurable_const

theorem lpm_nonneg (d : ℕ) (t : ℝ) (μ : Measure ℝ) : 0 ≤ lpm d t μ :=
  integral_nonneg (lowerDev_nonneg d t)

theorem upm_nonneg (d : ℕ) (t : ℝ) (μ : Measure ℝ) : 0 ≤ upm d t μ :=
  integral_nonneg (upperDev_nonneg d t)

theorem rlpm_nonneg (α t : ℝ) (μ : Measure ℝ) : 0 ≤ rlpm α t μ :=
  integral_nonneg (rLowerDev_nonneg α t)

theorem pvar_nonneg (μ : Measure ℝ) : 0 ≤ pvar μ :=
  integral_nonneg fun _ => sq_nonneg _

theorem semivar_nonneg (μ : Measure ℝ) : 0 ≤ semivar μ := lpm_nonneg _ _ _

theorem rlpm_natCast (d : ℕ) (t : ℝ) (μ : Measure ℝ) : rlpm d t μ = lpm d t μ := by
  unfold rlpm lpm
  simp_rw [rLowerDev_natCast]

theorem lpm_zero (t : ℝ) (μ : Measure ℝ) : lpm 0 t μ = μ.real (Iic t) := by
  rw [← integral_indicator_one measurableSet_Iic]
  unfold lpm
  congr 1

theorem lpm_zero_eq_cdf (t : ℝ) (μ : Measure ℝ) [IsProbabilityMeasure μ] :
    lpm 0 t μ = cdf μ t := by
  rw [lpm_zero, cdf_eq_real]

theorem pvar_eq_variance (μ : Measure ℝ) : pvar μ = variance id μ := by
  rw [variance_eq_integral measurable_id.aemeasurable]
  rfl

theorem integrable_empirical (x : Fin n → ℝ) (f : ℝ → ℝ) : Integrable f (empirical x) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [empirical]
  refine Integrable.smul_measure ?_ (by simp [hn.ne'])
  exact integrable_finsetSum_measure.2 fun i _ => integrable_dirac (by simp)

theorem integral_empirical (x : Fin n → ℝ) (f : ℝ → ℝ) :
    ∫ v, f v ∂(empirical x) = (∑ i, f (x i)) / n := by
  unfold empirical
  rw [integral_smul_measure, integral_finsetSum_measure fun i _ => integrable_dirac (by simp)]
  simp [integral_dirac, div_eq_inv_mul]

theorem isProbabilityMeasure_empirical (hn : 0 < n) (x : Fin n → ℝ) :
    IsProbabilityMeasure (empirical x) := by
  constructor
  simp only [empirical, Measure.smul_apply, Measure.coe_finsetSum, Finset.sum_apply,
    measure_univ, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul,
    mul_one, smul_eq_mul]
  exact ENNReal.inv_mul_cancel (by simp [hn.ne']) (by simp)

theorem memLp_empirical (x : Fin n → ℝ) : MemLp id 2 (empirical x) :=
  (memLp_two_iff_integrable_sq measurable_id.aestronglyMeasurable).2
    (integrable_empirical x _)

theorem lpm_empirical (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : lpm d t (empirical x) = LPM d t x := by
  unfold lpm LPM
  rw [integral_empirical]

theorem pmean_empirical (x : Fin n → ℝ) : pmean (empirical x) = mean x := by
  unfold pmean mean
  rw [integral_empirical]

theorem pvar_empirical (x : Fin n → ℝ) : pvar (empirical x) = popVariance x := by
  unfold pvar popVariance centralMoment
  rw [integral_empirical, pmean_empirical]

theorem semivar_empirical (x : Fin n → ℝ) : semivar (empirical x) = LPM 2 (mean x) x := by
  unfold semivar
  rw [lpm_empirical, pmean_empirical]

theorem real_Iic_empirical (t : ℝ) (x : Fin n → ℝ) : (empirical x).real (Iic t) = LPM 0 t x := by
  rw [← lpm_zero, lpm_empirical]

theorem real_le_integral_div_of_le {μ : Measure ℝ} [IsFiniteMeasure μ] {f : ℝ → ℝ}
    (hf0 : ∀ v, 0 ≤ f v) (hfi : Integrable f μ) {c : ℝ} (hc : 0 < c) {s : Set ℝ}
    (hs : ∀ v ∈ s, c ≤ f v) : μ.real s ≤ (∫ v, f v ∂μ) / c := by
  rw [le_div_iff₀ hc, mul_comm]
  calc c * μ.real s ≤ c * μ.real {v | c ≤ f v} :=
        mul_le_mul_of_nonneg_left (measureReal_mono fun v hv => hs v hv) hc.le
    _ ≤ ∫ v, f v ∂μ := mul_meas_ge_le_integral_of_nonneg (.of_forall hf0) hfi c

variable {μ : Measure ℝ} {α t g : ℝ}

theorem real_Iic_le_rlpm_div [IsFiniteMeasure μ] (hα : 0 ≤ α) (hgt : g < t)
    (hint : Integrable (rLowerDev α t) μ) :
    μ.real (Iic g) ≤ rlpm α t μ / (t - g) ^ α := by
  refine real_le_integral_div_of_le (rLowerDev_nonneg α t) hint
    (Real.rpow_pos_of_pos (sub_pos.2 hgt) α) fun v hv => ?_
  have hvg : v ≤ g := hv
  simp only [rLowerDev, hvg.trans hgt.le, ↓reduceIte]
  exact Real.rpow_le_rpow (sub_pos.2 hgt).le (by linarith) hα

theorem real_Iic_le_theta [IsFiniteMeasure μ] (hα : 0 < α) (hgt : g < t)
    (hint : Integrable (rLowerDev α t) μ) :
    μ.real (Iic g) ≤ (theta t α μ / (t - g)) ^ α := by
  unfold theta
  rw [Real.div_rpow (Real.rpow_nonneg (rlpm_nonneg α t μ) _) (sub_pos.2 hgt).le,
    ← Real.rpow_mul (rlpm_nonneg α t μ), one_div_mul_cancel hα.ne', Real.rpow_one]
  exact real_Iic_le_rlpm_div hα.le hgt hint

theorem theta_mean_two_div (μ : Measure ℝ) (g : ℝ) :
    (theta (pmean μ) 2 μ / (pmean μ - g)) ^ (2 : ℝ) = semivar μ / (pmean μ - g) ^ 2 := by
  have h2 : ((2 : ℕ) : ℝ) = 2 := by norm_num
  have hl : rlpm 2 (pmean μ) μ = semivar μ := by
    rw [← h2, rlpm_natCast]; rfl
  unfold theta
  rw [Real.rpow_two, div_pow, ← Real.rpow_two, ← Real.rpow_mul (rlpm_nonneg _ _ _)]
  norm_num [hl]

theorem real_Iic_le_lpm_div [IsFiniteMeasure μ] (d : ℕ) (hgt : g < t)
    (hint : Integrable (lowerDev d t) μ) :
    μ.real (Iic g) ≤ lpm d t μ / (t - g) ^ d := by
  refine real_le_integral_div_of_le (lowerDev_nonneg d t) hint
    (pow_pos (sub_pos.2 hgt) d) fun v hv => ?_
  have hvg : v ≤ g := hv
  simp only [lowerDev, hvg.trans hgt.le, ↓reduceIte]
  exact pow_le_pow_left₀ (sub_pos.2 hgt).le (by linarith) d

theorem measure_Iic_le_lintegral_div (μ : Measure ℝ) (hα : 0 ≤ α) (hgt : g < t) :
    μ (Iic g) ≤
      (∫⁻ v, ENNReal.ofReal (rLowerDev α t v) ∂μ) / ENNReal.ofReal ((t - g) ^ α) := by
  have hc : 0 < (t - g) ^ α := Real.rpow_pos_of_pos (sub_pos.2 hgt) α
  refine le_trans (measure_mono fun v hv => ?_) (meas_ge_le_lintegral_div
    (measurable_rLowerDev α t).ennreal_ofReal.aemeasurable
    (by simpa using hc) ENNReal.ofReal_ne_top)
  have hvg : v ≤ g := hv
  change ENNReal.ofReal ((t - g) ^ α) ≤ ENNReal.ofReal (rLowerDev α t v)
  refine ENNReal.ofReal_le_ofReal ?_
  simp only [rLowerDev, hvg.trans hgt.le, ↓reduceIte]
  exact Real.rpow_le_rpow (sub_pos.2 hgt).le (by linarith) hα

theorem rlpm_dirac_self (hα : 0 < α) (t : ℝ) : rlpm α t (Measure.dirac t) = 0 := by
  simp [rlpm, rLowerDev, Real.zero_rpow hα.ne']

section L2

variable (hμ : MemLp id 2 μ)
include hμ

theorem integrable_id_of_memLp [IsFiniteMeasure μ] : Integrable (fun v : ℝ => v) μ :=
  hμ.integrable (by norm_num)

theorem integrable_sq_sub_of_memLp [IsFiniteMeasure μ] (m : ℝ) :
    Integrable (fun v : ℝ => (v - m) ^ 2) μ :=
  (hμ.sub (memLp_const m)).integrable_sq

theorem integrable_lowerDev_two [IsFiniteMeasure μ] (m : ℝ) :
    Integrable (lowerDev 2 m) μ := by
  refine (integrable_sq_sub_of_memLp hμ m).mono'
    (measurable_lowerDev 2 m).aestronglyMeasurable (.of_forall fun v => ?_)
  rw [Real.norm_of_nonneg (lowerDev_nonneg 2 m v)]
  unfold lowerDev
  split_ifs <;> nlinarith [sq_nonneg (v - m)]

theorem integrable_upperDev_two [IsFiniteMeasure μ] (m : ℝ) :
    Integrable (upperDev 2 m) μ := by
  refine (integrable_sq_sub_of_memLp hμ m).mono'
    (measurable_upperDev 2 m).aestronglyMeasurable (.of_forall fun v => ?_)
  rw [Real.norm_of_nonneg (upperDev_nonneg 2 m v)]
  unfold upperDev
  split_ifs <;> nlinarith [sq_nonneg (v - m)]

theorem lpm_two_add_upm_two [IsFiniteMeasure μ] (m : ℝ) :
    lpm 2 m μ + upm 2 m μ = ∫ v, (v - m) ^ 2 ∂μ := by
  unfold lpm upm
  rw [← integral_add (integrable_lowerDev_two hμ m) (integrable_upperDev_two hμ m)]
  congr 1
  ext v
  unfold lowerDev upperDev
  split_ifs with h1 h2 h2
  · linarith
  · ring
  · ring
  · linarith

theorem semivar_add_upm_two [IsFiniteMeasure μ] :
    semivar μ + upm 2 (pmean μ) μ = pvar μ :=
  lpm_two_add_upm_two hμ _

theorem semivar_le_pvar [IsFiniteMeasure μ] : semivar μ ≤ pvar μ := by
  rw [← semivar_add_upm_two hμ]
  linarith [upm_nonneg 2 (pmean μ) μ]

theorem real_Iic_le_semivar_div [IsFiniteMeasure μ] (hg : g < pmean μ) :
    μ.real (Iic g) ≤ semivar μ / (pmean μ - g) ^ 2 :=
  real_Iic_le_lpm_div 2 hg (integrable_lowerDev_two hμ _)

theorem semivar_div_le_pvar_div [IsFiniteMeasure μ] :
    semivar μ / (pmean μ - g) ^ 2 ≤ pvar μ / (pmean μ - g) ^ 2 :=
  div_le_div_of_nonneg_right (semivar_le_pvar hμ) (sq_nonneg _)

theorem cantelli [IsProbabilityMeasure μ] (hg : g < pmean μ) :
    μ.real (Iic g) ≤ pvar μ / (pvar μ + (pmean μ - g) ^ 2) := by
  set m := pmean μ with hm
  set V := pvar μ with hV
  have hl : 0 < m - g := sub_pos.2 hg
  have hV0 : 0 ≤ V := pvar_nonneg μ
  set u := V / (m - g) with hu
  have hu0 : 0 ≤ u := div_nonneg hV0 hl.le
  have hid := integrable_id_of_memLp hμ
  have hsq := integrable_sq_sub_of_memLp hμ m
  have hexp : (fun v : ℝ => (m - v + u) ^ 2) =
      fun v => (v - m) ^ 2 + (2 * u * m + u ^ 2) - 2 * u * v := by
    ext v; ring
  have h1 : Integrable (fun v : ℝ => (v - m) ^ 2 + (2 * u * m + u ^ 2)) μ :=
    hsq.add (integrable_const _)
  have h2 : Integrable (fun v : ℝ => 2 * u * v) μ := hid.const_mul _
  have hf : Integrable (fun v : ℝ => (m - v + u) ^ 2) μ := by
    rw [hexp]
    exact h1.sub h2
  have hint : ∫ v, (m - v + u) ^ 2 ∂μ = V + u ^ 2 := by
    rw [hexp, integral_sub h1 h2,
      integral_add hsq (integrable_const _), integral_const, integral_const_mul]
    simp only [probReal_univ, one_smul]
    have e1 : ∫ v, (v - m) ^ 2 ∂μ = V := rfl
    have e2 : ∫ v, v ∂μ = m := rfl
    rw [e1, e2, hu]
    ring
  have key := real_le_integral_div_of_le (fun v => sq_nonneg (m - v + u)) hf
    (c := (m - g + u) ^ 2) (by positivity) (s := Iic g) fun v hv => by
      have hvg : v ≤ g := hv
      exact pow_le_pow_left₀ (by linarith) (by linarith) 2
  rw [hint] at key
  refine key.trans (le_of_eq ?_)
  rw [hu, div_eq_div_iff (by positivity) (by positivity)]
  field_simp
  ring

end L2

theorem upm_two_eq_semivar_of_symm (hsym : μ.map (fun v => 2 * pmean μ - v) = μ) :
    upm 2 (pmean μ) μ = semivar μ := by
  unfold semivar upm lpm
  generalize pmean μ = m at hsym ⊢
  have hφ : Measurable fun v : ℝ => 2 * m - v := by fun_prop
  calc ∫ v, upperDev 2 m v ∂μ = ∫ v, upperDev 2 m v ∂(μ.map fun v => 2 * m - v) := by
        rw [hsym]
    _ = ∫ v, upperDev 2 m (2 * m - v) ∂μ :=
        integral_map hφ.aemeasurable (measurable_upperDev 2 m).aestronglyMeasurable
    _ = ∫ v, lowerDev 2 m v ∂μ := by
        congr 1
        ext v
        unfold upperDev lowerDev
        split_ifs with h1 h2 h2
        · ring
        · linarith
        · have : v = m := by linarith
          subst this; ring
        · rfl

section L2

variable (hμ : MemLp id 2 μ)
include hμ

theorem semivar_eq_half_pvar_of_symm [IsFiniteMeasure μ]
    (hsym : μ.map (fun v => 2 * pmean μ - v) = μ) : semivar μ = pvar μ / 2 := by
  have h := semivar_add_upm_two hμ
  rw [upm_two_eq_semivar_of_symm hsym] at h
  linarith

theorem halved_chebyshev_of_semivar_le [IsFiniteMeasure μ] (hg : g < pmean μ)
    (hs : semivar μ ≤ pvar μ / 2) :
    μ.real (Iic g) ≤ pvar μ / (2 * (pmean μ - g) ^ 2) := by
  refine (real_Iic_le_semivar_div hμ hg).trans ?_
  rw [div_le_div_iff₀ (pow_pos (sub_pos.2 hg) 2) (by have := sub_pos.2 hg; positivity)]
  nlinarith [sq_nonneg (pmean μ - g)]

theorem halved_chebyshev_of_symm [IsFiniteMeasure μ] (hg : g < pmean μ)
    (hsym : μ.map (fun v => 2 * pmean μ - v) = μ) :
    μ.real (Iic g) ≤ pvar μ / (2 * (pmean μ - g) ^ 2) :=
  halved_chebyshev_of_semivar_le hμ hg (semivar_eq_half_pvar_of_symm hμ hsym).le

end L2

section RandomVariable

variable {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {X : Ω → ℝ}

theorem pmean_map (hX : AEMeasurable X P) : pmean (P.map X) = ∫ ω, X ω ∂P :=
  integral_map hX measurable_id.aestronglyMeasurable

theorem pvar_map (hX : AEMeasurable X P) : pvar (P.map X) = variance X P := by
  rw [variance_eq_integral hX]
  unfold pvar
  rw [pmean_map hX]
  exact integral_map hX (by fun_prop)

theorem real_map_Iic (hX : AEMeasurable X P) (g : ℝ) :
    (P.map X).real (Iic g) = P.real {ω | X ω ≤ g} := by
  rw [measureReal_def, Measure.map_apply_of_aemeasurable hX measurableSet_Iic]
  rfl

theorem rv_cantelli [IsProbabilityMeasure P] (hX : MemLp X 2 P) (hg : g < ∫ ω, X ω ∂P) :
    P.real {ω | X ω ≤ g} ≤ variance X P / (variance X P + ((∫ ω, X ω ∂P) - g) ^ 2) := by
  have hXm : AEMeasurable X P := hX.aemeasurable
  have hL : MemLp id 2 (P.map X) :=
    (memLp_map_measure_iff measurable_id.aestronglyMeasurable hXm).2 hX
  have := (Measure.isProbabilityMeasure_map_iff hXm).2 inferInstance
  rw [← real_map_Iic hXm, ← pvar_map hXm, ← pmean_map hXm]
  exact cantelli hL (by rwa [pmean_map hXm])

end RandomVariable

def zSample : Fin 4 → ℝ := ![-1, 1 / 3, 1 / 3, 1 / 3]

theorem zSample_facts :
    mean zSample = 0 ∧ popVariance zSample = 1 / 3 ∧ LPM 0 (-1) zSample = 1 / 4 ∧
      LPM 2 (mean zSample) zSample = 1 / 4 := by
  have hm : mean zSample = 0 := by
    simp [mean, zSample, Fin.sum_univ_four]; norm_num
  refine ⟨hm, ?_, ?_, ?_⟩
  · rw [popVariance, hm, centralMoment, Fin.sum_univ_four]
    simp [zSample]; norm_num
  · rw [LPM, Fin.sum_univ_four]
    simp [lowerDev, zSample]; norm_num
  · rw [hm, LPM, Fin.sum_univ_four]
    simp [lowerDev, zSample]; norm_num

def wSample : Fin 100 → ℝ := fun i => if i = 0 then -99 else 1

theorem wSample_facts :
    mean wSample = 0 ∧ popVariance wSample = 99 ∧ LPM 0 (-50) wSample = 1 / 100 ∧
      LPM 2 (mean wSample) wSample = 9801 / 100 := by
  have hm : mean wSample = 0 := by
    simp only [mean]
    rw [Fin.sum_univ_succ]
    simp [wSample, Fin.succ_ne_zero]
  refine ⟨hm, ?_, ?_, ?_⟩
  · simp only [popVariance, centralMoment, hm]
    rw [Fin.sum_univ_succ]
    simp [wSample, Fin.succ_ne_zero]
    norm_num
  · simp only [LPM]
    rw [Fin.sum_univ_succ]
    simp [wSample, Fin.succ_ne_zero, lowerDev]
    norm_num
  · simp only [LPM, hm]
    rw [Fin.sum_univ_succ]
    simp [wSample, Fin.succ_ne_zero, lowerDev]
    norm_num

theorem zSample_bounds :
    (empirical zSample).real (Iic (-1)) = 1 / 4 ∧
      pvar (empirical zSample) / (2 * (pmean (empirical zSample) - (-1)) ^ 2) = 1 / 6 ∧
      pvar (empirical zSample) /
        (pvar (empirical zSample) + (pmean (empirical zSample) - (-1)) ^ 2) = 1 / 4 ∧
      semivar (empirical zSample) / (pmean (empirical zSample) - (-1)) ^ 2 = 1 / 4 := by
  obtain ⟨hm, hv, hp, hs⟩ := zSample_facts
  rw [hm] at hs
  simp only [real_Iic_empirical, pvar_empirical, pmean_empirical, semivar_empirical, hm, hv,
    hp, hs]
  norm_num

theorem wSample_bounds :
    (empirical wSample).real (Iic (-50)) = 1 / 100 ∧
      pvar (empirical wSample) /
        (pvar (empirical wSample) + (pmean (empirical wSample) - (-50)) ^ 2) = 99 / 2599 ∧
      semivar (empirical wSample) / (pmean (empirical wSample) - (-50)) ^ 2 =
        9801 / 250000 ∧
      (99 : ℝ) / 2599 < 9801 / 250000 := by
  obtain ⟨hm, hv, hp, hs⟩ := wSample_facts
  rw [hm] at hs
  simp only [real_Iic_empirical, pvar_empirical, pmean_empirical, semivar_empirical, hm, hv,
    hp, hs]
  norm_num

def vSample : Fin 4 → ℝ := ![1, 1, 1, 13]

theorem vSample_facts :
    mean vSample = 4 ∧ popVariance vSample = 27 ∧ LPM 2 (mean vSample) vSample = 27 / 4 := by
  have hm : mean vSample = 4 := by
    norm_num [mean, vSample, Fin.sum_univ_four]
  refine ⟨hm, ?_, ?_⟩
  · rw [popVariance, hm, centralMoment, Fin.sum_univ_four]
    norm_num [vSample]
  · rw [hm, LPM, Fin.sum_univ_four]
    norm_num [lowerDev, vSample]

theorem vSample_bounds :
    pvar (empirical vSample) /
        (pvar (empirical vSample) + (pmean (empirical vSample) - (-1)) ^ 2) = 27 / 52 ∧
      semivar (empirical vSample) / (pmean (empirical vSample) - (-1)) ^ 2 = 27 / 100 := by
  obtain ⟨hm, hv, hs⟩ := vSample_facts
  rw [hm] at hs
  simp only [pvar_empirical, pmean_empirical, semivar_empirical, hm, hv, hs]
  norm_num

end

end NNS
