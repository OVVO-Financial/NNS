-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.MeasurePM
public import NNSProofs.Bounds
public import Mathlib.Probability.Moments.Basic

namespace NNS.GapOperators

open MeasureTheory ProbabilityTheory Finset

@[expose] public noncomputable section

variable {n : ℕ} {μ : Measure ℝ}

theorem posPart_eq_max (x : ℝ) : x⁺ = max x 0 := posPart_def x

theorem pmax_zero_apply (x : Fin n → ℝ) (i : Fin n) : x⁺ i = max (x i) 0 := by
  rw [Pi.posPart_apply, posPart_def]

theorem upperDev_one_eq_posPart (t x : ℝ) :
    upperDev 1 t x = (x - t)⁺ ∧ (x - t)⁺ = max (x - t) 0 :=
  ⟨by rw [posPart_def, upperDev_one], posPart_def _⟩

theorem lowerDev_one_eq_posPart (t x : ℝ) :
    lowerDev 1 t x = (t - x)⁺ ∧ (t - x)⁺ = max (t - x) 0 :=
  ⟨by rw [posPart_def, lowerDev_one], posPart_def _⟩

theorem upperDev_zero (t x : ℝ) : upperDev 0 t x = if t < x then 1 else 0 := by
  simp [upperDev]

theorem literal_degree_zero_sample (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    (∑ i, max (t - x i) 0 ^ 0) / n = 1 ∧ (∑ i, max (x i - t) 0 ^ 0) / n = 1 := by
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simp [hn']

theorem literal_degree_zero_integral [IsProbabilityMeasure μ] (t : ℝ) :
    ∫ x, max (t - x) 0 ^ 0 ∂μ = 1 ∧ ∫ x, max (x - t) 0 ^ 0 ∂μ = 1 := by
  simp

theorem literal_degree_zero_sum_sample (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    (∑ i, max (t - x i) 0 ^ 0) / n + (∑ i, max (x i - t) 0 ^ 0) / n = 2 ∧
      LPM 0 t x + UPM 0 t x = 1 := by
  obtain ⟨h1, h2⟩ := literal_degree_zero_sample hn t x
  exact ⟨by rw [h1, h2]; norm_num, LPM_zero_add_UPM_zero hn t x⟩

theorem literal_degree_zero_sum_integral [IsProbabilityMeasure μ] (t : ℝ) :
    ∫ x, max (t - x) 0 ^ 0 ∂μ + ∫ x, max (x - t) 0 ^ 0 ∂μ = 2 ∧
      mLPM 0 t μ + mUPM 0 t μ = 1 := by
  obtain ⟨h1, h2⟩ := literal_degree_zero_integral (μ := μ) t
  exact ⟨by rw [h1, h2]; norm_num, mLPM_zero_add_mUPM_zero t⟩

theorem LPM_one_eq_sum_max (t : ℝ) (x : Fin n → ℝ) :
    LPM 1 t x = 1 / n * ∑ i, max (t - x i) 0 := by
  rw [one_div, ← div_eq_inv_mul]
  simp only [LPM, lowerDev_one]

theorem UPM_one_eq_sum_max (t : ℝ) (x : Fin n → ℝ) :
    UPM 1 t x = 1 / n * ∑ i, max (x i - t) 0 := by
  rw [one_div, ← div_eq_inv_mul]
  simp only [UPM, upperDev_one]

theorem mLPM_one_eq_integral_max (t : ℝ) (μ : Measure ℝ) :
    mLPM 1 t μ = ∫ x, max (t - x) 0 ∂μ := by
  simp only [mLPM, lowerDev_one]

theorem mUPM_one_eq_integral_max (t : ℝ) (μ : Measure ℝ) :
    mUPM 1 t μ = ∫ x, max (x - t) 0 ∂μ := by
  simp only [mUPM, upperDev_one]

theorem rlpm_eq_integral_max_rpow {r : ℝ} (hr : 0 < r) (τ : ℝ) (μ : Measure ℝ) :
    rlpm r τ μ = ∫ x, max (τ - x) 0 ^ r ∂μ := by
  simp only [rlpm, rLowerDev_eq_max hr.ne']

theorem rLPM_eq_sum_max_rpow {r : ℝ} (hr : 0 < r) (τ : ℝ) (x : Fin n → ℝ) :
    rLPM r τ x = (∑ i, max (τ - x i) 0 ^ r) / n :=
  rLPM_eq_max hr.ne' τ x

theorem rlpm_zero (τ : ℝ) (μ : Measure ℝ) : rlpm 0 τ μ = μ.real (Set.Iic τ) := by
  have h := rlpm_natCast 0 τ μ
  rw [Nat.cast_zero] at h
  rw [h, lpm_zero]

theorem rLPM_zero_eq_ecdf (τ : ℝ) (x : Fin n → ℝ) :
    rLPM 0 τ x = (#{i | x i ≤ τ} : ℝ) / n := by
  simp only [rLPM, rLowerDev_zero]
  rw [sum_boole]

theorem rlpm_zero_ne_literal :
    rlpm 0 0 (Measure.dirac (1 : ℝ)) = 0 ∧
      ∫ x, max (0 - x) 0 ^ (0 : ℝ) ∂(Measure.dirac (1 : ℝ)) = 1 := by
  constructor
  · simp [rlpm, rLowerDev]
  · simp

def rupm (d t : ℝ) (μ : Measure ℝ) : ℝ := ∫ v, rUpperDev d t v ∂μ

theorem rlpm_eq_mLPM (d : ℕ) (t : ℝ) (μ : Measure ℝ) : rlpm d t μ = mLPM d t μ := by
  simp only [rlpm, mLPM, rLowerDev_natCast]

theorem rupm_eq_mUPM (d : ℕ) (t : ℝ) (μ : Measure ℝ) : rupm d t μ = mUPM d t μ := by
  simp only [rupm, mUPM, rUpperDev_natCast]

theorem rupm_eq_integral_max_rpow {r : ℝ} (hr : r ≠ 0) (t : ℝ) (μ : Measure ℝ) :
    rupm r t μ = ∫ x, max (x - t) 0 ^ r ∂μ := by
  simp only [rupm, rUpperDev_eq_max hr]

theorem rupm_zero (t : ℝ) (μ : Measure ℝ) : rupm 0 t μ = μ.real (Set.Ioi t) := by
  have h : rUpperDev 0 t = (Set.Ioi t).indicator 1 := by
    ext v; simp [rUpperDev_zero, Set.indicator]
  rw [rupm, h, integral_indicator_one measurableSet_Ioi]

theorem rUPM_zero_eq_card (t : ℝ) (x : Fin n → ℝ) :
    rUPM 0 t x = (#{i | t < x i} : ℝ) / n := by
  simp only [rUPM, rUpperDev_zero]
  rw [sum_boole]

theorem centralMoment_eq_sum (r : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    centralMoment r t x = (∑ i, (x i - t) ^ r) / n := rfl

theorem integral_sub_pow_mean_eq_centralMoment (r : ℕ) (μ : Measure ℝ) :
    ∫ x, (x - ∫ y, y ∂μ) ^ r ∂μ = ProbabilityTheory.centralMoment id r μ := rfl

theorem integral_sub_pow_empMeasure (r : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    ∫ v, (v - t) ^ r ∂empMeasure x = centralMoment r t x :=
  integral_empMeasure x _

theorem variance_eq_integral_sq {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {X : Ω → ℝ}
    (hX : AEMeasurable X P) : variance X P = ∫ ω, (X ω - ∫ ω', X ω' ∂P) ^ 2 ∂P :=
  variance_eq_integral hX

theorem mLPM_zero_eq_real (t : ℝ) (μ : Measure ℝ) : mLPM 0 t μ = μ.real (Set.Iic t) := by
  rw [← integral_indicator_one measurableSet_Iic, mLPM, lowerDev_eq_indicator]
  simp only [pow_zero]
  rfl

theorem setIntegral_eq_real_mul_cond [IsFiniteMeasure μ] (s : Set ℝ) (f : ℝ → ℝ) :
    ∫ x in s, f x ∂μ = μ.real s * ∫ x, f x ∂μ[|s] := by
  rw [ProbabilityTheory.cond, integral_smul_measure, smul_eq_mul, ← mul_assoc]
  by_cases h : μ s = 0
  · simp [Measure.restrict_eq_zero.mpr h]
  · rw [Measure.real, ENNReal.toReal_inv,
      mul_inv_cancel₀ (ENNReal.toReal_ne_zero.mpr ⟨h, measure_ne_top μ s⟩), one_mul]

theorem mLPM_one_eq_real_mul_cond [IsFiniteMeasure μ] (t : ℝ) :
    mLPM 1 t μ = μ.real (Set.Iic t) * ∫ x, (t - x) ∂μ[|Set.Iic t] := by
  rw [mLPM_eq_setIntegral, ← setIntegral_eq_real_mul_cond]
  simp only [pow_one]

theorem mUPM_one_eq_real_mul_cond [IsFiniteMeasure μ] (t : ℝ) :
    mUPM 1 t μ = μ.real (Set.Ioi t) * ∫ x, (x - t) ∂μ[|Set.Ioi t] := by
  rw [mUPM_eq_setIntegral, ← setIntegral_eq_real_mul_cond]
  simp only [pow_one]

theorem mLPM_one_le_cond [IsProbabilityMeasure μ] (t : ℝ) :
    mLPM 1 t μ ≤ ∫ x, (t - x) ∂μ[|Set.Iic t] := by
  rw [mLPM_one_eq_real_mul_cond]
  refine mul_le_of_le_one_left ?_ measureReal_le_one
  rw [ProbabilityTheory.cond, integral_smul_measure]
  exact smul_nonneg ENNReal.toReal_nonneg
    (setIntegral_nonneg measurableSet_Iic fun x hx => sub_nonneg.mpr hx)

theorem mUPM_one_le_cond [IsProbabilityMeasure μ] (t : ℝ) :
    mUPM 1 t μ ≤ ∫ x, (x - t) ∂μ[|Set.Ioi t] := by
  rw [mUPM_one_eq_real_mul_cond]
  refine mul_le_of_le_one_left ?_ measureReal_le_one
  rw [ProbabilityTheory.cond, integral_smul_measure]
  exact smul_nonneg ENNReal.toReal_nonneg
    (setIntegral_nonneg measurableSet_Ioi fun x hx => sub_nonneg.mpr (le_of_lt hx))

theorem mLPM_one_div_mLPM_zero [IsFiniteMeasure μ] (t : ℝ) (h : μ (Set.Iic t) ≠ 0) :
    mLPM 1 t μ / mLPM 0 t μ = ∫ x, (t - x) ∂μ[|Set.Iic t] := by
  have hne : μ.real (Set.Iic t) ≠ 0 := ENNReal.toReal_ne_zero.mpr ⟨h, measure_ne_top μ _⟩
  rw [mLPM_one_eq_real_mul_cond, mLPM_zero_eq_real, mul_div_cancel_left₀ _ hne]

theorem LPM_one_eq_LPM_zero_mul_condMean (t : ℝ) (x : Fin n → ℝ) (h : 0 < #{i | x i ≤ t}) :
    LPM 1 t x = LPM 0 t x * ((∑ i with x i ≤ t, (t - x i)) / #{i | x i ≤ t}) := by
  have hc : (#{i | x i ≤ t} : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr h.ne'
  have hs : ∑ i, lowerDev 1 t (x i) = ∑ i with x i ≤ t, (t - x i) := by
    rw [sum_filter]
    simp only [lowerDev, pow_one]
  rw [LPM_zero_eq_ecdf, LPM, hs, div_mul_div_comm, mul_comm, mul_div_mul_right _ _ hc]

theorem UPM_one_eq_UPM_zero_mul_condMean (t : ℝ) (x : Fin n → ℝ) (h : 0 < #{i | t < x i}) :
    UPM 1 t x = UPM 0 t x * ((∑ i with t < x i, (x i - t)) / #{i | t < x i}) := by
  have hc : (#{i | t < x i} : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr h.ne'
  have hs : ∑ i, upperDev 1 t (x i) = ∑ i with t < x i, (x i - t) := by
    rw [sum_filter]
    simp only [upperDev, pow_one]
  have h0 : UPM 0 t x = (#{i | t < x i} : ℝ) / n := by
    simp only [UPM, upperDev_zero]
    rw [sum_boole]
  rw [h0, UPM, hs, div_mul_div_comm, mul_comm, mul_div_mul_right _ _ hc]

theorem centralMoment_eq_mUPM_add_mLPM [IsProbabilityMeasure μ] {k : ℕ}
    (h : MemLp (fun x : ℝ => x) k μ) :
    ProbabilityTheory.centralMoment id k μ =
      mUPM k (∫ x, x ∂μ) μ + (-1) ^ k * mLPM k (∫ x, x ∂μ) μ := by
  rw [mUPM_add_neg_one_pow_mul_mLPM h]
  rfl

theorem rpow_mul_rpow_le {a b r s : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) (hr : 0 < r) (hrs : r ≤ s) :
    b ^ r * a ^ s ≤ b ^ s * a ^ r := by
  rcases ha.eq_or_lt with rfl | ha'
  · simp [Real.zero_rpow hr.ne', Real.zero_rpow (hr.trans_le hrs).ne']
  · have hb : 0 < b := ha'.trans_le hab
    have e : s = r + (s - r) := by ring
    have h1 : a ^ (s - r) ≤ b ^ (s - r) := Real.rpow_le_rpow ha hab (sub_nonneg.mpr hrs)
    have h2 : 0 ≤ b ^ r * a ^ r := by positivity
    calc b ^ r * a ^ s = b ^ r * a ^ r * a ^ (s - r) := by
          conv_lhs => rw [e, Real.rpow_add ha']
          ring
      _ ≤ b ^ r * a ^ r * b ^ (s - r) := mul_le_mul_of_nonneg_left h1 h2
      _ = b ^ s * a ^ r := by
          conv_rhs => rw [e, Real.rpow_add hb]
          ring

theorem rpow_mul_rpow_lt {a b r s : ℝ} (ha : 0 < a) (hab : a < b) (hrs : r < s) :
    b ^ r * a ^ s < b ^ s * a ^ r := by
  have hb : 0 < b := ha.trans hab
  have e : s = r + (s - r) := by ring
  have h1 : a ^ (s - r) < b ^ (s - r) := Real.rpow_lt_rpow ha.le hab (sub_pos.mpr hrs)
  have h2 : 0 < b ^ r * a ^ r := by positivity
  calc b ^ r * a ^ s = b ^ r * a ^ r * a ^ (s - r) := by
        conv_lhs => rw [e, Real.rpow_add ha]
        ring
    _ < b ^ r * a ^ r * b ^ (s - r) := mul_lt_mul_of_pos_left h1 h2
    _ = b ^ s * a ^ r := by
        conv_rhs => rw [e, Real.rpow_add hb]
        ring

theorem share_monotoneOn (d : Fin n → ℝ) (k : Fin n) (h0 : ∀ j, 0 ≤ d j)
    (hmax : ∀ j, d j ≤ d k) (hk : 0 < d k) :
    MonotoneOn (fun r : ℝ => d k ^ r / ∑ j, d j ^ r) (Set.Ioi 0) := by
  intro r hr s _ hrs
  have hS : ∀ u : ℝ, 0 < ∑ j, d j ^ u := fun u =>
    (Real.rpow_pos_of_pos hk u).trans_le
      (single_le_sum (fun j _ => Real.rpow_nonneg (h0 j) u) (mem_univ k))
  change d k ^ r / ∑ j, d j ^ r ≤ d k ^ s / ∑ j, d j ^ s
  rw [div_le_div_iff₀ (hS r) (hS s), mul_sum, mul_sum]
  exact sum_le_sum fun j _ => rpow_mul_rpow_le (h0 j) (hmax j) hr hrs

theorem share_strictMonoOn (d : Fin n → ℝ) (k : Fin n) (h0 : ∀ j, 0 ≤ d j)
    (hmax : ∀ j, d j ≤ d k) (hj : ∃ j, 0 < d j ∧ d j < d k) :
    StrictMonoOn (fun r : ℝ => d k ^ r / ∑ j, d j ^ r) (Set.Ioi 0) := by
  obtain ⟨j, hj0, hjk⟩ := hj
  have hk : 0 < d k := hj0.trans hjk
  intro r hr s _ hrs
  have hS : ∀ u : ℝ, 0 < ∑ j, d j ^ u := fun u =>
    (Real.rpow_pos_of_pos hk u).trans_le
      (single_le_sum (fun j _ => Real.rpow_nonneg (h0 j) u) (mem_univ k))
  change d k ^ r / ∑ j, d j ^ r < d k ^ s / ∑ j, d j ^ s
  rw [div_lt_div_iff₀ (hS r) (hS s), mul_sum, mul_sum]
  exact sum_lt_sum (fun i _ => rpow_mul_rpow_le (h0 i) (hmax i) hr hrs.le)
    ⟨j, mem_univ j, rpow_mul_rpow_lt hj0 hjk hrs⟩

def ch04Sample : Fin 6 → ℝ := ![-4, -2, -1, 1, 3, 5]

def ch07Sample : Fin 5 → ℝ := ![-2, -1, 0, 3, 5]

def ch08Sample : Fin 5 → ℝ := ![-3, -1, 0, 2, 4]

theorem ch04_deviations :
    (fun i => lowerDev 1 0 (ch04Sample i)) = ![4, 2, 1, 0, 0, 0] := by
  funext i
  fin_cases i <;> simp [ch04Sample, lowerDev]

theorem ch04_LPM_one : LPM 1 0 ch04Sample = 7 / 6 := by
  simp [LPM, lowerDev, ch04Sample, Fin.sum_univ_six]
  norm_num

theorem ch07_LPM_one : LPM 1 1 ch07Sample = 6 / 5 := by
  simp [LPM, lowerDev, ch07Sample, Fin.sum_univ_five]
  norm_num

theorem ch07_UPM_one : UPM 1 1 ch07Sample = 6 / 5 := by
  simp [UPM, upperDev, ch07Sample, Fin.sum_univ_five]
  norm_num

theorem ch07_LPM_zero : LPM 0 1 ch07Sample = 3 / 5 := by
  simp only [LPM, Fin.sum_univ_five]
  simp [lowerDev, ch07Sample]
  norm_num

theorem ch07_UPM_zero : UPM 0 1 ch07Sample = 2 / 5 := by
  simp only [UPM, Fin.sum_univ_five]
  simp [upperDev, ch07Sample]
  norm_num

theorem ch07_no_tie (i : Fin 5) : ch07Sample i ≠ 1 := by
  fin_cases i <;> simp [ch07Sample] <;> norm_num

theorem ch07_L1_ne_conditional :
    mLPM 1 1 (empMeasure ch07Sample) = 6 / 5 ∧
      ∫ v, (1 - v) ∂(empMeasure ch07Sample)[|Set.Iic 1] = 2 := by
  have := isProbabilityMeasure_empMeasure (by norm_num : 0 < 5) ch07Sample
  have h1 : mLPM 1 1 (empMeasure ch07Sample) = 6 / 5 := by
    rw [mLPM_empMeasure, ch07_LPM_one]
  have h0 : (empMeasure ch07Sample).real (Set.Iic 1) = 3 / 5 := by
    rw [← mLPM_zero_eq_real, mLPM_empMeasure, ch07_LPM_zero]
  have h := mLPM_one_eq_real_mul_cond (μ := empMeasure ch07Sample) 1
  rw [h1, h0] at h
  exact ⟨h1, by linarith⟩

theorem ch08_at_one :
    LPM 0 1 ch08Sample = 3 / 5 ∧ LPM 1 1 ch08Sample = 7 / 5 ∧ UPM 1 1 ch08Sample = 4 / 5 := by
  refine ⟨?_, ?_, ?_⟩
  · simp only [LPM, Fin.sum_univ_five]
    simp [lowerDev, ch08Sample]
    norm_num
  · simp [LPM, lowerDev, ch08Sample, Fin.sum_univ_five]
    norm_num
  · simp [UPM, upperDev, ch08Sample, Fin.sum_univ_five]
    norm_num

theorem ch08_no_tie_at_one (i : Fin 5) : ch08Sample i ≠ 1 := by
  fin_cases i <;> simp [ch08Sample] <;> norm_num

theorem ch08_at_zero :
    LPM 0 0 ch08Sample = 3 / 5 ∧ LPM 1 0 ch08Sample = 4 / 5 ∧ UPM 1 0 ch08Sample = 6 / 5 := by
  refine ⟨?_, ?_, ?_⟩
  · simp only [LPM, Fin.sum_univ_five]
    simp [lowerDev, ch08Sample]
    norm_num
  · simp [LPM, lowerDev, ch08Sample, Fin.sum_univ_five]
    norm_num
  · simp [UPM, upperDev, ch08Sample, Fin.sum_univ_five]
    norm_num

theorem ch08_comparisons :
    UPM 1 1 ch08Sample < LPM 1 1 ch08Sample ∧ LPM 1 0 ch08Sample < UPM 1 0 ch08Sample := by
  rw [ch08_at_one.2.1, ch08_at_one.2.2, ch08_at_zero.2.1, ch08_at_zero.2.2]
  constructor <;> norm_num

theorem ch08_mean : mean ch08Sample = 2 / 5 := by
  simp [mean, ch08Sample, Fin.sum_univ_five]
  norm_num

theorem ch08_UPM_one_sub_LPM_one (t : ℝ) :
    UPM 1 t ch08Sample - LPM 1 t ch08Sample = 2 / 5 - t := by
  rw [UPM_one_sub_LPM_one (by norm_num), ch08_mean]

theorem ch08_shortfall_dominates_iff (t : ℝ) :
    UPM 1 t ch08Sample < LPM 1 t ch08Sample ↔ 2 / 5 < t := by
  have h := ch08_UPM_one_sub_LPM_one t
  constructor <;> intro h' <;> linarith

end

end NNS.GapOperators
