-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.CoPartialMoments
public import Mathlib.Probability.CDF
public import Mathlib.Probability.Moments.Covariance
public import Mathlib.MeasureTheory.Integral.Prod
public import Mathlib.MeasureTheory.Integral.DominatedConvergence
public import Mathlib.Analysis.Calculus.Deriv.Slope

namespace NNS

open MeasureTheory ProbabilityTheory Set Filter Topology
open scoped ENNReal

@[expose] public noncomputable section

def mCoLPM (r s : ℕ) (tx ty : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ :=
  ∫ p, lowerDev r tx p.1 * lowerDev s ty p.2 ∂μ

def mCoUPM (r s : ℕ) (tx ty : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ :=
  ∫ p, upperDev r tx p.1 * upperDev s ty p.2 ∂μ

def mDLPM (dl du : ℕ) (tx ty : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ :=
  ∫ p, upperDev du tx p.1 * lowerDev dl ty p.2 ∂μ

def mDUPM (dl du : ℕ) (tx ty : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ :=
  ∫ p, lowerDev dl tx p.1 * upperDev du ty p.2 ∂μ

def pairEmp {n : ℕ} (x y : Fin n → ℝ) : Measure (ℝ × ℝ) :=
  (n : ℝ≥0∞)⁻¹ • ∑ i, Measure.dirac (x i, y i)

lemma measurable_lowerDev_co (r : ℕ) (t : ℝ) : Measurable (lowerDev r t) := by
  have : lowerDev r t = (Iic t).indicator fun v => (t - v) ^ r := by
    ext v; simp [lowerDev, indicator]
  rw [this]; exact (by fun_prop : Measurable _).indicator measurableSet_Iic

lemma measurable_upperDev_co (r : ℕ) (t : ℝ) : Measurable (upperDev r t) := by
  have : upperDev r t = (Ioi t).indicator fun v => (v - t) ^ r := by
    ext v; simp [upperDev, indicator]
  rw [this]; exact (by fun_prop : Measurable _).indicator measurableSet_Ioi

lemma lowerDev_zero_eq (t v : ℝ) : lowerDev 0 t v = (Iic t).indicator 1 v := by
  simp [lowerDev, indicator]

lemma upperDev_zero_eq (t v : ℝ) : upperDev 0 t v = (Ioi t).indicator 1 v := by
  simp [upperDev, indicator]

lemma indicator_mul_indicator_prod (s : Set ℝ) (u : Set ℝ) (p : ℝ × ℝ) :
    s.indicator (1 : ℝ → ℝ) p.1 * u.indicator 1 p.2 = (s ×ˢ u).indicator 1 p := by
  by_cases h1 : p.1 ∈ s <;> by_cases h2 : p.2 ∈ u <;> simp [indicator, h1, h2, mem_prod]

theorem mCoLPM_zero (tx ty : ℝ) (μ : Measure (ℝ × ℝ)) :
    mCoLPM 0 0 tx ty μ = μ.real (Iic tx ×ˢ Iic ty) := by
  simp_rw [mCoLPM, lowerDev_zero_eq, indicator_mul_indicator_prod]
  exact integral_indicator_one (measurableSet_Iic.prod measurableSet_Iic)

theorem mCoUPM_zero (tx ty : ℝ) (μ : Measure (ℝ × ℝ)) :
    mCoUPM 0 0 tx ty μ = μ.real (Ioi tx ×ˢ Ioi ty) := by
  simp_rw [mCoUPM, upperDev_zero_eq, indicator_mul_indicator_prod]
  exact integral_indicator_one (measurableSet_Ioi.prod measurableSet_Ioi)

theorem mDLPM_zero (tx ty : ℝ) (μ : Measure (ℝ × ℝ)) :
    mDLPM 0 0 tx ty μ = μ.real (Ioi tx ×ˢ Iic ty) := by
  simp_rw [mDLPM, upperDev_zero_eq, lowerDev_zero_eq, indicator_mul_indicator_prod]
  exact integral_indicator_one (measurableSet_Ioi.prod measurableSet_Iic)

theorem mDUPM_zero (tx ty : ℝ) (μ : Measure (ℝ × ℝ)) :
    mDUPM 0 0 tx ty μ = μ.real (Iic tx ×ˢ Ioi ty) := by
  simp_rw [mDUPM, upperDev_zero_eq, lowerDev_zero_eq, indicator_mul_indicator_prod]
  exact integral_indicator_one (measurableSet_Iic.prod measurableSet_Ioi)

theorem integral_pairEmp {n : ℕ} (x y : Fin n → ℝ) (f : ℝ × ℝ → ℝ) :
    ∫ p, f p ∂pairEmp x y = (∑ i, f (x i, y i)) / n := by
  rw [pairEmp, integral_smul_measure,
    integral_finsetSum_measure fun i _ => integrable_dirac enorm_lt_top]
  simp [integral_dirac, div_eq_inv_mul]

theorem mCoLPM_pairEmp {n : ℕ} (r s : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) :
    mCoLPM r s tx ty (pairEmp x y) = CoLPM r s tx ty x y := integral_pairEmp x y _

theorem mCoUPM_pairEmp {n : ℕ} (r s : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) :
    mCoUPM r s tx ty (pairEmp x y) = CoUPM r s tx ty x y := integral_pairEmp x y _

theorem mDLPM_pairEmp {n : ℕ} (dl du : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) :
    mDLPM dl du tx ty (pairEmp x y) = DLPM dl du tx ty x y := integral_pairEmp x y _

theorem mDUPM_pairEmp {n : ℕ} (dl du : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) :
    mDUPM dl du tx ty (pairEmp x y) = DUPM dl du tx ty x y := integral_pairEmp x y _

theorem isProbabilityMeasure_pairEmp {n : ℕ} (hn : 0 < n) (x y : Fin n → ℝ) :
    IsProbabilityMeasure (pairEmp x y) := by
  constructor
  simp only [pairEmp, Measure.smul_apply, Measure.coe_finsetSum, Finset.sum_apply,
    measure_univ, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one,
    smul_eq_mul]
  exact ENNReal.inv_mul_cancel (by exact_mod_cast hn.ne') (ENNReal.natCast_ne_top n)

def mLPMx (r : ℕ) (tx : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ := ∫ p, lowerDev r tx p.1 ∂μ

def mUPMx (r : ℕ) (tx : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ := ∫ p, upperDev r tx p.1 ∂μ

def mLPMy (s : ℕ) (ty : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ := ∫ p, lowerDev s ty p.2 ∂μ

def mUPMy (s : ℕ) (ty : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ := ∫ p, upperDev s ty p.2 ∂μ

abbrev xLe (tx : ℝ) : Set (ℝ × ℝ) := Iic tx ×ˢ univ

abbrev xGt (tx : ℝ) : Set (ℝ × ℝ) := Ioi tx ×ˢ univ

abbrev yLe (ty : ℝ) : Set (ℝ × ℝ) := univ ×ˢ Iic ty

abbrev yGt (ty : ℝ) : Set (ℝ × ℝ) := univ ×ˢ Ioi ty

theorem mLPMx_zero (tx : ℝ) (μ : Measure (ℝ × ℝ)) : mLPMx 0 tx μ = μ.real (xLe tx) := by
  have : (fun p : ℝ × ℝ => lowerDev 0 tx p.1) = (xLe tx).indicator 1 := by
    ext p; simp [lowerDev, indicator, mem_prod]
  rw [mLPMx, this]; exact integral_indicator_one (measurableSet_Iic.prod MeasurableSet.univ)

theorem mUPMx_zero (tx : ℝ) (μ : Measure (ℝ × ℝ)) : mUPMx 0 tx μ = μ.real (xGt tx) := by
  have : (fun p : ℝ × ℝ => upperDev 0 tx p.1) = (xGt tx).indicator 1 := by
    ext p; simp [upperDev, indicator, mem_prod]
  rw [mUPMx, this]; exact integral_indicator_one (measurableSet_Ioi.prod MeasurableSet.univ)

theorem mLPMy_zero (ty : ℝ) (μ : Measure (ℝ × ℝ)) : mLPMy 0 ty μ = μ.real (yLe ty) := by
  have : (fun p : ℝ × ℝ => lowerDev 0 ty p.2) = (yLe ty).indicator 1 := by
    ext p; simp [lowerDev, indicator, mem_prod]
  rw [mLPMy, this]; exact integral_indicator_one (MeasurableSet.univ.prod measurableSet_Iic)

theorem mUPMy_zero (ty : ℝ) (μ : Measure (ℝ × ℝ)) : mUPMy 0 ty μ = μ.real (yGt ty) := by
  have : (fun p : ℝ × ℝ => upperDev 0 ty p.2) = (yGt ty).indicator 1 := by
    ext p; simp [upperDev, indicator, mem_prod]
  rw [mUPMy, this]; exact integral_indicator_one (MeasurableSet.univ.prod measurableSet_Ioi)

theorem mLPMx_zero_eq_cdf (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] (tx : ℝ) :
    mLPMx 0 tx μ = cdf (μ.map Prod.fst) tx := by
  rw [mLPMx_zero, cdf_eq_real, map_measureReal_apply measurable_fst measurableSet_Iic, xLe,
    prod_univ]

theorem mLPMy_zero_eq_cdf (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] (ty : ℝ) :
    mLPMy 0 ty μ = cdf (μ.map Prod.snd) ty := by
  rw [mLPMy_zero, cdf_eq_real, map_measureReal_apply measurable_snd measurableSet_Iic, yLe,
    univ_prod]

theorem mLPMx_zero_add_mUPMx_zero (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] (tx : ℝ) :
    mLPMx 0 tx μ + mUPMx 0 tx μ = 1 := by
  rw [mLPMx_zero, mUPMx_zero, ← measureReal_union, ← union_prod, Iic_union_Ioi, univ_prod_univ,
    probReal_univ]
  · exact disjoint_prod.mpr (Or.inl (Iic_disjoint_Ioi le_rfl))
  · exact measurableSet_Ioi.prod MeasurableSet.univ

theorem mLPMy_zero_add_mUPMy_zero (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] (ty : ℝ) :
    mLPMy 0 ty μ + mUPMy 0 ty μ = 1 := by
  rw [mLPMy_zero, mUPMy_zero, ← measureReal_union, ← prod_union, Iic_union_Ioi, univ_prod_univ,
    probReal_univ]
  · exact disjoint_prod.mpr (Or.inr (Iic_disjoint_Ioi le_rfl))
  · exact MeasurableSet.univ.prod measurableSet_Ioi

section Sums

variable (μ : Measure (ℝ × ℝ)) [IsFiniteMeasure μ] (tx ty : ℝ)

lemma measureReal_prod_Iic_add_Ioi {a : Set ℝ} (ha : MeasurableSet a) (t : ℝ) :
    μ.real (a ×ˢ Iic t) + μ.real (a ×ˢ Ioi t) = μ.real (a ×ˢ univ) := by
  rw [← measureReal_union, ← prod_union, Iic_union_Ioi]
  · exact disjoint_prod.mpr (Or.inr (Iic_disjoint_Ioi le_rfl))
  · exact ha.prod measurableSet_Ioi

lemma measureReal_Iic_prod_add_Ioi {b : Set ℝ} (hb : MeasurableSet b) (t : ℝ) :
    μ.real (Iic t ×ˢ b) + μ.real (Ioi t ×ˢ b) = μ.real (univ ×ˢ b) := by
  rw [← measureReal_union, ← union_prod, Iic_union_Ioi]
  · exact disjoint_prod.mpr (Or.inl (Iic_disjoint_Ioi le_rfl))
  · exact measurableSet_Ioi.prod hb

theorem mCoLPM_zero_add_mDUPM_zero :
    mCoLPM 0 0 tx ty μ + mDUPM 0 0 tx ty μ = mLPMx 0 tx μ := by
  rw [mCoLPM_zero, mDUPM_zero, mLPMx_zero, measureReal_prod_Iic_add_Ioi μ measurableSet_Iic]

theorem mDLPM_zero_add_mCoUPM_zero :
    mDLPM 0 0 tx ty μ + mCoUPM 0 0 tx ty μ = mUPMx 0 tx μ := by
  rw [mDLPM_zero, mCoUPM_zero, mUPMx_zero, measureReal_prod_Iic_add_Ioi μ measurableSet_Ioi]

theorem mCoLPM_zero_add_mDLPM_zero :
    mCoLPM 0 0 tx ty μ + mDLPM 0 0 tx ty μ = mLPMy 0 ty μ := by
  rw [mCoLPM_zero, mDLPM_zero, mLPMy_zero, measureReal_Iic_prod_add_Ioi μ measurableSet_Iic]

end Sums

theorem quadrants_sum_eq_one_measure (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (tx ty : ℝ) :
    mCoLPM 0 0 tx ty μ + mCoUPM 0 0 tx ty μ + mDUPM 0 0 tx ty μ + mDLPM 0 0 tx ty μ = 1 := by
  have h1 := mCoLPM_zero_add_mDUPM_zero μ tx ty
  have h2 := mDLPM_zero_add_mCoUPM_zero μ tx ty
  have h3 := mLPMx_zero_add_mUPMx_zero μ tx
  linarith

lemma cond_real_eq {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω} {s : Set Ω}
    (hs : MeasurableSet s) (t : Set Ω) : (μ[|s]).real t = μ.real (s ∩ t) / μ.real s := by
  rw [measureReal_def, cond_apply hs, ENNReal.toReal_mul, ENNReal.toReal_inv, div_eq_inv_mul]
  rfl

lemma cond_real_mul {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω} [IsFiniteMeasure μ]
    {s : Set Ω} (hs : MeasurableSet s) (t : Set Ω) :
    (μ[|s]).real t * μ.real s = μ.real (s ∩ t) := by
  rw [measureReal_def, measureReal_def, measureReal_def, ← ENNReal.toReal_mul,
    cond_mul_eq_inter hs]

lemma xLe_inter_yLe (tx ty : ℝ) : xLe tx ∩ yLe ty = Iic tx ×ˢ Iic ty := by
  simp [prod_inter_prod]
lemma xGt_inter_yGt (tx ty : ℝ) : xGt tx ∩ yGt ty = Ioi tx ×ˢ Ioi ty := by
  simp [prod_inter_prod]
lemma xLe_inter_yGt (tx ty : ℝ) : xLe tx ∩ yGt ty = Iic tx ×ˢ Ioi ty := by
  simp [prod_inter_prod]
lemma xGt_inter_yLe (tx ty : ℝ) : xGt tx ∩ yLe ty = Ioi tx ×ˢ Iic ty := by
  simp [prod_inter_prod]

lemma measurableSet_xLe (tx : ℝ) : MeasurableSet (xLe tx) := measurableSet_Iic.prod .univ
lemma measurableSet_xGt (tx : ℝ) : MeasurableSet (xGt tx) := measurableSet_Ioi.prod .univ
lemma measurableSet_yLe (ty : ℝ) : MeasurableSet (yLe ty) :=
  MeasurableSet.univ.prod measurableSet_Iic
lemma measurableSet_yGt (ty : ℝ) : MeasurableSet (yGt ty) :=
  MeasurableSet.univ.prod measurableSet_Ioi

section Conditional

variable (μ : Measure (ℝ × ℝ)) (tx ty : ℝ)

theorem cond_yLe_xLe : (μ[|xLe tx]).real (yLe ty) = mCoLPM 0 0 tx ty μ / mLPMx 0 tx μ := by
  rw [cond_real_eq (measurableSet_xLe tx), xLe_inter_yLe, mCoLPM_zero, mLPMx_zero]

theorem cond_xLe_yLe : (μ[|yLe ty]).real (xLe tx) = mCoLPM 0 0 tx ty μ / mLPMy 0 ty μ := by
  rw [cond_real_eq (measurableSet_yLe ty), inter_comm, xLe_inter_yLe, mCoLPM_zero, mLPMy_zero]

theorem cond_yGt_xGt : (μ[|xGt tx]).real (yGt ty) = mCoUPM 0 0 tx ty μ / mUPMx 0 tx μ := by
  rw [cond_real_eq (measurableSet_xGt tx), xGt_inter_yGt, mCoUPM_zero, mUPMx_zero]

theorem cond_xGt_yGt : (μ[|yGt ty]).real (xGt tx) = mCoUPM 0 0 tx ty μ / mUPMy 0 ty μ := by
  rw [cond_real_eq (measurableSet_yGt ty), inter_comm, xGt_inter_yGt, mCoUPM_zero, mUPMy_zero]

theorem cond_yGt_xLe : (μ[|xLe tx]).real (yGt ty) = mDUPM 0 0 tx ty μ / mLPMx 0 tx μ := by
  rw [cond_real_eq (measurableSet_xLe tx), xLe_inter_yGt, mDUPM_zero, mLPMx_zero]

theorem cond_yLe_xGt : (μ[|xGt tx]).real (yLe ty) = mDLPM 0 0 tx ty μ / mUPMx 0 tx μ := by
  rw [cond_real_eq (measurableSet_xGt tx), xGt_inter_yLe, mDLPM_zero, mUPMx_zero]

theorem cond_xLe_yGt : (μ[|yGt ty]).real (xLe tx) = mDUPM 0 0 tx ty μ / mUPMy 0 ty μ := by
  rw [cond_real_eq (measurableSet_yGt ty), inter_comm, xLe_inter_yGt, mDUPM_zero, mUPMy_zero]

theorem cond_xGt_yLe : (μ[|yLe ty]).real (xGt tx) = mDLPM 0 0 tx ty μ / mLPMy 0 ty μ := by
  rw [cond_real_eq (measurableSet_yLe ty), inter_comm, xGt_inter_yLe, mDLPM_zero, mLPMy_zero]

variable [IsFiniteMeasure μ]

theorem bayes_lower :
    (μ[|xLe tx]).real (yLe ty) * mLPMx 0 tx μ = mCoLPM 0 0 tx ty μ ∧
      (μ[|yLe ty]).real (xLe tx) * mLPMy 0 ty μ = mCoLPM 0 0 tx ty μ := by
  rw [mLPMx_zero, mLPMy_zero, cond_real_mul (measurableSet_xLe tx),
    cond_real_mul (measurableSet_yLe ty), inter_comm (yLe ty), xLe_inter_yLe, mCoLPM_zero]
  exact ⟨rfl, rfl⟩

theorem bayes_upper :
    (μ[|xGt tx]).real (yGt ty) * mUPMx 0 tx μ = mCoUPM 0 0 tx ty μ ∧
      (μ[|yGt ty]).real (xGt tx) * mUPMy 0 ty μ = mCoUPM 0 0 tx ty μ := by
  rw [mUPMx_zero, mUPMy_zero, cond_real_mul (measurableSet_xGt tx),
    cond_real_mul (measurableSet_yGt ty), inter_comm (yGt ty), xGt_inter_yGt, mCoUPM_zero]
  exact ⟨rfl, rfl⟩

theorem bayes_DUPM :
    (μ[|xLe tx]).real (yGt ty) * mLPMx 0 tx μ = mDUPM 0 0 tx ty μ ∧
      (μ[|yGt ty]).real (xLe tx) * mUPMy 0 ty μ = mDUPM 0 0 tx ty μ := by
  rw [mLPMx_zero, mUPMy_zero, cond_real_mul (measurableSet_xLe tx),
    cond_real_mul (measurableSet_yGt ty), inter_comm (yGt ty), xLe_inter_yGt, mDUPM_zero]
  exact ⟨rfl, rfl⟩

theorem bayes_DLPM :
    (μ[|xGt tx]).real (yLe ty) * mUPMx 0 tx μ = mDLPM 0 0 tx ty μ ∧
      (μ[|yLe ty]).real (xGt tx) * mLPMy 0 ty μ = mDLPM 0 0 tx ty μ := by
  rw [mUPMx_zero, mLPMy_zero, cond_real_mul (measurableSet_xGt tx),
    cond_real_mul (measurableSet_yLe ty), inter_comm (yLe ty), xGt_inter_yLe, mDLPM_zero]
  exact ⟨rfl, rfl⟩

end Conditional

theorem mDLPM_zero_sub_mDUPM_zero (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ] (t : ℝ) :
    mDLPM 0 0 t t μ - mDUPM 0 0 t t μ = cdf (μ.map Prod.snd) t - cdf (μ.map Prod.fst) t := by
  rw [← mLPMx_zero_eq_cdf, ← mLPMy_zero_eq_cdf, ← mCoLPM_zero_add_mDLPM_zero μ t t,
    ← mCoLPM_zero_add_mDUPM_zero μ t t]
  ring

theorem mDLPM_zero_sub_mDUPM_zero_eq_upper (μ : Measure (ℝ × ℝ)) [IsProbabilityMeasure μ]
    (t : ℝ) : mDLPM 0 0 t t μ - mDUPM 0 0 t t μ = mUPMx 0 t μ - mUPMy 0 t μ := by
  have h1 := mCoLPM_zero_add_mDLPM_zero μ t t
  have h2 := mCoLPM_zero_add_mDUPM_zero μ t t
  have h3 := mLPMx_zero_add_mUPMx_zero μ t
  have h4 := mLPMy_zero_add_mUPMy_zero μ t
  linarith

lemma abs_lowerDev_one_le (t v : ℝ) : |lowerDev 1 t v| ≤ |v - t| := by
  unfold lowerDev
  split_ifs with h
  · rw [pow_one, abs_sub_comm]
  · simp

lemma abs_upperDev_one_le (t v : ℝ) : |upperDev 1 t v| ≤ |v - t| := by
  unfold upperDev
  split_ifs with h
  · rw [pow_one]
  · simp

lemma upperDev_one_add_lowerDev_one (t v : ℝ) : upperDev 1 t v + lowerDev 1 t v = |v - t| := by
  unfold upperDev lowerDev
  by_cases h : t < v
  · simp [h, not_le.mpr h, abs_of_pos (sub_pos.mpr h)]
  · have h' : v ≤ t := not_lt.mp h
    simp [h, h', abs_of_nonpos (sub_nonpos.mpr h')]

section DegreeOne

variable {μ : Measure (ℝ × ℝ)} {tx ty : ℝ}

lemma integrable_of_abs_le_coMoment {f : ℝ × ℝ → ℝ} (hf : Measurable f)
    (h : Integrable (fun p : ℝ × ℝ => (p.1 - tx) * (p.2 - ty)) μ)
    (hb : ∀ a b : ℝ, |f (a, b)| ≤ |a - tx| * |b - ty|) : Integrable f μ :=
  h.mono hf.aestronglyMeasurable (ae_of_all _ fun p => by
    simpa [Real.norm_eq_abs, abs_mul] using hb p.1 p.2)

lemma integrable_uu (h : Integrable (fun p : ℝ × ℝ => (p.1 - tx) * (p.2 - ty)) μ) :
    Integrable (fun p : ℝ × ℝ => upperDev 1 tx p.1 * upperDev 1 ty p.2) μ :=
  integrable_of_abs_le_coMoment
    (((measurable_upperDev_co 1 tx).comp measurable_fst).mul
      ((measurable_upperDev_co 1 ty).comp measurable_snd)) h fun a b => by
    rw [abs_mul]
    exact mul_le_mul (abs_upperDev_one_le _ _) (abs_upperDev_one_le _ _) (abs_nonneg _)
      (abs_nonneg _)

lemma integrable_ll (h : Integrable (fun p : ℝ × ℝ => (p.1 - tx) * (p.2 - ty)) μ) :
    Integrable (fun p : ℝ × ℝ => lowerDev 1 tx p.1 * lowerDev 1 ty p.2) μ :=
  integrable_of_abs_le_coMoment
    (((measurable_lowerDev_co 1 tx).comp measurable_fst).mul
      ((measurable_lowerDev_co 1 ty).comp measurable_snd)) h fun a b => by
    rw [abs_mul]
    exact mul_le_mul (abs_lowerDev_one_le _ _) (abs_lowerDev_one_le _ _) (abs_nonneg _)
      (abs_nonneg _)

lemma integrable_ul (h : Integrable (fun p : ℝ × ℝ => (p.1 - tx) * (p.2 - ty)) μ) :
    Integrable (fun p : ℝ × ℝ => upperDev 1 tx p.1 * lowerDev 1 ty p.2) μ :=
  integrable_of_abs_le_coMoment
    (((measurable_upperDev_co 1 tx).comp measurable_fst).mul
      ((measurable_lowerDev_co 1 ty).comp measurable_snd)) h fun a b => by
    rw [abs_mul]
    exact mul_le_mul (abs_upperDev_one_le _ _) (abs_lowerDev_one_le _ _) (abs_nonneg _)
      (abs_nonneg _)

lemma integrable_lu (h : Integrable (fun p : ℝ × ℝ => (p.1 - tx) * (p.2 - ty)) μ) :
    Integrable (fun p : ℝ × ℝ => lowerDev 1 tx p.1 * upperDev 1 ty p.2) μ :=
  integrable_of_abs_le_coMoment
    (((measurable_lowerDev_co 1 tx).comp measurable_fst).mul
      ((measurable_upperDev_co 1 ty).comp measurable_snd)) h fun a b => by
    rw [abs_mul]
    exact mul_le_mul (abs_lowerDev_one_le _ _) (abs_upperDev_one_le _ _) (abs_nonneg _)
      (abs_nonneg _)

theorem coMoment_decomposition_measure
    (h : Integrable (fun p : ℝ × ℝ => (p.1 - tx) * (p.2 - ty)) μ) :
    mCoUPM 1 1 tx ty μ + mCoLPM 1 1 tx ty μ - mDLPM 1 1 tx ty μ - mDUPM 1 1 tx ty μ
      = ∫ p, (p.1 - tx) * (p.2 - ty) ∂μ := by
  have hfun : (fun p : ℝ × ℝ => (p.1 - tx) * (p.2 - ty)) = fun p =>
      upperDev 1 tx p.1 * upperDev 1 ty p.2 + lowerDev 1 tx p.1 * lowerDev 1 ty p.2
        - upperDev 1 tx p.1 * lowerDev 1 ty p.2 - lowerDev 1 tx p.1 * upperDev 1 ty p.2 := by
    funext p
    rw [← upperDev_sub_lowerDev_one tx p.1, ← upperDev_sub_lowerDev_one ty p.2]
    ring
  rw [hfun, integral_sub, integral_sub, integral_add]
  · rfl
  all_goals first
    | exact integrable_uu h | exact integrable_ll h | exact integrable_ul h
    | exact integrable_lu h | exact (integrable_uu h).add (integrable_ll h)
    | exact ((integrable_uu h).add (integrable_ll h)).sub (integrable_ul h)

theorem covariance_decomposition_measure
    (h : Integrable (fun p : ℝ × ℝ => (p.1 - ∫ q, q.1 ∂μ) * (p.2 - ∫ q, q.2 ∂μ)) μ) :
    mCoUPM 1 1 (∫ q, q.1 ∂μ) (∫ q, q.2 ∂μ) μ + mCoLPM 1 1 (∫ q, q.1 ∂μ) (∫ q, q.2 ∂μ) μ
      - mDLPM 1 1 (∫ q, q.1 ∂μ) (∫ q, q.2 ∂μ) μ - mDUPM 1 1 (∫ q, q.1 ∂μ) (∫ q, q.2 ∂μ) μ
      = cov[Prod.fst, Prod.snd; μ] :=
  coMoment_decomposition_measure h

lemma integrable_coMoment_of_memLp [IsFiniteMeasure μ] (hX : MemLp Prod.fst 2 μ)
    (hY : MemLp Prod.snd 2 μ) (tx ty : ℝ) :
    Integrable (fun p : ℝ × ℝ => (p.1 - tx) * (p.2 - ty)) μ :=
  (hX.sub (memLp_const tx)).integrable_mul (hY.sub (memLp_const ty))

theorem hinge_total_eq (h : Integrable (fun p : ℝ × ℝ => (p.1 - tx) * (p.2 - ty)) μ) :
    mCoLPM 1 1 tx ty μ + mCoUPM 1 1 tx ty μ + mDLPM 1 1 tx ty μ + mDUPM 1 1 tx ty μ
      = ∫ p, |(p.1 - tx) * (p.2 - ty)| ∂μ := by
  have hfun : (fun p : ℝ × ℝ => |(p.1 - tx) * (p.2 - ty)|) = fun p =>
      lowerDev 1 tx p.1 * lowerDev 1 ty p.2 + upperDev 1 tx p.1 * upperDev 1 ty p.2
        + upperDev 1 tx p.1 * lowerDev 1 ty p.2 + lowerDev 1 tx p.1 * upperDev 1 ty p.2 := by
    funext p
    rw [abs_mul, ← upperDev_one_add_lowerDev_one tx p.1,
      ← upperDev_one_add_lowerDev_one ty p.2]
    ring
  rw [hfun, integral_add, integral_add, integral_add]
  · rfl
  all_goals first
    | exact integrable_uu h | exact integrable_ll h | exact integrable_ul h
    | exact integrable_lu h | exact (integrable_ll h).add (integrable_uu h)
    | exact ((integrable_ll h).add (integrable_uu h)).add (integrable_ul h)

theorem hinge_total_eq_zero_iff
    (h : Integrable (fun p : ℝ × ℝ => (p.1 - tx) * (p.2 - ty)) μ) :
    mCoLPM 1 1 tx ty μ + mCoUPM 1 1 tx ty μ + mDLPM 1 1 tx ty μ + mDUPM 1 1 tx ty μ = 0
      ↔ ∀ᵐ p ∂μ, p.1 = tx ∨ p.2 = ty := by
  rw [hinge_total_eq h, integral_eq_zero_iff_of_nonneg (fun p => abs_nonneg _) h.abs]
  refine Filter.eventually_congr (ae_of_all _ fun p => ?_)
  simp [mul_eq_zero, sub_eq_zero]

theorem hinge_weights {S : ℝ}
    (h : Integrable (fun p : ℝ × ℝ => (p.1 - tx) * (p.2 - ty)) μ)
    (hSdef : S = mCoLPM 1 1 tx ty μ + mCoUPM 1 1 tx ty μ + mDLPM 1 1 tx ty μ
      + mDUPM 1 1 tx ty μ) (hS : S ≠ 0) :
    S = ∫ p, |(p.1 - tx) * (p.2 - ty)| ∂μ ∧
      mCoLPM 1 1 tx ty μ / S + mCoUPM 1 1 tx ty μ / S + mDLPM 1 1 tx ty μ / S
        + mDUPM 1 1 tx ty μ / S = 1 ∧
      0 ≤ mCoLPM 1 1 tx ty μ ∧ 0 ≤ mCoUPM 1 1 tx ty μ ∧ 0 ≤ mDLPM 1 1 tx ty μ ∧
        0 ≤ mDUPM 1 1 tx ty μ := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · rw [hSdef, hinge_total_eq h]
  · rw [← add_div, ← add_div, ← add_div, ← hSdef, div_self hS]
  all_goals refine integral_nonneg fun p => ?_
  · exact mul_nonneg (lowerDev_nonneg _ _ _) (lowerDev_nonneg _ _ _)
  · exact mul_nonneg (upperDev_nonneg _ _ _) (upperDev_nonneg _ _ _)
  · exact mul_nonneg (upperDev_nonneg _ _ _) (lowerDev_nonneg _ _ _)
  · exact mul_nonneg (lowerDev_nonneg _ _ _) (upperDev_nonneg _ _ _)

end DegreeOne

section Independence

variable (r s : ℕ) (tx ty : ℝ) (ν₁ ν₂ : Measure ℝ) [SFinite ν₁] [SFinite ν₂]

theorem mCoLPM_prod :
    mCoLPM r s tx ty (ν₁.prod ν₂) = (∫ x, lowerDev r tx x ∂ν₁) * ∫ y, lowerDev s ty y ∂ν₂ :=
  integral_prod_mul _ _

theorem mCoUPM_prod :
    mCoUPM r s tx ty (ν₁.prod ν₂) = (∫ x, upperDev r tx x ∂ν₁) * ∫ y, upperDev s ty y ∂ν₂ :=
  integral_prod_mul _ _

theorem mDLPM_prod :
    mDLPM r s tx ty (ν₁.prod ν₂) = (∫ x, upperDev s tx x ∂ν₁) * ∫ y, lowerDev r ty y ∂ν₂ :=
  integral_prod_mul _ _

theorem mDUPM_prod :
    mDUPM r s tx ty (ν₁.prod ν₂) = (∫ x, lowerDev r tx x ∂ν₁) * ∫ y, upperDev s ty y ∂ν₂ :=
  integral_prod_mul _ _

end Independence

theorem mCoLPM_map_of_indepFun {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    [IsFiniteMeasure P] {X Y : Ω → ℝ} (hX : Measurable X) (hY : Measurable Y)
    (hXY : IndepFun X Y P) (r s : ℕ) (tx ty : ℝ) :
    mCoLPM r s tx ty (P.map fun ω => (X ω, Y ω))
      = (∫ ω, lowerDev r tx (X ω) ∂P) * ∫ ω, lowerDev s ty (Y ω) ∂P := by
  rw [(indepFun_iff_map_prod_eq_prod_map_map hX.aemeasurable hY.aemeasurable).mp hXY,
    mCoLPM_prod, integral_map hX.aemeasurable
      (measurable_lowerDev_co r tx).aestronglyMeasurable,
    integral_map hY.aemeasurable (measurable_lowerDev_co s ty).aestronglyMeasurable]

section Empirical

variable {n : ℕ} (x y : Fin n → ℝ)

theorem mLPMx_pairEmp (r : ℕ) (tx : ℝ) : mLPMx r tx (pairEmp x y) = LPM r tx x :=
  integral_pairEmp x y _

theorem mUPMx_pairEmp (r : ℕ) (tx : ℝ) : mUPMx r tx (pairEmp x y) = UPM r tx x :=
  integral_pairEmp x y _

theorem mLPMy_pairEmp (s : ℕ) (ty : ℝ) : mLPMy s ty (pairEmp x y) = LPM s ty y :=
  integral_pairEmp x y _

theorem mUPMy_pairEmp (s : ℕ) (ty : ℝ) : mUPMy s ty (pairEmp x y) = UPM s ty y :=
  integral_pairEmp x y _

end Empirical

lemma tendsto_hinge_slope (a x : ℝ) :
    Tendsto (fun a' => (lowerDev 1 a' x - lowerDev 1 a x) / (a' - a)) (𝓝[>] a)
      (𝓝 (lowerDev 0 a x)) := by
  by_cases hx : x ≤ a
  · refine tendsto_const_nhds.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with a' ha'
    have ha'x : x ≤ a' := hx.trans (le_of_lt ha')
    simp only [lowerDev, hx, ha'x, ↓reduceIte, pow_one, pow_zero]
    rw [show a' - x - (a - x) = a' - a by ring, div_self (sub_pos.mpr ha').ne']
  · refine tendsto_const_nhds.congr' ?_
    filter_upwards [Ioo_mem_nhdsGT (not_le.mp hx)] with a' ha'
    have ha'x : ¬ x ≤ a' := not_le.mpr ha'.2
    simp [lowerDev, hx, ha'x]

lemma abs_hinge_slope_le {a a' : ℝ} (h1 : a < a') (h2 : a' ≤ a + 1 / 2) (x : ℝ) :
    |(lowerDev 1 a' x - lowerDev 1 a x) / (a' - a)| ≤ 2 * lowerDev 1 (a + 1) x := by
  have hpos : 0 < a' - a := sub_pos.mpr h1
  rw [abs_div, abs_of_pos hpos, div_le_iff₀ hpos]
  simp only [lowerDev, pow_one]
  split_ifs <;> (rw [abs_of_nonneg (by linarith)]; nlinarith)

theorem hasDerivWithinAt_integral_hinge {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω}
    {X g : Ω → ℝ} (hX : Measurable X) (hg : AEStronglyMeasurable g μ)
    (hint : ∀ b, Integrable (fun ω => lowerDev 1 b (X ω) * g ω) μ) (a : ℝ) :
    HasDerivWithinAt (fun b => ∫ ω, lowerDev 1 b (X ω) * g ω ∂μ)
      (∫ ω, lowerDev 0 a (X ω) * g ω ∂μ) (Ici a) a := by
  rw [← hasDerivWithinAt_Ioi_iff_Ici, hasDerivWithinAt_iff_tendsto_slope,
    sdiff_singleton_eq_self (s := Ioi a) (lt_irrefl a)]
  have hslope : ∀ᶠ b in 𝓝[>] a, slope (fun b => ∫ ω, lowerDev 1 b (X ω) * g ω ∂μ) a b
      = ∫ ω, (lowerDev 1 b (X ω) - lowerDev 1 a (X ω)) / (b - a) * g ω ∂μ := by
    filter_upwards with b
    rw [slope_def_field, ← integral_sub (hint b) (hint a), ← integral_div]
    congr 1; funext ω; ring
  refine Tendsto.congr' (EventuallyEq.symm hslope) ?_
  refine tendsto_integral_filter_of_dominated_convergence
    (fun ω => 2 * |lowerDev 1 (a + 1) (X ω) * g ω|) ?_ ?_ ((hint (a + 1)).abs.const_mul 2) ?_
  · filter_upwards with b
    exact ((((measurable_lowerDev_co 1 b).comp hX).sub
      ((measurable_lowerDev_co 1 a).comp hX)).div_const _).aestronglyMeasurable.mul hg
  · filter_upwards [Ioo_mem_nhdsGT (show a < a + 1 / 2 by norm_num)] with b hb
    refine ae_of_all _ fun ω => ?_
    rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_of_nonneg (lowerDev_nonneg 1 (a + 1) (X ω)),
      ← mul_assoc]
    exact mul_le_mul_of_nonneg_right (abs_hinge_slope_le hb.1 hb.2.le (X ω)) (abs_nonneg _)
  · exact ae_of_all _ fun ω => (tendsto_hinge_slope a (X ω)).mul_const (g ω)

lemma lowerDev_zero_le_lowerDev_one_add_one (a x : ℝ) : lowerDev 0 a x ≤ lowerDev 1 (a + 1) x := by
  unfold lowerDev
  split_ifs <;> linarith

section Recovery

variable {μ : Measure (ℝ × ℝ)}

def HingeFinite (μ : Measure (ℝ × ℝ)) : Prop :=
  ∀ a b : ℝ, Integrable (fun p : ℝ × ℝ => lowerDev 1 a p.1 * lowerDev 1 b p.2) μ

theorem hingeFinite_of_integrable (μ : Measure (ℝ × ℝ))
    (h : Integrable (fun p : ℝ × ℝ => (1 + |p.1|) * (1 + |p.2|)) μ) : HingeFinite μ := by
  intro a b
  refine (h.const_mul ((1 + |a|) * (1 + |b|))).mono
    ((((measurable_lowerDev_co 1 a).comp measurable_fst).mul
      ((measurable_lowerDev_co 1 b).comp measurable_snd)).aestronglyMeasurable)
    (ae_of_all _ fun p => ?_)
  have hl : ∀ t v : ℝ, lowerDev 1 t v ≤ (1 + |t|) * (1 + |v|) := fun t v => by
    unfold lowerDev
    split_ifs
    · rw [pow_one]
      nlinarith [abs_nonneg t, abs_nonneg v, le_abs_self t, neg_abs_le v]
    · positivity
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg (mul_nonneg (lowerDev_nonneg _ _ _)
      (lowerDev_nonneg _ _ _)), abs_of_nonneg (by positivity)]
  calc lowerDev 1 a p.1 * lowerDev 1 b p.2
      ≤ ((1 + |a|) * (1 + |p.1|)) * ((1 + |b|) * (1 + |p.2|)) :=
        mul_le_mul (hl a p.1) (hl b p.2) (lowerDev_nonneg _ _ _) (by positivity)
    _ = (1 + |a|) * (1 + |b|) * ((1 + |p.1|) * (1 + |p.2|)) := by ring

theorem hasDerivWithinAt_hinge_x (hH : HingeFinite μ) (a b : ℝ) :
    HasDerivWithinAt (fun a' => mCoLPM 1 1 a' b μ) (mCoLPM 0 1 a b μ) (Ici a) a :=
  hasDerivWithinAt_integral_hinge measurable_fst
    ((measurable_lowerDev_co 1 b).comp measurable_snd).aestronglyMeasurable (fun a' => hH a' b) a

lemma integrable_mixed (hH : HingeFinite μ) (a b : ℝ) :
    Integrable (fun p : ℝ × ℝ => lowerDev 1 b p.2 * lowerDev 0 a p.1) μ :=
  (hH (a + 1) b).mono
    ((((measurable_lowerDev_co 1 b).comp measurable_snd).mul
      ((measurable_lowerDev_co 0 a).comp measurable_fst)).aestronglyMeasurable)
    (ae_of_all _ fun p => by
      rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg (mul_nonneg (lowerDev_nonneg _ _ _)
        (lowerDev_nonneg _ _ _)), abs_of_nonneg (mul_nonneg (lowerDev_nonneg _ _ _)
        (lowerDev_nonneg _ _ _)), mul_comm]
      exact mul_le_mul_of_nonneg_right (lowerDev_zero_le_lowerDev_one_add_one a p.1)
        (lowerDev_nonneg _ _ _))

theorem hasDerivWithinAt_hinge_y (hH : HingeFinite μ) (a b : ℝ) :
    HasDerivWithinAt (fun b' => mCoLPM 0 1 a b' μ) (mCoLPM 0 0 a b μ) (Ici b) b := by
  have h := hasDerivWithinAt_integral_hinge (μ := μ) measurable_snd
    ((measurable_lowerDev_co 0 a).comp measurable_fst).aestronglyMeasurable
    (integrable_mixed hH a) b
  simp only [mul_comm (lowerDev _ _ (Prod.snd _))] at h
  exact h

theorem hinge_recovers_joint_cdf (hH : HingeFinite μ) (a b : ℝ) :
    derivWithin (fun b' => derivWithin (fun a' => mCoLPM 1 1 a' b' μ) (Ici a) a) (Ici b) b
      = μ.real (Iic a ×ˢ Iic b) := by
  have hin : (fun b' => derivWithin (fun a' => mCoLPM 1 1 a' b' μ) (Ici a) a)
      = fun b' => mCoLPM 0 1 a b' μ := by
    funext b'
    exact (hasDerivWithinAt_hinge_x hH a b').derivWithin (uniqueDiffWithinAt_Ici a)
  rw [hin, (hasDerivWithinAt_hinge_y hH a b).derivWithin (uniqueDiffWithinAt_Ici b),
    mCoLPM_zero]

theorem hinge_mixed_deriv_eq (hH : HingeFinite μ) (a b : ℝ)
    (hin : ∀ᶠ b' in 𝓝 b, DifferentiableAt ℝ (fun a' => mCoLPM 1 1 a' b' μ) a)
    (hout : DifferentiableAt ℝ (fun b' => deriv (fun a' => mCoLPM 1 1 a' b' μ) a) b) :
    deriv (fun b' => deriv (fun a' => mCoLPM 1 1 a' b' μ) a) b = μ.real (Iic a ×ˢ Iic b) := by
  have heq : (fun b' => deriv (fun a' => mCoLPM 1 1 a' b' μ) a) =ᶠ[𝓝 b]
      fun b' => mCoLPM 0 1 a b' μ :=
    hin.mono fun b' hb' => (uniqueDiffWithinAt_Ici a).eq_deriv _
      hb'.hasDerivAt.hasDerivWithinAt (hasDerivWithinAt_hinge_x hH a b')
  rw [← mCoLPM_zero]
  exact (uniqueDiffWithinAt_Ici b).eq_deriv _
    (hout.hasDerivAt.congr_of_eventuallyEq heq.symm).hasDerivWithinAt
    (hasDerivWithinAt_hinge_y hH a b)

end Recovery

lemma isCountablySpanning_range_Iic : IsCountablySpanning (range (Iic : ℝ → Set ℝ)) := by
  refine ⟨fun n : ℕ => Iic (n : ℝ), fun n => mem_range_self _, ?_⟩
  exact eq_univ_of_forall fun x => mem_iUnion.mpr ⟨⌈x⌉₊, Nat.le_ceil x⟩

lemma generateFrom_Iic_prod_Iic :
    MeasurableSpace.generateFrom
        (image2 (· ×ˢ ·) (range (Iic : ℝ → Set ℝ)) (range (Iic : ℝ → Set ℝ)))
      = (inferInstance : MeasurableSpace (ℝ × ℝ)) := by
  have h : MeasurableSpace.generateFrom (range (Iic : ℝ → Set ℝ))
      = (inferInstance : MeasurableSpace ℝ) :=
    (borel_eq_generateFrom_Iic ℝ).symm.trans BorelSpace.measurable_eq.symm
  exact generateFrom_eq_prod h h isCountablySpanning_range_Iic isCountablySpanning_range_Iic

theorem ext_of_mCoLPM_zero {μ ν : Measure (ℝ × ℝ)} [IsProbabilityMeasure μ]
    [IsProbabilityMeasure ν] (h : ∀ a b, mCoLPM 0 0 a b μ = mCoLPM 0 0 a b ν) : μ = ν := by
  refine ext_of_generate_finite _ generateFrom_Iic_prod_Iic.symm
    (isPiSystem_Iic.prod isPiSystem_Iic) ?_ (by simp)
  rintro _ ⟨_, ⟨a, rfl⟩, _, ⟨b, rfl⟩, rfl⟩
  have hab := h a b
  rw [mCoLPM_zero, mCoLPM_zero, measureReal_def, measureReal_def] at hab
  exact (ENNReal.toReal_eq_toReal_iff' (measure_ne_top _ _) (measure_ne_top _ _)).mp hab

theorem ext_of_hinge {μ ν : Measure (ℝ × ℝ)} [IsProbabilityMeasure μ] [IsProbabilityMeasure ν]
    (hμ : HingeFinite μ) (hν : HingeFinite ν) (h : ∀ a b, mCoLPM 1 1 a b μ = mCoLPM 1 1 a b ν) :
    μ = ν := by
  refine ext_of_mCoLPM_zero fun a b => ?_
  rw [mCoLPM_zero, mCoLPM_zero, ← hinge_recovers_joint_cdf hμ, ← hinge_recovers_joint_cdf hν]
  simp only [h]

end

end NNS
