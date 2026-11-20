-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.CoConditional
public import NNSProofs.CoMomentsND
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
public import Mathlib.Analysis.Calculus.ParametricIntegral

namespace NNS.GapCoMoments

@[expose] public noncomputable section

section Pointwise

lemma lowerDev_eq_hinge_pow {d : ℕ} (hd : d ≠ 0) (t v : ℝ) :
    lowerDev d t v = max (t - v) 0 ^ d := by
  unfold lowerDev
  split_ifs with h
  · rw [max_eq_left (sub_nonneg.mpr h)]
  · rw [max_eq_right (sub_nonpos.mpr (not_le.mp h).le), zero_pow hd]

lemma upperDev_eq_hinge_pow {d : ℕ} (hd : d ≠ 0) (t v : ℝ) :
    upperDev d t v = max (v - t) 0 ^ d := by
  unfold upperDev
  split_ifs with h
  · rw [max_eq_left (sub_nonneg.mpr h.le)]
  · rw [max_eq_right (sub_nonpos.mpr (not_lt.mp h)), zero_pow hd]

lemma lowerDev_zero_ite (t v : ℝ) : lowerDev 0 t v = if v ≤ t then 1 else 0 := by
  simp [lowerDev]

lemma upperDev_zero_ite (t v : ℝ) : upperDev 0 t v = if t < v then 1 else 0 := by
  simp [upperDev]

lemma lowerDev_one_eq_hinge (t v : ℝ) : lowerDev 1 t v = max (t - v) 0 := by
  rw [lowerDev_eq_hinge_pow one_ne_zero, pow_one]

lemma upperDev_one_eq_hinge (t v : ℝ) : upperDev 1 t v = max (v - t) 0 := by
  rw [upperDev_eq_hinge_pow one_ne_zero, pow_one]

end Pointwise

section SampleConventions

open Finset

variable {n : ℕ}

lemma sum_ite_mul_ite (p q : Fin n → Prop) [DecidablePred p] [DecidablePred q] :
    (∑ i, (if p i then (1 : ℝ) else 0) * (if q i then 1 else 0)) = (#{i | p i ∧ q i} : ℝ) := by
  rw [← sum_boole]
  refine sum_congr rfl fun i _ => ?_
  by_cases hp : p i <;> by_cases hq : q i <;> simp [hp, hq]

theorem CoLPM_one_one_eq_hinge (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoLPM 1 1 tx ty x y = (∑ i, max (tx - x i) 0 * max (ty - y i) 0) / n := by
  simp only [CoLPM, lowerDev_one_eq_hinge]

theorem CoUPM_one_one_eq_hinge (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoUPM 1 1 tx ty x y = (∑ i, max (x i - tx) 0 * max (y i - ty) 0) / n := by
  simp only [CoUPM, upperDev_one_eq_hinge]

theorem DLPM_one_one_eq_hinge (tx ty : ℝ) (x y : Fin n → ℝ) :
    DLPM 1 1 tx ty x y = (∑ i, max (x i - tx) 0 * max (ty - y i) 0) / n := by
  simp only [DLPM, upperDev_one_eq_hinge, lowerDev_one_eq_hinge]

theorem DUPM_one_one_eq_hinge (tx ty : ℝ) (x y : Fin n → ℝ) :
    DUPM 1 1 tx ty x y = (∑ i, max (tx - x i) 0 * max (y i - ty) 0) / n := by
  simp only [DUPM, upperDev_one_eq_hinge, lowerDev_one_eq_hinge]

theorem CoLPM_zero_eq_card (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoLPM 0 0 tx ty x y = (#{i | x i ≤ tx ∧ y i ≤ ty} : ℝ) / n := by
  simp only [CoLPM, lowerDev_zero_ite, sum_ite_mul_ite]

theorem CoUPM_zero_eq_card (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoUPM 0 0 tx ty x y = (#{i | tx < x i ∧ ty < y i} : ℝ) / n := by
  simp only [CoUPM, upperDev_zero_ite, sum_ite_mul_ite]

theorem DLPM_zero_eq_card (tx ty : ℝ) (x y : Fin n → ℝ) :
    DLPM 0 0 tx ty x y = (#{i | tx < x i ∧ y i ≤ ty} : ℝ) / n := by
  simp only [DLPM, upperDev_zero_ite, lowerDev_zero_ite, sum_ite_mul_ite]

theorem DUPM_zero_eq_card (tx ty : ℝ) (x y : Fin n → ℝ) :
    DUPM 0 0 tx ty x y = (#{i | x i ≤ tx ∧ ty < y i} : ℝ) / n := by
  simp only [DUPM, upperDev_zero_ite, lowerDev_zero_ite, sum_ite_mul_ite]

end SampleConventions

section WorkedExample

def ch10X : Fin 5 → ℝ := ![-3, -1, 0, 2, 3]

def ch10Y : Fin 5 → ℝ := ![-2, -1, 1, 4, 5]

theorem ch10_worked_example :
    CoLPM 1 1 0 0 ch10X ch10Y = 7 / 5 ∧ CoUPM 1 1 0 0 ch10X ch10Y = 23 / 5 ∧
      DLPM 1 1 0 0 ch10X ch10Y = 0 ∧ DUPM 1 1 0 0 ch10X ch10Y = 0 := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · simp [CoLPM, lowerDev, ch10X, ch10Y, Fin.sum_univ_five]
    norm_num
  · simp [CoUPM, upperDev, ch10X, ch10Y, Fin.sum_univ_five]
    norm_num
  · simp [DLPM, lowerDev, upperDev, ch10X, ch10Y, Fin.sum_univ_five]
    norm_num
  · simp [DUPM, lowerDev, upperDev, ch10X, ch10Y, Fin.sum_univ_five]
    norm_num

theorem DLPM_degree_order :
    DLPM 1 2 (-2) 2 ch10X ch10Y = 7 / 5 ∧ DLPM 2 1 (-2) 2 ch10X ch10Y = 11 / 5 := by
  constructor
  · simp [DLPM, lowerDev, upperDev, ch10X, ch10Y, Fin.sum_univ_five]
    norm_num
  · simp [DLPM, lowerDev, upperDev, ch10X, ch10Y, Fin.sum_univ_five]
    norm_num

end WorkedExample

section PopulationConventions

open MeasureTheory

variable (μ : Measure (ℝ × ℝ)) (tx ty : ℝ)

theorem mCoLPM_one_one_eq_hinge :
    mCoLPM 1 1 tx ty μ = ∫ p, max (tx - p.1) 0 * max (ty - p.2) 0 ∂μ := by
  simp only [mCoLPM, lowerDev_one_eq_hinge]

theorem mCoUPM_one_one_eq_hinge :
    mCoUPM 1 1 tx ty μ = ∫ p, max (p.1 - tx) 0 * max (p.2 - ty) 0 ∂μ := by
  simp only [mCoUPM, upperDev_one_eq_hinge]

theorem mDLPM_one_one_eq_hinge :
    mDLPM 1 1 tx ty μ = ∫ p, max (p.1 - tx) 0 * max (ty - p.2) 0 ∂μ := by
  simp only [mDLPM, upperDev_one_eq_hinge, lowerDev_one_eq_hinge]

theorem mDUPM_one_one_eq_hinge :
    mDUPM 1 1 tx ty μ = ∫ p, max (tx - p.1) 0 * max (p.2 - ty) 0 ∂μ := by
  simp only [mDUPM, upperDev_one_eq_hinge, lowerDev_one_eq_hinge]

theorem hinge_pairEmp {n : ℕ} (x y : Fin n → ℝ) :
    mCoLPM 1 1 tx ty (pairEmp x y) = (∑ i, max (tx - x i) 0 * max (ty - y i) 0) / n := by
  rw [mCoLPM_pairEmp, CoLPM_one_one_eq_hinge]

end PopulationConventions

section Copula

open Finset

variable {n : ℕ}

theorem CoLPM_zero_pseudo_eq_empirical_copula (u v : ℝ) (x y : Fin n → ℝ) :
    CoLPM 0 0 u v (fun i => LPM 0 (x i) x) (fun i => LPM 0 (y i) y)
      = (#{i | (#{j | x j ≤ x i} : ℝ) / n ≤ u ∧ (#{j | y j ≤ y i} : ℝ) / n ≤ v} : ℝ) / n := by
  simp only [CoLPM_zero_eq_card, LPM_zero_eq_ecdf]

lemma pseudo_le_one (x : Fin n → ℝ) (i : Fin n) : LPM 0 (x i) x ≤ 1 := by
  rw [LPM_zero_eq_ecdf]
  have hn : (0 : ℝ) < n := by exact_mod_cast Fin.pos i
  rw [div_le_one hn]
  exact_mod_cast (card_le_univ _).trans_eq (Fintype.card_fin n)

lemma one_div_le_pseudo (x : Fin n → ℝ) (i : Fin n) : 1 / (n : ℝ) ≤ LPM 0 (x i) x := by
  rw [LPM_zero_eq_ecdf]
  have hn : (0 : ℝ) < n := by exact_mod_cast Fin.pos i
  have hc : 1 ≤ #{j | x j ≤ x i} := card_pos.mpr ⟨i, by simp⟩
  gcongr
  exact_mod_cast hc

theorem empirical_copula_of_one_le {u : ℝ} (hu : 1 ≤ u) (v : ℝ) (x y : Fin n → ℝ) :
    CoLPM 0 0 u v (fun i => LPM 0 (x i) x) (fun i => LPM 0 (y i) y)
      = LPM 0 v fun i => LPM 0 (y i) y := by
  have h : (univ.filter fun i => LPM 0 (x i) x ≤ u ∧ LPM 0 (y i) y ≤ v)
      = univ.filter fun i => LPM 0 (y i) y ≤ v := by
    ext i
    simp only [mem_filter, mem_univ, true_and, and_iff_right_iff_imp]
    exact fun _ => (pseudo_le_one x i).trans hu
  rw [CoLPM_zero_eq_card, LPM_zero_eq_ecdf, h]

theorem empirical_copula_of_lt {u : ℝ} (hu : u < 1 / (n : ℝ)) (v : ℝ) (x y : Fin n → ℝ) :
    CoLPM 0 0 u v (fun i => LPM 0 (x i) x) (fun i => LPM 0 (y i) y) = 0 := by
  have h0 : #{i | LPM 0 (x i) x ≤ u ∧ LPM 0 (y i) y ≤ v} = 0 := by
    rw [card_eq_zero, filter_eq_empty_iff]
    exact fun i _ h => absurd (hu.trans_le (one_div_le_pseudo x i)) (not_lt.mpr h.1)
  rw [CoLPM_zero_eq_card, h0, Nat.cast_zero, zero_div]

theorem empirical_copula_one_one (hn : 0 < n) (x y : Fin n → ℝ) :
    CoLPM 0 0 1 1 (fun i => LPM 0 (x i) x) (fun i => LPM 0 (y i) y) = 1 := by
  rw [CoLPM_zero_eq_card, filter_true_of_mem, card_univ, Fintype.card_fin,
    div_self (Nat.cast_ne_zero.mpr hn.ne')]
  exact fun i _ => ⟨pseudo_le_one x i, pseudo_le_one y i⟩

end Copula

section Completeness

open MeasureTheory

theorem coPartialMoments_eq_of_degree_zero_eq {μ ν : Measure (ℝ × ℝ)} [IsProbabilityMeasure μ]
    [IsProbabilityMeasure ν] (h : ∀ a b, mCoLPM 0 0 a b μ = mCoLPM 0 0 a b ν) (r s : ℕ)
    (tx ty : ℝ) :
    mCoLPM r s tx ty μ = mCoLPM r s tx ty ν ∧ mCoUPM r s tx ty μ = mCoUPM r s tx ty ν ∧
      mDLPM r s tx ty μ = mDLPM r s tx ty ν ∧ mDUPM r s tx ty μ = mDUPM r s tx ty ν := by
  obtain rfl := ext_of_mCoLPM_zero h
  exact ⟨rfl, rfl, rfl, rfl⟩

end Completeness

section RealDegreeCo

open Finset

variable {n : ℕ}

def rCoLPM (dx dy tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, rLowerDev dx tx (x i) * rLowerDev dy ty (y i)) / n

def rCoUPM (dx dy tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, rUpperDev dx tx (x i) * rUpperDev dy ty (y i)) / n

def rDLPM (dLower dUpper tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, rUpperDev dUpper tx (x i) * rLowerDev dLower ty (y i)) / n

def rDUPM (dLower dUpper tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, rLowerDev dLower tx (x i) * rUpperDev dUpper ty (y i)) / n

theorem rCo_natCast (a b : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) :
    rCoLPM a b tx ty x y = CoLPM a b tx ty x y ∧ rCoUPM a b tx ty x y = CoUPM a b tx ty x y ∧
      rDLPM a b tx ty x y = DLPM a b tx ty x y ∧ rDUPM a b tx ty x y = DUPM a b tx ty x y := by
  simp only [rCoLPM, rCoUPM, rDLPM, rDUPM, CoLPM, CoUPM, DLPM, DUPM, rLowerDev_natCast,
    rUpperDev_natCast, and_self]

def lowKernel (d t v : ℝ) : ℝ := if d = 0 then (if v ≤ t then 1 else 0) else max (t - v) 0 ^ d

def upKernel (d t v : ℝ) : ℝ := if d = 0 then (if t < v then 1 else 0) else max (v - t) 0 ^ d

lemma rLowerDev_eq_lowKernel (d t v : ℝ) : rLowerDev d t v = lowKernel d t v := by
  by_cases hd : d = 0
  · simp [lowKernel, hd, rLowerDev_zero]
  · simp [lowKernel, hd, rLowerDev_eq_max hd]

lemma rUpperDev_eq_upKernel (d t v : ℝ) : rUpperDev d t v = upKernel d t v := by
  by_cases hd : d = 0
  · simp [upKernel, hd, rUpperDev_zero]
  · simp [upKernel, hd, rUpperDev_eq_max hd]

theorem rCo_eq_kernels (dx dy dl du tx ty : ℝ) (x y : Fin n → ℝ) :
    rCoLPM dx dy tx ty x y = (∑ i, lowKernel dx tx (x i) * lowKernel dy ty (y i)) / n ∧
    rCoUPM dx dy tx ty x y = (∑ i, upKernel dx tx (x i) * upKernel dy ty (y i)) / n ∧
    rDLPM dl du tx ty x y = (∑ i, upKernel du tx (x i) * lowKernel dl ty (y i)) / n ∧
    rDUPM dl du tx ty x y = (∑ i, lowKernel dl tx (x i) * upKernel du ty (y i)) / n := by
  simp only [rCoLPM, rCoUPM, rDLPM, rDUPM, rLowerDev_eq_lowKernel, rUpperDev_eq_upKernel,
    and_self]

end RealDegreeCo

section FiniteDifference

open Finset

variable {n : ℕ}

lemma lowerDev_one_add_sub {a h v : ℝ} (hh : 0 ≤ h) (hv : v ∉ Set.Ioo a (a + h)) :
    lowerDev 1 (a + h) v - lowerDev 1 a v = h * lowerDev 0 a v := by
  simp only [Set.mem_Ioo, not_and, not_lt] at hv
  unfold lowerDev
  by_cases h1 : v ≤ a
  · rw [ite_eq_left (h1.trans (le_add_of_nonneg_right hh)), ite_eq_left h1, ite_eq_left h1]
    ring
  · rw [ite_eq_right h1, ite_eq_right h1]
    have h2 := hv (not_le.mp h1)
    split_ifs with h3
    · rw [le_antisymm h3 h2]
      ring
    · ring

theorem hinge_finite_difference {a b h : ℝ} (hh : 0 < h) (x y : Fin n → ℝ)
    (hx : ∀ i, x i ∉ Set.Ioo a (a + h)) (hy : ∀ i, y i ∉ Set.Ioo b (b + h)) :
    (CoLPM 1 1 (a + h) (b + h) x y - CoLPM 1 1 (a + h) b x y - CoLPM 1 1 a (b + h) x y
      + CoLPM 1 1 a b x y) / h ^ 2 = CoLPM 0 0 a b x y := by
  have key : ∀ i, lowerDev 1 (a + h) (x i) * lowerDev 1 (b + h) (y i)
      - lowerDev 1 (a + h) (x i) * lowerDev 1 b (y i)
      - lowerDev 1 a (x i) * lowerDev 1 (b + h) (y i) + lowerDev 1 a (x i) * lowerDev 1 b (y i)
      = h ^ 2 * (lowerDev 0 a (x i) * lowerDev 0 b (y i)) := fun i => by
    have e1 := lowerDev_one_add_sub hh.le (hx i)
    have e2 := lowerDev_one_add_sub hh.le (hy i)
    linear_combination (lowerDev 1 (b + h) (y i) - lowerDev 1 b (y i)) * e1
      + h * lowerDev 0 a (x i) * e2
  simp only [CoLPM]
  rw [← sub_div, ← sub_div, ← add_div, ← sum_sub_distrib, ← sum_sub_distrib,
    ← sum_add_distrib, sum_congr rfl fun i _ => key i, ← mul_sum]
  field_simp

open Filter Topology in
lemma eventually_no_observation_inside (a : ℝ) (x : Fin n → ℝ) :
    ∀ᶠ h in 𝓝[>] (0 : ℝ), ∀ i, x i ∉ Set.Ioo a (a + h) := by
  rw [Filter.eventually_all]
  intro i
  by_cases hi : a < x i
  · filter_upwards [(eventually_lt_nhds (sub_pos.mpr hi)).filter_mono nhdsWithin_le_nhds]
      with h hlt hp
    linarith [hp.2]
  · exact Filter.Eventually.of_forall fun _ hp => hi hp.1

open Filter Topology in
theorem eventually_hinge_finite_difference (a b : ℝ) (x y : Fin n → ℝ) :
    ∀ᶠ h in 𝓝[>] (0 : ℝ), (CoLPM 1 1 (a + h) (b + h) x y - CoLPM 1 1 (a + h) b x y -
      CoLPM 1 1 a (b + h) x y + CoLPM 1 1 a b x y) / h ^ 2 = CoLPM 0 0 a b x y := by
  filter_upwards [self_mem_nhdsWithin, eventually_no_observation_inside a x,
    eventually_no_observation_inside b y] with h hh hx hy
  exact hinge_finite_difference hh x y hx hy

end FiniteDifference

section DegreeZeroLimit

open MeasureTheory Filter Topology Set

variable {μ : Measure (ℝ × ℝ)}

def quadLL (q a b : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ := ∫ p, rLowerDev q a p.1 * rLowerDev q b p.2 ∂μ

def quadUU (q a b : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ := ∫ p, rUpperDev q a p.1 * rUpperDev q b p.2 ∂μ

def quadUL (q a b : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ := ∫ p, rUpperDev q a p.1 * rLowerDev q b p.2 ∂μ

def quadLU (q a b : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ := ∫ p, rLowerDev q a p.1 * rUpperDev q b p.2 ∂μ

def quadTotal (q a b : ℝ) (μ : Measure (ℝ × ℝ)) : ℝ :=
  quadLL q a b μ + quadUU q a b μ + quadUL q a b μ + quadLU q a b μ

lemma measurable_rLowerDev' (q a : ℝ) : Measurable (rLowerDev q a) :=
  Measurable.ite (measurableSet_le measurable_id measurable_const)
    ((measurable_const.sub measurable_id).pow_const q) measurable_const

lemma measurable_rUpperDev' (q a : ℝ) : Measurable (rUpperDev q a) :=
  Measurable.ite (measurableSet_lt measurable_const measurable_id)
    ((measurable_id.sub measurable_const).pow_const q) measurable_const

lemma tendsto_rLowerDev_zero (a x : ℝ) :
    Tendsto (fun q => rLowerDev q a x) (𝓝[>] 0) (𝓝 (if x < a then 1 else 0)) := by
  rcases lt_trichotomy x a with h | rfl | h
  · simp only [rLowerDev, h.le, h, ↓reduceIte]
    have := (Real.continuousAt_const_rpow (a := a - x) (b := 0) (sub_pos.mpr h).ne').tendsto
    rw [Real.rpow_zero] at this
    exact this.mono_left nhdsWithin_le_nhds
  · simp only [rLowerDev, le_refl, sub_self, lt_irrefl, ↓reduceIte]
    refine tendsto_const_nhds.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with q hq
    rw [Real.zero_rpow (ne_of_gt hq)]
  · simp [rLowerDev, not_le.mpr h, not_lt.mpr h.le]

lemma tendsto_rUpperDev_zero (a x : ℝ) :
    Tendsto (fun q => rUpperDev q a x) (𝓝[>] 0) (𝓝 (if a < x then 1 else 0)) := by
  rcases lt_trichotomy a x with h | rfl | h
  · simp only [rUpperDev, h, ↓reduceIte]
    have := (Real.continuousAt_const_rpow (a := x - a) (b := 0) (sub_pos.mpr h).ne').tendsto
    rw [Real.rpow_zero] at this
    exact this.mono_left nhdsWithin_le_nhds
  · simp [rUpperDev]
  · simp [rUpperDev, not_lt.mpr h.le]

lemma ite_mul_ite_eq_indicator (P Q : ℝ → Prop) [DecidablePred P] [DecidablePred Q] :
    (fun p : ℝ × ℝ => (if P p.1 then (1 : ℝ) else 0) * (if Q p.2 then 1 else 0)) =
      {p : ℝ × ℝ | P p.1 ∧ Q p.2}.indicator 1 := by
  ext p
  by_cases h1 : P p.1 <;> by_cases h2 : Q p.2 <;> simp [Set.indicator, h1, h2]

lemma integral_ite_mul_ite (P Q : ℝ → Prop) [DecidablePred P] [DecidablePred Q]
    (hP : MeasurableSet {x | P x}) (hQ : MeasurableSet {y | Q y}) :
    ∫ p, (if P p.1 then (1 : ℝ) else 0) * (if Q p.2 then 1 else 0) ∂μ =
      μ.real {p : ℝ × ℝ | P p.1 ∧ Q p.2} := by
  rw [ite_mul_ite_eq_indicator]
  exact integral_indicator_one ((hP.preimage measurable_fst).inter (hQ.preimage measurable_snd))

lemma integrable_ite_mul_ite [IsFiniteMeasure μ] (P Q : ℝ → Prop) [DecidablePred P]
    [DecidablePred Q] (hP : MeasurableSet {x | P x}) (hQ : MeasurableSet {y | Q y}) :
    Integrable (fun p : ℝ × ℝ => (if P p.1 then (1 : ℝ) else 0) * (if Q p.2 then 1 else 0)) μ := by
  rw [ite_mul_ite_eq_indicator]
  exact (integrable_const (1 : ℝ)).indicator
    ((hP.preimage measurable_fst).inter (hQ.preimage measurable_snd))

lemma measureReal_off_axis [IsFiniteMeasure μ] (a b : ℝ) :
    μ.real {p : ℝ × ℝ | p.1 < a ∧ p.2 < b} + μ.real {p : ℝ × ℝ | a < p.1 ∧ b < p.2}
      + μ.real {p : ℝ × ℝ | a < p.1 ∧ p.2 < b} + μ.real {p : ℝ × ℝ | p.1 < a ∧ b < p.2}
      = μ.real {p : ℝ × ℝ | p.1 ≠ a ∧ p.2 ≠ b} := by
  have hl : MeasurableSet {x : ℝ | x < a} := measurableSet_Iio
  have hu : MeasurableSet {x : ℝ | a < x} := measurableSet_Ioi
  have hl' : MeasurableSet {y : ℝ | y < b} := measurableSet_Iio
  have hu' : MeasurableSet {y : ℝ | b < y} := measurableSet_Ioi
  have hne : MeasurableSet {x : ℝ | x ≠ a} := (measurableSet_singleton a).compl
  have hne' : MeasurableSet {y : ℝ | y ≠ b} := (measurableSet_singleton b).compl
  have i1 := integrable_ite_mul_ite (μ := μ) _ _ hl hl'
  have i2 := integrable_ite_mul_ite (μ := μ) _ _ hu hu'
  have i3 := integrable_ite_mul_ite (μ := μ) _ _ hu hl'
  have i4 := integrable_ite_mul_ite (μ := μ) _ _ hl hu'
  rw [← integral_ite_mul_ite _ _ hl hl', ← integral_ite_mul_ite _ _ hu hu',
    ← integral_ite_mul_ite _ _ hu hl', ← integral_ite_mul_ite _ _ hl hu',
    ← integral_ite_mul_ite _ _ hne hne']
  rw [← integral_add i1 i2]
  rw [← integral_add (f := fun p => _ + _) (by exact i1.add i2) i3]
  rw [← integral_add (f := fun p => _ + _ + _) (by exact (i1.add i2).add i3) i4]
  congr 1
  ext ⟨x, y⟩
  rcases lt_trichotomy x a with h1 | rfl | h1 <;> rcases lt_trichotomy y b with h2 | rfl | h2
  · simp [h1, h2, h1.ne, h2.ne, h1.not_gt, h2.not_gt]
  · simp [h1]
  · simp [h1, h2, h1.ne, h2.ne', h1.not_gt, h2.not_gt]
  · simp
  · simp
  · simp
  · simp [h1, h2, h1.ne', h2.ne, h1.not_gt, h2.not_gt]
  · simp [h1]
  · simp [h1, h2, h1.ne', h2.ne', h1.not_gt, h2.not_gt]

lemma measureReal_eq_of_subset_union_null {S T N : Set (ℝ × ℝ)} [IsFiniteMeasure μ]
    (hST : S ⊆ T) (hTS : T ⊆ S ∪ N) (hN : μ N = 0) : μ.real S = μ.real T := by
  refine le_antisymm (measureReal_mono hST) ?_
  calc μ.real T ≤ μ.real (S ∪ N) := measureReal_mono hTS
    _ ≤ μ.real S + μ.real N := measureReal_union_le _ _
    _ = μ.real S := by simp [measureReal_def, hN]

lemma quad_sum_eq_rpow (q a b : ℝ) (p : ℝ × ℝ) :
    rLowerDev q a p.1 * rLowerDev q b p.2 + rUpperDev q a p.1 * rUpperDev q b p.2 +
      rUpperDev q a p.1 * rLowerDev q b p.2 + rLowerDev q a p.1 * rUpperDev q b p.2 =
      (|p.1 - a| * |p.2 - b|) ^ q := by
  have h1 := rLowerDev_add_rUpperDev q a p.1
  have h2 := rLowerDev_add_rUpperDev q b p.2
  rw [abs_sub_comm] at h1 h2
  rw [Real.mul_rpow (abs_nonneg _) (abs_nonneg _), ← h1, ← h2]
  ring

theorem integrable_rpow_iff_quad (q a b : ℝ) :
    Integrable (fun p : ℝ × ℝ => (|p.1 - a| * |p.2 - b|) ^ q) μ ↔
      Integrable (fun p : ℝ × ℝ => rLowerDev q a p.1 * rLowerDev q b p.2) μ ∧
      Integrable (fun p : ℝ × ℝ => rUpperDev q a p.1 * rUpperDev q b p.2) μ ∧
      Integrable (fun p : ℝ × ℝ => rUpperDev q a p.1 * rLowerDev q b p.2) μ ∧
      Integrable (fun p : ℝ × ℝ => rLowerDev q a p.1 * rUpperDev q b p.2) μ := by
  have mL := measurable_rLowerDev' q
  have mU := measurable_rUpperDev' q
  have hL := rLowerDev_nonneg q
  have hU := rUpperDev_nonneg q
  constructor
  · intro h
    have dom : ∀ {g : ℝ × ℝ → ℝ}, Measurable g → (∀ p, 0 ≤ g p) →
        (∀ p, g p ≤ (|p.1 - a| * |p.2 - b|) ^ q) → Integrable g μ := fun hg h0 hle =>
      h.mono' hg.aestronglyMeasurable (ae_of_all _ fun p => by
        rw [Real.norm_of_nonneg (h0 p)]; exact hle p)
    refine ⟨dom (((mL a).comp measurable_fst).mul ((mL b).comp measurable_snd))
        (fun p => mul_nonneg (hL a _) (hL b _)) fun p => ?_,
      dom (((mU a).comp measurable_fst).mul ((mU b).comp measurable_snd))
        (fun p => mul_nonneg (hU a _) (hU b _)) fun p => ?_,
      dom (((mU a).comp measurable_fst).mul ((mL b).comp measurable_snd))
        (fun p => mul_nonneg (hU a _) (hL b _)) fun p => ?_,
      dom (((mL a).comp measurable_fst).mul ((mU b).comp measurable_snd))
        (fun p => mul_nonneg (hL a _) (hU b _)) fun p => ?_⟩ <;>
    · rw [← quad_sum_eq_rpow q a b p]
      nlinarith [mul_nonneg (hL a p.1) (hL b p.2), mul_nonneg (hU a p.1) (hU b p.2),
        mul_nonneg (hU a p.1) (hL b p.2), mul_nonneg (hL a p.1) (hU b p.2)]
  · rintro ⟨h1, h2, h3, h4⟩
    have := ((h1.add h2).add h3).add h4
    refine this.congr (ae_of_all _ fun p => ?_)
    exact quad_sum_eq_rpow q a b p

lemma rpow_le_one_add_rpow {v q ε : ℝ} (hv : 0 ≤ v) (hq : 0 ≤ q) (hqε : q ≤ ε) :
    v ^ q ≤ 1 + v ^ ε := by
  rcases le_total v 1 with h | h
  · linarith [Real.rpow_le_one hv h hq, Real.rpow_nonneg hv ε]
  · linarith [Real.rpow_le_rpow_of_exponent_le h hqε]

lemma rLowerDev_le_abs_rpow (q a x : ℝ) : rLowerDev q a x ≤ |x - a| ^ q := by
  have := rLowerDev_add_rUpperDev q a x
  rw [abs_sub_comm] at this
  linarith [rUpperDev_nonneg q a x]

lemma rUpperDev_le_abs_rpow (q a x : ℝ) : rUpperDev q a x ≤ |x - a| ^ q := by
  have := rLowerDev_add_rUpperDev q a x
  rw [abs_sub_comm] at this
  linarith [rLowerDev_nonneg q a x]

lemma tendsto_integral_kernel_mul_rpow [IsFiniteMeasure μ] {K₁ K₂ : ℝ → ℝ → ℝ} {g₁ g₂ : ℝ → ℝ}
    {a b ε : ℝ} (hε : 0 < ε)
    (hm₁ : ∀ q, Measurable (K₁ q)) (hm₂ : ∀ q, Measurable (K₂ q))
    (hK₁ : ∀ x, Tendsto (fun q => K₁ q x) (𝓝[>] 0) (𝓝 (g₁ x)))
    (hK₂ : ∀ y, Tendsto (fun q => K₂ q y) (𝓝[>] 0) (𝓝 (g₂ y)))
    (hb₁ : ∀ q ∈ Ioo (0 : ℝ) ε, ∀ x, 0 ≤ K₁ q x ∧ K₁ q x ≤ |x - a| ^ q)
    (hb₂ : ∀ q ∈ Ioo (0 : ℝ) ε, ∀ y, 0 ≤ K₂ q y ∧ K₂ q y ≤ |y - b| ^ q)
    (hμ : Integrable (fun p : ℝ × ℝ => (|p.1 - a| * |p.2 - b|) ^ ε) μ) :
    Tendsto (fun q => ∫ p, K₁ q p.1 * K₂ q p.2 ∂μ) (𝓝[>] 0)
      (𝓝 (∫ p, g₁ p.1 * g₂ p.2 ∂μ)) := by
  refine tendsto_integral_filter_of_dominated_convergence
    (fun p => 1 + (|p.1 - a| * |p.2 - b|) ^ ε) ?_ ?_
    ((integrable_const 1).add hμ) (ae_of_all _ fun p => (hK₁ p.1).mul (hK₂ p.2))
  · exact Eventually.of_forall fun q =>
      (((hm₁ q).comp measurable_fst).mul ((hm₂ q).comp measurable_snd)).aestronglyMeasurable
  · filter_upwards [Ioo_mem_nhdsGT hε] with q hq
    refine ae_of_all _ fun p => ?_
    obtain ⟨h1, h1'⟩ := hb₁ q hq p.1
    obtain ⟨h2, h2'⟩ := hb₂ q hq p.2
    rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg h1 h2)]
    calc K₁ q p.1 * K₂ q p.2 ≤ |p.1 - a| ^ q * |p.2 - b| ^ q :=
          mul_le_mul h1' h2' h2 (Real.rpow_nonneg (abs_nonneg _) _)
      _ = (|p.1 - a| * |p.2 - b|) ^ q := (Real.mul_rpow (abs_nonneg _) (abs_nonneg _)).symm
      _ ≤ 1 + (|p.1 - a| * |p.2 - b|) ^ ε :=
          rpow_le_one_add_rpow (by positivity) hq.1.le hq.2.le

theorem tendsto_quad_degree_zero_of_rpow [IsProbabilityMeasure μ] {a b ε : ℝ} (hε : 0 < ε)
    (hμ : Integrable (fun p : ℝ × ℝ => (|p.1 - a| * |p.2 - b|) ^ ε) μ) :
    Tendsto (fun q => quadLL q a b μ) (𝓝[>] 0) (𝓝 (μ.real {p | p.1 < a ∧ p.2 < b})) ∧
    Tendsto (fun q => quadUU q a b μ) (𝓝[>] 0) (𝓝 (μ.real {p | a < p.1 ∧ b < p.2})) ∧
    Tendsto (fun q => quadUL q a b μ) (𝓝[>] 0) (𝓝 (μ.real {p | a < p.1 ∧ p.2 < b})) ∧
    Tendsto (fun q => quadLU q a b μ) (𝓝[>] 0) (𝓝 (μ.real {p | p.1 < a ∧ b < p.2})) ∧
    Tendsto (fun q => quadTotal q a b μ) (𝓝[>] 0) (𝓝 (μ.real {p | p.1 ≠ a ∧ p.2 ≠ b})) := by
  have bl : ∀ c, ∀ q ∈ Ioo (0 : ℝ) ε, ∀ x, 0 ≤ rLowerDev q c x ∧
      rLowerDev q c x ≤ |x - c| ^ q := fun c q hq x =>
    ⟨rLowerDev_nonneg q c x, rLowerDev_le_abs_rpow q c x⟩
  have bu : ∀ c, ∀ q ∈ Ioo (0 : ℝ) ε, ∀ x, 0 ≤ rUpperDev q c x ∧
      rUpperDev q c x ≤ |x - c| ^ q := fun c q hq x =>
    ⟨rUpperDev_nonneg q c x, rUpperDev_le_abs_rpow q c x⟩
  have hLL : Tendsto (fun q => quadLL q a b μ) (𝓝[>] 0)
      (𝓝 (μ.real {p | p.1 < a ∧ p.2 < b})) := by
    rw [← integral_ite_mul_ite _ _ measurableSet_Iio measurableSet_Iio]
    exact tendsto_integral_kernel_mul_rpow hε (measurable_rLowerDev' · a)
      (measurable_rLowerDev' · b) (tendsto_rLowerDev_zero a) (tendsto_rLowerDev_zero b)
      (bl a) (bl b) hμ
  have hUU : Tendsto (fun q => quadUU q a b μ) (𝓝[>] 0)
      (𝓝 (μ.real {p | a < p.1 ∧ b < p.2})) := by
    rw [← integral_ite_mul_ite _ _ measurableSet_Ioi measurableSet_Ioi]
    exact tendsto_integral_kernel_mul_rpow hε (measurable_rUpperDev' · a)
      (measurable_rUpperDev' · b) (tendsto_rUpperDev_zero a) (tendsto_rUpperDev_zero b)
      (bu a) (bu b) hμ
  have hUL : Tendsto (fun q => quadUL q a b μ) (𝓝[>] 0)
      (𝓝 (μ.real {p | a < p.1 ∧ p.2 < b})) := by
    rw [← integral_ite_mul_ite _ _ measurableSet_Ioi measurableSet_Iio]
    exact tendsto_integral_kernel_mul_rpow hε (measurable_rUpperDev' · a)
      (measurable_rLowerDev' · b) (tendsto_rUpperDev_zero a) (tendsto_rLowerDev_zero b)
      (bu a) (bl b) hμ
  have hLU : Tendsto (fun q => quadLU q a b μ) (𝓝[>] 0)
      (𝓝 (μ.real {p | p.1 < a ∧ b < p.2})) := by
    rw [← integral_ite_mul_ite _ _ measurableSet_Iio measurableSet_Ioi]
    exact tendsto_integral_kernel_mul_rpow hε (measurable_rLowerDev' · a)
      (measurable_rUpperDev' · b) (tendsto_rLowerDev_zero a) (tendsto_rUpperDev_zero b)
      (bl a) (bu b) hμ
  refine ⟨hLL, hUU, hUL, hLU, ?_⟩
  rw [← measureReal_off_axis]
  exact ((hLL.add hUU).add hUL).add hLU

theorem tendsto_quad_weights_degree_zero_of_rpow [IsProbabilityMeasure μ] {a b ε : ℝ}
    (hε : 0 < ε) (hμ : Integrable (fun p : ℝ × ℝ => (|p.1 - a| * |p.2 - b|) ^ ε) μ)
    (hax : μ {p | p.1 = a} = 0) (hay : μ {p | p.2 = b} = 0) :
    Tendsto (fun q => quadLL q a b μ / quadTotal q a b μ) (𝓝[>] 0)
      (𝓝 (μ.real {p | p.1 ≤ a ∧ p.2 ≤ b})) ∧
    Tendsto (fun q => quadUU q a b μ / quadTotal q a b μ) (𝓝[>] 0)
      (𝓝 (μ.real {p | a < p.1 ∧ b < p.2})) ∧
    Tendsto (fun q => quadUL q a b μ / quadTotal q a b μ) (𝓝[>] 0)
      (𝓝 (μ.real {p | a < p.1 ∧ p.2 ≤ b})) ∧
    Tendsto (fun q => quadLU q a b μ / quadTotal q a b μ) (𝓝[>] 0)
      (𝓝 (μ.real {p | p.1 ≤ a ∧ b < p.2})) := by
  obtain ⟨hLL, hUU, hUL, hLU, hS⟩ := tendsto_quad_degree_zero_of_rpow hε hμ
  have hN : μ ({p : ℝ × ℝ | p.1 = a} ∪ {p | p.2 = b}) = 0 := measure_union_null hax hay
  have hoff : μ.real {p : ℝ × ℝ | p.1 ≠ a ∧ p.2 ≠ b} = 1 := by
    rw [measureReal_eq_of_subset_union_null (subset_univ _) (fun p _ => by
      by_cases h1 : p.1 = a
      · exact Or.inr (Or.inl h1)
      · by_cases h2 : p.2 = b
        · exact Or.inr (Or.inr h2)
        · exact Or.inl ⟨h1, h2⟩) hN]
    simp
  have key : ∀ {S T : Set (ℝ × ℝ)}, S ⊆ T → T ⊆ S ∪ ({p | p.1 = a} ∪ {p | p.2 = b}) →
      ∀ {f : ℝ → ℝ}, Tendsto f (𝓝[>] 0) (𝓝 (μ.real S)) →
      Tendsto (fun q => f q / quadTotal q a b μ) (𝓝[>] 0) (𝓝 (μ.real T)) := by
    intro S T hST hTS f hf
    rw [← measureReal_eq_of_subset_union_null hST hTS hN, ← div_one (μ.real S), ← hoff]
    exact hf.div hS (by rw [hoff]; exact one_ne_zero)
  refine ⟨key (S := {p | p.1 < a ∧ p.2 < b})
      (T := {p | p.1 ≤ a ∧ p.2 ≤ b}) (fun p hp => ⟨hp.1.le, hp.2.le⟩) ?_ hLL,
    key subset_rfl subset_union_left hUU,
    key (S := {p | a < p.1 ∧ p.2 < b})
      (T := {p | a < p.1 ∧ p.2 ≤ b}) (fun p hp => ⟨hp.1, hp.2.le⟩) ?_ hUL,
    key (S := {p | p.1 < a ∧ b < p.2})
      (T := {p | p.1 ≤ a ∧ b < p.2}) (fun p hp => ⟨hp.1.le, hp.2⟩) ?_ hLU⟩
  · intro p hp
    rcases hp.1.lt_or_eq with h1 | h1
    · rcases hp.2.lt_or_eq with h2 | h2
      · exact Or.inl ⟨h1, h2⟩
      · exact Or.inr (Or.inr h2)
    · exact Or.inr (Or.inl h1)
  · intro p hp
    rcases hp.2.lt_or_eq with h2 | h2
    · exact Or.inl ⟨hp.1, h2⟩
    · exact Or.inr (Or.inr h2)
  · intro p hp
    rcases hp.1.lt_or_eq with h1 | h1
    · exact Or.inl ⟨h1, hp.2⟩
    · exact Or.inr (Or.inl h1)

theorem off_axis_eq_one_sub [IsProbabilityMeasure μ] (a b : ℝ) :
    μ.real {p : ℝ × ℝ | p.1 ≠ a ∧ p.2 ≠ b} = 1 - μ.real {p | p.1 = a ∨ p.2 = b} := by
  have hm : MeasurableSet {p : ℝ × ℝ | p.1 = a ∨ p.2 = b} :=
    ((measurableSet_singleton a).preimage measurable_fst).union
      ((measurableSet_singleton b).preimage measurable_snd)
  have hc : {p : ℝ × ℝ | p.1 ≠ a ∧ p.2 ≠ b} = {p | p.1 = a ∨ p.2 = b}ᶜ := by
    ext p
    simp [not_or]
  rw [hc, measureReal_compl hm, probReal_univ]

end DegreeZeroLimit

section DensityRecovery

open MeasureTheory Filter Topology Set

lemma hasDerivAt_integral_max_mul {g : ℝ → ℝ} (hg : Continuous g) (hgc : HasCompactSupport g)
    (a : ℝ) :
    HasDerivAt (fun a' => ∫ x, max (a' - x) 0 * g x)
      (∫ x, (if x < a then (1 : ℝ) else 0) * g x) a := by
  have hcF : ∀ a' : ℝ, Continuous fun x => max (a' - x) 0 * g x := fun a' =>
    ((continuous_const.sub continuous_id).max continuous_const).mul hg
  refine (hasDerivAt_integral_of_dominated_loc_of_lip (bound := fun x => |g x|) (s := univ)
    univ_mem (Eventually.of_forall fun a' => (hcF a').aestronglyMeasurable)
    ((hcF a).integrable_of_hasCompactSupport hgc.mul_left)
    ((Measurable.ite measurableSet_Iio measurable_const measurable_const).mul
      hg.measurable).aestronglyMeasurable
    (ae_of_all _ fun x => ?_) (hg.integrable_of_hasCompactSupport hgc).abs ?_).2
  · refine (LipschitzWith.of_dist_le_mul fun a₁ a₂ => ?_).lipschitzOnWith
    rw [Real.dist_eq, Real.dist_eq, ← sub_mul, abs_mul, Real.coe_nnabs, abs_abs, mul_comm]
    exact mul_le_mul_of_nonneg_left
      ((abs_max_sub_max_le_abs _ _ _).trans_eq (by ring_nf)) (abs_nonneg _)
  · have hne : ∀ᵐ x : ℝ, x ≠ a := by simp [ae_iff, measure_singleton]
    filter_upwards [hne] with x hx
    rcases hx.lt_or_gt with h | h
    · simp only [Pi.mul_apply, h, ↓reduceIte, one_mul]
      have : HasDerivAt (fun a' => (a' - x) * g x) (1 * g x) a :=
        ((hasDerivAt_id a).sub_const x).mul_const (g x)
      rw [one_mul] at this
      refine this.congr_of_eventuallyEq ?_
      filter_upwards [lt_mem_nhds h] with a' ha'
      rw [max_eq_left (sub_nonneg.mpr ha'.le)]
    · simp only [Pi.mul_apply, not_lt.mpr h.le, ↓reduceIte, zero_mul]
      refine (hasDerivAt_const a (0 : ℝ)).congr_of_eventuallyEq ?_
      filter_upwards [gt_mem_nhds h] with a' ha'
      rw [max_eq_right (sub_nonpos.mpr ha'.le), zero_mul]

lemma hasDerivAt_integral_ite_mul {g : ℝ → ℝ} (hg : Continuous g) (hgi : Integrable g) (a : ℝ) :
    HasDerivAt (fun a' => ∫ x, (if x < a' then (1 : ℝ) else 0) * g x) (g a) a := by
  have e : ∀ a', ∫ x, (if x < a' then (1 : ℝ) else 0) * g x =
      (∫ x in Iic 0, g x) + ∫ x in (0)..a', g x := fun a' => by
    have : (fun x => (if x < a' then (1 : ℝ) else 0) * g x) = (Iio a').indicator g := by
      ext x; by_cases h : x < a' <;> simp [indicator, h]
    rw [this, integral_indicator measurableSet_Iio, ← integral_Iic_eq_integral_Iio,
      ← intervalIntegral.integral_Iic_sub_Iic hgi.integrableOn hgi.integrableOn, add_sub_cancel]
  simp_rw [e]
  exact (intervalIntegral.integral_hasDerivAt_right (hgi.intervalIntegrable)
    (hg.stronglyMeasurableAtFilter _ _) hg.continuousAt).const_add _

variable {f : ℝ × ℝ → ℝ}

lemma exists_square (hfc : HasCompactSupport f) :
    ∃ R, ∀ x y : ℝ, (R < |x| ∨ R < |y|) → f (x, y) = 0 := by
  obtain ⟨R, hR⟩ := hfc.isCompact.isBounded.subset_closedBall 0
  refine ⟨R, fun x y hxy => image_eq_zero_of_notMem_tsupport fun hmem => ?_⟩
  have := hR hmem
  rw [Metric.mem_closedBall, dist_zero_right, Prod.norm_def, Real.norm_eq_abs,
    Real.norm_eq_abs, max_le_iff] at this
  rcases hxy with h | h <;> linarith [this.1, this.2]

def partialHinge (f : ℝ × ℝ → ℝ) (b x : ℝ) : ℝ := ∫ y, max (b - y) 0 * f (x, y)

lemma continuous_partialHinge (hf : Continuous f) (hfc : HasCompactSupport f) (b : ℝ) :
    Continuous (partialHinge f b) := by
  obtain ⟨R, hR⟩ := exists_square hfc
  have e : partialHinge f b = fun x => ∫ y in Icc (-R) R, max (b - y) 0 * f (x, y) := by
    ext x
    refine (setIntegral_eq_integral_of_forall_compl_eq_zero fun y hy => ?_).symm
    rw [hR x y (Or.inr ?_), mul_zero]
    rw [mem_Icc, not_and_or, not_le, not_le] at hy
    exact lt_abs.mpr (hy.elim (fun h => Or.inr (by linarith)) Or.inl)
  rw [e]
  exact continuous_parametric_integral_of_continuous
    (((continuous_const.sub continuous_snd).max continuous_const).mul
      (hf.comp (continuous_fst.prodMk continuous_snd))) isCompact_Icc

lemma hasCompactSupport_partialHinge (hfc : HasCompactSupport f) (b : ℝ) :
    HasCompactSupport (partialHinge f b) := by
  obtain ⟨R, hR⟩ := exists_square hfc
  refine HasCompactSupport.intro (isCompact_Icc (a := -R) (b := R)) fun x hx => ?_
  have hx' : R < |x| := by
    rw [mem_Icc, not_and_or, not_le, not_le] at hx
    exact lt_abs.mpr (hx.elim (fun h => Or.inr (by linarith)) Or.inl)
  simp [partialHinge, hR x _ (Or.inl hx')]

lemma hasCompactSupport_section (hfc : HasCompactSupport f) (a : ℝ) :
    HasCompactSupport fun y => f (a, y) := by
  obtain ⟨R, hR⟩ := exists_square hfc
  refine HasCompactSupport.intro (isCompact_Icc (a := -R) (b := R)) fun y hy => hR a y (Or.inr ?_)
  rw [mem_Icc, not_and_or, not_le, not_le] at hy
  exact lt_abs.mpr (hy.elim (fun h => Or.inr (by linarith)) Or.inl)

lemma mCoLPM_withDensity (hf : Continuous f) (hfc : HasCompactSupport f) (hf0 : ∀ p, 0 ≤ f p)
    (a b : ℝ) :
    mCoLPM 1 1 a b (volume.withDensity fun p => ENNReal.ofReal (f p)) =
      ∫ x, max (a - x) 0 * partialHinge f b x := by
  have hint : Integrable (fun p : ℝ × ℝ => f p * (max (a - p.1) 0 * max (b - p.2) 0)) :=
    (hf.mul (((continuous_const.sub continuous_fst).max continuous_const).mul
      ((continuous_const.sub continuous_snd).max continuous_const))).integrable_of_hasCompactSupport
      hfc.mul_right
  rw [mCoLPM, integral_withDensity_eq_integral_toReal_smul hf.measurable.ennreal_ofReal
    (Eventually.of_forall fun _ => ENNReal.ofReal_lt_top)]
  simp only [ENNReal.toReal_ofReal (hf0 _), smul_eq_mul, lowerDev_one_eq_hinge]
  rw [Measure.volume_eq_prod, integral_prod _ hint]
  refine integral_congr_ae (Eventually.of_forall fun x => ?_)
  simp only [partialHinge, ← integral_const_mul]
  congr 1
  ext y
  ring

theorem hinge_fourth_deriv_eq_density (hf : Continuous f) (hfc : HasCompactSupport f)
    (hf0 : ∀ p, 0 ≤ f p) :
    ∃ H₁ H₂ H₃ : ℝ → ℝ → ℝ,
      (∀ a b, HasDerivAt (fun a' => mCoLPM 1 1 a' b
        (volume.withDensity fun p => ENNReal.ofReal (f p))) (H₁ a b) a) ∧
      (∀ a b, HasDerivAt (fun a' => H₁ a' b) (H₂ a b) a) ∧
      (∀ a b, HasDerivAt (fun b' => H₂ a b') (H₃ a b) b) ∧
      (∀ a b, HasDerivAt (fun b' => H₃ a b') (f (a, b)) b) := by
  refine ⟨fun a b => ∫ x, (if x < a then (1 : ℝ) else 0) * partialHinge f b x,
    fun a b => partialHinge f b a,
    fun a b => ∫ y, (if y < b then (1 : ℝ) else 0) * f (a, y), fun a b => ?_, fun a b => ?_,
    fun a b => ?_, fun a b => ?_⟩
  · simp_rw [mCoLPM_withDensity hf hfc hf0]
    exact hasDerivAt_integral_max_mul (continuous_partialHinge hf hfc b)
      (hasCompactSupport_partialHinge hfc b) a
  · exact hasDerivAt_integral_ite_mul (continuous_partialHinge hf hfc b)
      ((continuous_partialHinge hf hfc b).integrable_of_hasCompactSupport
        (hasCompactSupport_partialHinge hfc b)) a
  · exact hasDerivAt_integral_max_mul (hf.comp (continuous_const.prodMk continuous_id))
      (hasCompactSupport_section hfc a) b
  · have hc : Continuous fun y => f (a, y) := hf.comp (continuous_const.prodMk continuous_id)
    exact hasDerivAt_integral_ite_mul hc
      (hc.integrable_of_hasCompactSupport (hasCompactSupport_section hfc a)) b

theorem deriv_hinge_fourth_eq_density (hf : Continuous f) (hfc : HasCompactSupport f)
    (hf0 : ∀ p, 0 ≤ f p) (a b : ℝ) :
    deriv (fun b' => deriv (fun b'' => deriv (fun a' => deriv (fun a'' =>
      mCoLPM 1 1 a'' b'' (volume.withDensity fun p => ENNReal.ofReal (f p))) a') a) b') b
      = f (a, b) := by
  obtain ⟨H₁, H₂, H₃, h1, h2, h3, h4⟩ := hinge_fourth_deriv_eq_density hf hfc hf0
  have e1 : ∀ b'', (fun a' => deriv (fun a'' => mCoLPM 1 1 a'' b''
      (volume.withDensity fun p => ENNReal.ofReal (f p))) a') = fun a' => H₁ a' b'' :=
    fun b'' => funext fun a' => (h1 a' b'').deriv
  simp_rw [e1, (h2 _ _).deriv, (h3 _ _).deriv, (h4 a b).deriv]

lemma hasDerivAt_integral_max_mul_of_integrable {g : ℝ → ℝ} (hg : Continuous g)
    (hg0 : ∀ x, 0 ≤ g x) (hgi : ∀ c, Integrable fun x => max (c - x) 0 * g x) (a : ℝ) :
    HasDerivAt (fun a' => ∫ x, max (a' - x) 0 * g x)
      (∫ x, (if x < a then (1 : ℝ) else 0) * g x) a := by
  have hcF : ∀ a' : ℝ, Continuous fun x => max (a' - x) 0 * g x := fun a' =>
    ((continuous_const.sub continuous_id).max continuous_const).mul hg
  have hb : Integrable ((Iic (a + 1)).indicator g) :=
    (hgi (a + 2)).mono' (hg.measurable.indicator measurableSet_Iic).aestronglyMeasurable
      (ae_of_all _ fun x => by
        by_cases hx : x ≤ a + 1
        · rw [indicator_of_mem (mem_Iic.mpr hx), Real.norm_of_nonneg (hg0 x)]
          exact le_mul_of_one_le_left (hg0 x) (le_max_of_le_left (by linarith))
        · rw [indicator_of_notMem (fun h => hx (mem_Iic.mp h)), norm_zero]
          exact mul_nonneg (le_max_right _ _) (hg0 x))
  refine (hasDerivAt_integral_of_dominated_loc_of_lip (bound := (Iic (a + 1)).indicator g)
    (s := Metric.ball a 1) (Metric.ball_mem_nhds a one_pos)
    (Eventually.of_forall fun a' => (hcF a').aestronglyMeasurable) (hgi a)
    ((Measurable.ite measurableSet_Iio measurable_const measurable_const).mul
      hg.measurable).aestronglyMeasurable
    (ae_of_all _ fun x => ?_) hb ?_).2
  · refine LipschitzOnWith.of_dist_le_mul fun a₁ ha₁ a₂ ha₂ => ?_
    rw [Real.dist_eq, Real.dist_eq, ← sub_mul, abs_mul, abs_of_nonneg (hg0 x), Real.coe_nnabs]
    by_cases hx : x ≤ a + 1
    · rw [indicator_of_mem (mem_Iic.mpr hx), abs_of_nonneg (hg0 x), mul_comm]
      exact mul_le_mul_of_nonneg_left
        ((abs_max_sub_max_le_abs _ _ _).trans_eq (by ring_nf)) (hg0 x)
    · rw [Metric.mem_ball, Real.dist_eq, abs_lt] at ha₁ ha₂
      rw [not_le] at hx
      rw [max_eq_right (by linarith), max_eq_right (by linarith), sub_zero, abs_zero, zero_mul]
      exact mul_nonneg (abs_nonneg _) (abs_nonneg _)
  · have hne : ∀ᵐ x : ℝ, x ≠ a := by simp [ae_iff, measure_singleton]
    filter_upwards [hne] with x hx
    rcases hx.lt_or_gt with h | h
    · simp only [Pi.mul_apply, h, ↓reduceIte, one_mul]
      have : HasDerivAt (fun a' => (a' - x) * g x) (1 * g x) a :=
        ((hasDerivAt_id a).sub_const x).mul_const (g x)
      rw [one_mul] at this
      refine this.congr_of_eventuallyEq ?_
      filter_upwards [lt_mem_nhds h] with a' ha'
      rw [max_eq_left (sub_nonneg.mpr ha'.le)]
    · simp only [Pi.mul_apply, not_lt.mpr h.le, ↓reduceIte, zero_mul]
      refine (hasDerivAt_const a (0 : ℝ)).congr_of_eventuallyEq ?_
      filter_upwards [gt_mem_nhds h] with a' ha'
      rw [max_eq_right (sub_nonpos.mpr ha'.le), zero_mul]

lemma integrableOn_Iic_of_hinge {g : ℝ → ℝ} (hg : Continuous g) (hg0 : ∀ x, 0 ≤ g x)
    (hgi : ∀ c, Integrable fun x => max (c - x) 0 * g x) (c : ℝ) : IntegrableOn g (Iic c) := by
  rw [← integrable_indicator_iff measurableSet_Iic]
  refine (hgi (c + 1)).mono' (hg.measurable.indicator measurableSet_Iic).aestronglyMeasurable
    (ae_of_all _ fun x => ?_)
  by_cases hx : x ≤ c
  · rw [indicator_of_mem (mem_Iic.mpr hx), Real.norm_of_nonneg (hg0 x)]
    exact le_mul_of_one_le_left (hg0 x) (le_max_of_le_left (by linarith))
  · rw [indicator_of_notMem (fun h => hx (mem_Iic.mp h)), norm_zero]
    exact mul_nonneg (le_max_right _ _) (hg0 x)

lemma hasDerivAt_integral_ite_mul_of_integrableOn {g : ℝ → ℝ} (hg : Continuous g)
    (hgi : ∀ c, IntegrableOn g (Iic c)) (a : ℝ) :
    HasDerivAt (fun a' => ∫ x, (if x < a' then (1 : ℝ) else 0) * g x) (g a) a := by
  have e : ∀ a', ∫ x, (if x < a' then (1 : ℝ) else 0) * g x =
      (∫ x in Iic 0, g x) + ∫ x in (0)..a', g x := fun a' => by
    have : (fun x => (if x < a' then (1 : ℝ) else 0) * g x) = (Iio a').indicator g := by
      ext x; by_cases h : x < a' <;> simp [indicator, h]
    rw [this, integral_indicator measurableSet_Iio, ← integral_Iic_eq_integral_Iio,
      ← intervalIntegral.integral_Iic_sub_Iic (hgi 0) (hgi a'), add_sub_cancel]
  simp_rw [e]
  exact (intervalIntegral.integral_hasDerivAt_right (hg.intervalIntegrable _ _)
    (hg.stronglyMeasurableAtFilter _ _) hg.continuousAt).const_add _

lemma continuous_partialHinge_of_dominated (hf : Continuous f)
    (hdom : ∀ b x₀ : ℝ, ∃ φ : ℝ → ℝ, Integrable φ ∧
      ∀ᶠ x in 𝓝 x₀, ∀ y, max (b - y) 0 * f (x, y) ≤ φ y)
    (hf0 : ∀ p, 0 ≤ f p) (b : ℝ) : Continuous (partialHinge f b) := by
  refine continuous_iff_continuousAt.mpr fun x₀ => ?_
  obtain ⟨φ, hφ, hle⟩ := hdom b x₀
  have hc : Continuous fun p : ℝ × ℝ => max (b - p.2) 0 * f (p.1, p.2) :=
    ((continuous_const.sub continuous_snd).max continuous_const).mul
      (hf.comp (continuous_fst.prodMk continuous_snd))
  refine continuousAt_of_dominated (bound := φ)
    (Eventually.of_forall fun x =>
      (hc.comp (continuous_const.prodMk continuous_id)).aestronglyMeasurable)
    ?_ hφ (ae_of_all _ fun y => (hc.comp (continuous_id.prodMk continuous_const)).continuousAt)
  filter_upwards [hle] with x hx
  refine ae_of_all _ fun y => ?_
  rw [Real.norm_of_nonneg (mul_nonneg (le_max_right _ _) (hf0 _))]
  exact hx y

lemma partialHinge_nonneg (hf0 : ∀ p, 0 ≤ f p) (b x : ℝ) : 0 ≤ partialHinge f b x :=
  integral_nonneg fun _ => mul_nonneg (le_max_right _ _) (hf0 _)

lemma integrable_max_mul_partialHinge
    (hH : ∀ a b : ℝ, Integrable fun p : ℝ × ℝ => f p * (max (a - p.1) 0 * max (b - p.2) 0))
    (c b : ℝ) : Integrable fun x => max (c - x) 0 * partialHinge f b x := by
  rw [Measure.volume_eq_prod] at hH
  refine (hH c b).integral_prod_left.congr (ae_of_all _ fun x => ?_)
  simp only [partialHinge, ← integral_const_mul]
  congr 1
  ext y
  ring

lemma mCoLPM_withDensity_of_integrable (hf : Continuous f) (hf0 : ∀ p, 0 ≤ f p) {a b : ℝ}
    (hint : Integrable fun p : ℝ × ℝ => f p * (max (a - p.1) 0 * max (b - p.2) 0)) :
    mCoLPM 1 1 a b (volume.withDensity fun p => ENNReal.ofReal (f p)) =
      ∫ x, max (a - x) 0 * partialHinge f b x := by
  rw [mCoLPM, integral_withDensity_eq_integral_toReal_smul hf.measurable.ennreal_ofReal
    (Eventually.of_forall fun _ => ENNReal.ofReal_lt_top)]
  simp only [ENNReal.toReal_ofReal (hf0 _), smul_eq_mul, lowerDev_one_eq_hinge]
  rw [Measure.volume_eq_prod] at hint ⊢
  rw [integral_prod _ hint]
  refine integral_congr_ae (Eventually.of_forall fun x => ?_)
  simp only [partialHinge, ← integral_const_mul]
  congr 1
  ext y
  ring

theorem hinge_fourth_deriv_eq_density_of_dominated (hf : Continuous f) (hf0 : ∀ p, 0 ≤ f p)
    (hH : ∀ a b : ℝ, Integrable fun p : ℝ × ℝ => f p * (max (a - p.1) 0 * max (b - p.2) 0))
    (hdom : ∀ b x₀ : ℝ, ∃ φ : ℝ → ℝ, Integrable φ ∧
      ∀ᶠ x in 𝓝 x₀, ∀ y, max (b - y) 0 * f (x, y) ≤ φ y) :
    ∃ H₁ H₂ H₃ : ℝ → ℝ → ℝ,
      (∀ a b, HasDerivAt (fun a' => mCoLPM 1 1 a' b
        (volume.withDensity fun p => ENNReal.ofReal (f p))) (H₁ a b) a) ∧
      (∀ a b, HasDerivAt (fun a' => H₁ a' b) (H₂ a b) a) ∧
      (∀ a b, HasDerivAt (fun b' => H₂ a b') (H₃ a b) b) ∧
      (∀ a b, HasDerivAt (fun b' => H₃ a b') (f (a, b)) b) := by
  have hg := continuous_partialHinge_of_dominated hf hdom hf0
  have hg0 := partialHinge_nonneg hf0
  have hgi := integrable_max_mul_partialHinge hH
  have hs : ∀ a, Continuous fun y => f (a, y) := fun a =>
    hf.comp (continuous_const.prodMk continuous_id)
  have hsi : ∀ a c, Integrable fun y => max (c - y) 0 * f (a, y) := fun a c => by
    obtain ⟨φ, hφ, hle⟩ := hdom c a
    exact hφ.mono' (((continuous_const.sub continuous_id).max continuous_const).mul
      (hs a)).aestronglyMeasurable (ae_of_all _ fun y => by
        rw [Real.norm_of_nonneg (mul_nonneg (le_max_right _ _) (hf0 _))]
        exact hle.self_of_nhds y)
  refine ⟨fun a b => ∫ x, (if x < a then (1 : ℝ) else 0) * partialHinge f b x,
    fun a b => partialHinge f b a,
    fun a b => ∫ y, (if y < b then (1 : ℝ) else 0) * f (a, y), fun a b => ?_, fun a b => ?_,
    fun a b => ?_, fun a b => ?_⟩
  · simp_rw [mCoLPM_withDensity_of_integrable hf hf0 (hH _ b)]
    exact hasDerivAt_integral_max_mul_of_integrable (hg b) (hg0 b) (hgi · b) a
  · exact hasDerivAt_integral_ite_mul_of_integrableOn (hg b)
      (integrableOn_Iic_of_hinge (hg b) (hg0 b) (hgi · b)) a
  · exact hasDerivAt_integral_max_mul_of_integrable (hs a) (fun _ => hf0 _) (hsi a) b
  · exact hasDerivAt_integral_ite_mul_of_integrableOn (hs a)
      (integrableOn_Iic_of_hinge (hs a) (fun _ => hf0 _) (hsi a)) b

theorem deriv_hinge_fourth_eq_density_of_dominated (hf : Continuous f) (hf0 : ∀ p, 0 ≤ f p)
    (hH : ∀ a b : ℝ, Integrable fun p : ℝ × ℝ => f p * (max (a - p.1) 0 * max (b - p.2) 0))
    (hdom : ∀ b x₀ : ℝ, ∃ φ : ℝ → ℝ, Integrable φ ∧
      ∀ᶠ x in 𝓝 x₀, ∀ y, max (b - y) 0 * f (x, y) ≤ φ y) (a b : ℝ) :
    deriv (fun b' => deriv (fun b'' => deriv (fun a' => deriv (fun a'' =>
      mCoLPM 1 1 a'' b'' (volume.withDensity fun p => ENNReal.ofReal (f p))) a') a) b') b
      = f (a, b) := by
  obtain ⟨H₁, H₂, H₃, h1, h2, h3, h4⟩ := hinge_fourth_deriv_eq_density_of_dominated hf hf0 hH hdom
  have e1 : ∀ b'', (fun a' => deriv (fun a'' => mCoLPM 1 1 a'' b''
      (volume.withDensity fun p => ENNReal.ofReal (f p))) a') = fun a' => H₁ a' b'' :=
    fun b'' => funext fun a' => (h1 a' b'').deriv
  simp_rw [e1, (h2 _ _).deriv, (h3 _ _).deriv, (h4 a b).deriv]

end DensityRecovery

section Rectangle

open Finset

variable {n : ℕ}

lemma clpmND_two_one (X : Matrix (Fin n) (Fin 2) ℝ) (t : Fin 2 → ℝ) :
    clpmND X t 1 = (∑ i, max (t 0 - X i 0) 0 * max (t 1 - X i 1) 0) / n := by
  simp only [clpmND, lowRow, Fin.prod_univ_two, rLowerDev_eq_max one_ne_zero, Real.rpow_one]

lemma totalND_two_one (X : Matrix (Fin n) (Fin 2) ℝ) (t : Fin 2 → ℝ) :
    totalND X t 1 = (∑ i, |X i 0 - t 0| * |X i 1 - t 1|) / n := by
  rw [totalND_eq_absND (by norm_num)]
  simp only [absND, absRow, Fin.prod_univ_two, Real.rpow_one]

theorem coLPMnD_R_two_one_norm (X : Matrix (Fin n) (Fin 2) ℝ) (t : Fin 2 → ℝ)
    (htot : 0 < ∑ i, |X i 0 - t 0| * |X i 1 - t 1|) :
    coLPMnD_R X t 1 true = (∑ i, max (t 0 - X i 0) 0 * max (t 1 - X i 1) 0)
      / ∑ i, |X i 0 - t 0| * |X i 1 - t 1| := by
  have hn : (0 : ℝ) < n := by
    rcases Nat.eq_zero_or_pos n with rfl | h
    · simp at htot
    · exact_mod_cast h
  have ht : 0 < totalND X t 1 := by rw [totalND_two_one]; exact div_pos htot hn
  rw [coLPMnD_R_eq zero_le_one, ite_eq_left ⟨one_ne_zero, rfl⟩, normND, ite_eq_left ht]
  change clpmND X t 1 / totalND X t 1 = _
  rw [clpmND_two_one, totalND_two_one]
  field_simp

def rectX : Matrix (Fin 3) (Fin 2) ℝ := !![0, 0; 1, 1; 4, 4]

lemma coLPMnD_R_rectX (s t : ℝ) (h : 0 < |0 - s| * |0 - t| + |1 - s| * |1 - t|
    + |4 - s| * |4 - t|) :
    coLPMnD_R rectX ![s, t] 1 true = (max (s - 0) 0 * max (t - 0) 0 + max (s - 1) 0 * max (t - 1) 0
      + max (s - 4) 0 * max (t - 4) 0)
      / (|0 - s| * |0 - t| + |1 - s| * |1 - t| + |4 - s| * |4 - t|) := by
  rw [coLPMnD_R_two_one_norm]
  · simp [rectX, Fin.sum_univ_three]
  · simpa [rectX, Fin.sum_univ_three] using h

theorem coLPMnD_R_rectX_values :
    coLPMnD_R rectX ![3, 3] 1 true = 13 / 14 ∧ coLPMnD_R rectX ![2, 3] 1 true = 4 / 5 ∧
      coLPMnD_R rectX ![3, 2] 1 true = 4 / 5 ∧ coLPMnD_R rectX ![2, 2] 1 true = 5 / 9 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;>
  · rw [coLPMnD_R_rectX _ _ (by norm_num [abs_of_nonneg, abs_of_nonpos])]
    norm_num [abs_of_nonneg, abs_of_nonpos]

theorem coLPMnD_R_rectangle_eq :
    coLPMnD_R rectX ![3, 3] 1 true - coLPMnD_R rectX ![2, 3] 1 true
      - coLPMnD_R rectX ![3, 2] 1 true + coLPMnD_R rectX ![2, 2] 1 true = -73 / 630 := by
  obtain ⟨h1, h2, h3, h4⟩ := coLPMnD_R_rectX_values
  rw [h1, h2, h3, h4]
  norm_num

lemma coLPMnD_R_rectX_v8898017 (s t : ℝ) (h : 0 < |0 - s| * |0 - t| + |1 - s| * |1 - t|
    + |4 - s| * |4 - t|) :
    coLPMnD_R_v8898017 rectX ![s, t] 1 true = coLPMnD_R rectX ![s, t] 1 true := by
  refine (coLPMnD_R_eq_v8898017_of zero_le_one _ _ _ (Or.inr (Or.inr ?_))).symm
  rw [totalND_two_one]
  exact div_pos (by simpa [rectX, Fin.sum_univ_three] using h) (by norm_num)

theorem coLPMnD_R_rectangle_eq_v8898017 :
    coLPMnD_R_v8898017 rectX ![3, 3] 1 true - coLPMnD_R_v8898017 rectX ![2, 3] 1 true
      - coLPMnD_R_v8898017 rectX ![3, 2] 1 true + coLPMnD_R_v8898017 rectX ![2, 2] 1 true
      = -73 / 630 := by
  rw [coLPMnD_R_rectX_v8898017 3 3 (by norm_num), coLPMnD_R_rectX_v8898017 2 3 (by norm_num),
    coLPMnD_R_rectX_v8898017 3 2 (by norm_num), coLPMnD_R_rectX_v8898017 2 2 (by norm_num)]
  exact coLPMnD_R_rectangle_eq

theorem coUPMnD_R_neg_v8898017 {d : ℝ} (hd : 0 < d) (X : Matrix (Fin n) (Fin 2) ℝ)
    (t : Fin 2 → ℝ) (norm : Bool) :
    coUPMnD_R_v8898017 (-X) (-t) d norm = coLPMnD_R_v8898017 X t d norm := by
  rw [coUPMnD_R_eq_v8898017 hd.le, coLPMnD_R_eq_v8898017 hd.le, totalND_neg hd.ne',
    cupmND_neg hd.ne']

theorem exists_negative_rectangle_batch_v8898017 :
    ∃ (X : Matrix (Fin 3) (Fin 2) ℝ) (a b : Fin 2 → ℝ), (∀ j, a j < b j) ∧
      coLPMnDBatch_v8898017 X (Matrix.of fun _ : Fin 1 => b) 1 true 0
        - coLPMnDBatch_v8898017 X (Matrix.of fun _ : Fin 1 => ![a 0, b 1]) 1 true 0
        - coLPMnDBatch_v8898017 X (Matrix.of fun _ : Fin 1 => ![b 0, a 1]) 1 true 0
        + coLPMnDBatch_v8898017 X (Matrix.of fun _ : Fin 1 => a) 1 true 0 < 0 := by
  refine ⟨rectX, ![2, 2], ![3, 3], fun j => by fin_cases j <;> norm_num, ?_⟩
  have hrow : ∀ t : Fin 2 → ℝ, (Matrix.of fun _ : Fin 1 => t) 0 = t := fun _ => rfl
  simp only [coLPMnDBatch_eq_v8898017 zero_le_one, hrow, Matrix.cons_val_zero,
    Matrix.cons_val_one]
  rw [coLPMnD_R_rectangle_eq_v8898017]
  norm_num

theorem exists_negative_survival_rectangle_v8898017 :
    ∃ (X : Matrix (Fin 3) (Fin 2) ℝ) (a b : Fin 2 → ℝ), (∀ j, a j < b j) ∧
      coLPMnDBatch_v8898017 (-X) (-(Matrix.of fun _ : Fin 1 => a)) 1 true 0
        - coLPMnDBatch_v8898017 (-X) (-(Matrix.of fun _ : Fin 1 => ![b 0, a 1])) 1 true 0
        - coLPMnDBatch_v8898017 (-X) (-(Matrix.of fun _ : Fin 1 => ![a 0, b 1])) 1 true 0
        + coLPMnDBatch_v8898017 (-X) (-(Matrix.of fun _ : Fin 1 => b)) 1 true 0 < 0 := by
  refine ⟨-rectX, ![-3, -3], ![-2, -2], fun j => by fin_cases j <;> norm_num, ?_⟩
  have hrow : ∀ t : Fin 2 → ℝ, (Matrix.of fun _ : Fin 1 => t) 0 = t := fun _ => rfl
  simp only [coLPMnDBatch_neg_eq_coUPMnD_R_v8898017 one_pos, hrow, Matrix.cons_val_zero,
    Matrix.cons_val_one]
  have e : ∀ s t : ℝ, coUPMnD_R_v8898017 (-rectX) ![-s, -t] 1 true =
      coLPMnD_R_v8898017 rectX ![s, t] 1 true := fun s t => by
    rw [← coUPMnD_R_neg_v8898017 one_pos rectX ![s, t] true]
    congr 1
    ext j
    fin_cases j <;> rfl
  have h := coLPMnD_R_rectangle_eq_v8898017
  rw [← e 3 3, ← e 2 3, ← e 3 2, ← e 2 2] at h
  norm_num at h ⊢
  linarith

end Rectangle

end

end NNS.GapCoMoments
