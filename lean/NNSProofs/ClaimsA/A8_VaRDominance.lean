-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.RealDegree
public import NNSProofs.GapDominance

namespace NNS.ClaimsA

open Finset Set MeasureTheory ProbabilityTheory GapDominance
@[expose] public noncomputable section

namespace A8_VaRDominance

variable {n m : ℕ} {x : Fin n → ℝ} {y : Fin m → ℝ} {p d : ℝ}

def balance (d p : ℝ) (x : Fin n → ℝ) (t : ℝ) : ℝ :=
  (1 - p) * rLPM d t x - p * rUPM d t x

def lower (d p : ℝ) (x : Fin n → ℝ) : ℝ :=
  if d = 0 then LPMVaR0 p x else
    sInf {t | t ∈ Icc (varMin x) (varMax x) ∧ balance d p x t = 0}

def upper (d p : ℝ) (x : Fin n → ℝ) : ℝ := lower d (1 - p) x

private theorem min_lower (hd : 0 < d) (x : Fin n → ℝ) : rLPM d (varMin x) x = 0 := by
  simp [rLPM_eq_max hd.ne', max_eq_right (sub_nonpos.mpr (varMin_le x _)), Real.zero_rpow hd.ne']

private theorem max_upper (hd : 0 < d) (x : Fin n → ℝ) : rUPM d (varMax x) x = 0 := by
  simp [rUPM_eq_max hd.ne', max_eq_right (sub_nonpos.mpr (le_varMax x _)), Real.zero_rpow hd.ne']

private theorem lower_strict (hn : 0 < n) (hd : 0 < d) (x : Fin n → ℝ) {s t : ℝ}
    (hs : varMin x ≤ s) (hst : s < t) : rLPM d s x < rLPM d t x := by
  obtain ⟨j, hj⟩ := exists_eq_varMin hn x
  simp only [rLPM_eq_max hd.ne']
  apply div_lt_div_of_pos_right _ (Nat.cast_pos.mpr hn)
  apply sum_lt_sum
  · intro i _
    exact Real.rpow_le_rpow (le_max_right _ _) (max_le_max_right _ (by linarith)) hd.le
  · refine ⟨j, mem_univ _, ?_⟩
    rw [hj, max_eq_left (sub_nonneg.mpr hs), max_eq_left (sub_nonneg.mpr (hs.trans hst.le))]
    exact Real.rpow_lt_rpow (sub_nonneg.mpr hs) (by linarith) hd

private theorem balance_strict (hn : 0 < n) (hd : 0 < d) (hp : p ∈ Set.Icc (0 : ℝ) 1) :
    StrictMonoOn (balance d p x) (Icc (varMin x) (varMax x)) := by
  intro s hs t ht hst
  have hL := lower_strict hn hd x hs.1 hst
  have hU : rUPM d t x < rUPM d s x := by
    rw [rUPM_eq_rLPM_neg hd.ne', rUPM_eq_rLPM_neg hd.ne']
    exact lower_strict hn hd _ (by rw [varMin_neg hn]; exact neg_le_neg ht.2) (neg_lt_neg hst)
  unfold balance
  by_cases hp1 : p = 1
  · simpa [hp1] using hU
  · have h1 := (mul_lt_mul_of_pos_left hL (sub_pos.mpr (lt_of_le_of_ne hp.2 hp1)))
    have h2 := mul_le_mul_of_nonneg_left hU.le hp.1
    linarith

theorem root_spec (hn : 0 < n) (hd : 0 < d) (hp : p ∈ Set.Icc (0 : ℝ) 1) :
    lower d p x ∈ Icc (varMin x) (varMax x) ∧ balance d p x (lower d p x) = 0 ∧
      ∀ t ∈ Icc (varMin x) (varMax x), balance d p x t = 0 → t = lower d p x := by
  have hmin := min_lower hd x
  have hmax := max_upper hd x
  have hle := (varMin_le x ⟨0, hn⟩).trans (le_varMax x _)
  have hc : Continuous (balance d p x) :=
    (continuous_const.mul (rLPM_continuous hd x)).sub (continuous_const.mul (rUPM_continuous hd x))
  have hzero : (0 : ℝ) ∈ Set.Icc (balance d p x (varMin x)) (balance d p x (varMax x)) := by
    constructor
    · simp only [balance, hmin, mul_zero, zero_sub, neg_nonpos]
      exact mul_nonneg hp.1 (rUPM_nonneg ..)
    · simp only [balance, hmax, mul_zero, sub_zero]
      exact mul_nonneg (sub_nonneg.mpr hp.2) (rLPM_nonneg ..)
  obtain ⟨q, hq, hz⟩ := intermediate_value_Icc hle hc.continuousOn hzero
  have unique : ∀ t ∈ Icc (varMin x) (varMax x), balance d p x t = 0 → t = q :=
    fun _ ht he => (balance_strict hn hd hp).injOn ht hq (he.trans hz.symm)
  have he : {t | t ∈ Icc (varMin x) (varMax x) ∧ balance d p x t = 0} = {q} := by
    ext t
    exact ⟨fun h => unique t h.1 h.2, fun h => (mem_singleton_iff.mp h) ▸ ⟨hq, hz⟩⟩
  simpa only [lower, ite_eq_right hd.ne', he, csInf_singleton] using
    And.intro hq (And.intro hz unique)

theorem positive_total (_hd : 0 < d) (hx : ∃ i j, x i ≠ x j) (t : ℝ) :
    0 < rLPM d t x + rUPM d t x := by
  obtain ⟨i,j,hij⟩ := hx
  obtain ⟨k,hk⟩ : ∃ k, x k ≠ t := by
    by_contra! hall
    exact hij ((hall i).trans (hall j).symm)
  rw [rLPM_add_rUPM]
  exact div_pos (sum_pos' (fun _ _ => Real.rpow_nonneg (abs_nonneg _) _)
    ⟨k, mem_univ _, Real.rpow_pos_of_pos (abs_pos.mpr (sub_ne_zero.mpr hk)) _⟩)
    (Nat.cast_pos.mpr (Fin.pos i))

theorem lower_share (hd : 0 < d) (hp : p ∈ Set.Icc (0 : ℝ) 1) (hx : ∃ i j, x i ≠ x j) :
    let t := lower d p x
    rLPM d t x / (rLPM d t x + rUPM d t x) = p := by
  dsimp only
  obtain ⟨i,_,_⟩ := id hx
  apply (div_eq_iff (positive_total hd hx _).ne').mpr
  have h := (root_spec (x := x) (Fin.pos i) hd hp).2.1
  unfold balance at h
  linarith

theorem upper_share (hd : 0 < d) (hp : p ∈ Set.Icc (0 : ℝ) 1) (hx : ∃ i j, x i ≠ x j) :
    let t := upper d p x
    rUPM d t x / (rLPM d t x + rUPM d t x) = p := by
  dsimp only [upper]
  have h := lower_share hd (p := 1 - p) ⟨sub_nonneg.mpr hp.2, sub_le_self _ hp.1⟩ hx
  dsimp only at h
  rw [div_eq_iff (positive_total hd hx _).ne'] at h ⊢
  linarith

theorem zero_pair (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    lower 0 p x = quantile7 p x ∧ upper 0 p x = quantile7 (1 - p) x := by
  simp [upper, lower, LPMVaR0, clampP_of_mem hp,
    clampP_of_mem (p := 1-p) ⟨sub_nonneg.mpr hp.2, sub_le_self _ hp.1⟩]

theorem expectile (hn : 0 < n) (hp : p ∈ Set.Icc (0 : ℝ) 1) (t : ℝ) :
    t = lower 1 p x ↔ t ∈ Icc (varMin x) (varMax x) ∧
      (1-p) * LPM 1 t x = p * UPM 1 t x := by
  have h := root_spec (x := x) hn (d := 1) zero_lt_one hp
  have hL (a : ℝ) : rLPM 1 a x = LPM 1 a x := by simpa using rLPM_natCast 1 a x
  have hU (a : ℝ) : rUPM 1 a x = UPM 1 a x := by simpa using rUPM_natCast 1 a x
  simp only [balance, hL, hU, sub_eq_zero] at h
  exact ⟨fun e => e ▸ ⟨h.1,h.2.1⟩, fun ⟨ht, he⟩ => h.2.2 t ht he⟩

theorem ratio_at_expectile (hx : ∃ i j, x i ≠ x j) (hp : p ∈ Set.Icc (0 : ℝ) 1) :
    LPMratio 1 (lower 1 p x) x = p := by
  obtain ⟨i, j, hij⟩ := id hx
  have h := lower_share zero_lt_one hp hx (d := 1)
  dsimp only at h
  have hpos : 0 < LPM 1 (lower 1 p x) x + UPM 1 (lower 1 p x) x := by
    apply ratio_denom_pos
    by_contra! he
    exact hij ((he i).trans (he j).symm)
  rw [LPMratio_one_of_pos hpos]
  simpa only [← Nat.cast_one (R := ℝ), rLPM_natCast, rUPM_natCast] using h

theorem expectile_le (hn : 0 < n) (hp : p ∈ Set.Ioo (0 : ℝ) 1) (t : ℝ) :
    lower 1 p x ≤ t ↔ p ≤ LPMratio 1 t x := by
  have hq := (root_spec (x := x) hn zero_lt_one (d := 1) ⟨hp.1.le, hp.2.le⟩).1
  by_cases h0 : t < varMin x
  · rw [LPMratio_eq_zero_of_forall_lt 1 (fun i => h0.trans_le (varMin_le x i))]
    exact iff_of_false (not_le.mpr (h0.trans_le hq.1)) (not_le.mpr hp.1)
  by_cases h1 : varMax x ≤ t
  · rw [LPMratio_eq_one_of_forall_le hn 1 (fun i => (le_varMax x i).trans h1)]
    exact iff_of_true (hq.2.trans h1) hp.2.le
  have ht : t ∈ Icc (varMin x) (varMax x) := ⟨le_of_not_gt h0, (lt_of_not_ge h1).le⟩
  have hx := (varMin_lt_varMax_iff hn x).mp (ht.1.trans_lt (lt_of_not_ge h1))
  have h := ((LPMratio_strictMonoOn_Icc one_ne_zero hn x).le_iff_le hq ht).symm
  rwa [ratio_at_expectile hx ⟨hp.1.le, hp.2.le⟩] at h

theorem expectile_order (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    (∀ t, LPMratio 1 t x ≤ LPMratio 1 t y) ↔
      ∀ p ∈ Set.Ioo (0 : ℝ) 1, lower 1 p y ≤ lower 1 p x := by
  constructor
  · intro h p hp
    exact (expectile_le hm hp _).mpr (((expectile_le hn hp _).mp le_rfl).trans (h _))
  · intro h t
    by_contra! hgt
    let p := (LPMratio 1 t x + LPMratio 1 t y) / 2
    have hy := (LPMratio_mem_Icc 1 t y).1
    have hx := (LPMratio_mem_Icc 1 t x).2
    have hp : p ∈ Set.Ioo (0 : ℝ) 1 := ⟨by dsimp [p]; linarith, by dsimp [p]; linarith⟩
    have hq := (expectile_le (x := x) hn hp t).mpr (by dsimp [p]; linarith)
    have := (expectile_le hm hp t).mp ((h p hp).trans hq)
    dsimp [p] at this
    linarith

theorem strict_expectile_order (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    sdDomPairCentred x y 1 false = 1 ↔
      (∀ p ∈ Set.Ioo (0 : ℝ) 1, lower 1 p y ≤ lower 1 p x) ∧
      ∃ p ∈ Set.Ioo (0 : ℝ) 1, lower 1 p y < lower 1 p x := by
  rw [sdDomPairCentred_omega_eq hn hm, sdDomPair_omega_iff hn hm]
  constructor
  · rintro ⟨h, t, ht⟩
    refine ⟨(expectile_order hn hm x y).mp h, ?_⟩
    by_contra! rev
    exact ht (le_antisymm (h t) ((expectile_order hm hn y x).mpr rev t))
  · rintro ⟨h, p, hp, hlt⟩
    refine ⟨(expectile_order hn hm x y).mpr h, ?_⟩
    by_contra! rev
    exact hlt.not_ge ((expectile_order hm hn y x).mp (fun t => (rev t).ge) p hp)

theorem omega_identity (x : Fin n → ℝ) (t : ℝ) (h : 0 < LPM 1 t x) :
    LPMratio 1 t x = 1 / (1 + UPM 1 t x / LPM 1 t x) := by
  have hu := UPM_nonneg 1 t x
  rw [LPMratio_one_of_pos (by positivity)]
  field_simp

theorem tied_median :
    LPMVaR0 (1/2) ![(0 : ℝ), 0, 1] = 0 ∧ LPMratio 0 0 ![(0 : ℝ), 0, 1] = 2/3 := by
  have hs : Monotone ![(0 : ℝ),0,1] := by
    apply Fin.monotone_iff_le_succ.mpr
    intro i
    fin_cases i <;> norm_num
  constructor
  · rw [LPMVaR0, clampP_of_mem ⟨by norm_num, by norm_num⟩, quantile7_eq ⟨by norm_num, by norm_num⟩]
    norm_num [sortedAt, Tuple.sort_eq_refl_iff_monotone.mpr hs]
  · simp only [LPMratio, ↓reduceIte, LPM, lowerDev, Fin.sum_univ_three]
    norm_num

end A8_VaRDominance

open A8_VaRDominance
variable {n m : ℕ} {x : Fin n → ℝ} {y : Fin m → ℝ} {p d : ℝ}

theorem claim_C275 (hn : 0 < n) (x : Fin n → ℝ) :
    (∀ p ∈ Set.Ioo (0 : ℝ) 1, ∀ t, GapQuantile.quantile (empMeasure x) p ≤ t ↔
      p ≤ LPMratio 0 t x) ∧
    ((∃ i j, x i ≠ x j) → ∀ p ∈ Set.Icc (0 : ℝ) 1,
      (∀ t, t = lower 1 p x ↔ t ∈ Icc (varMin x) (varMax x) ∧ LPMratio 1 t x = p) ∧
      (∀ t, t = lower 1 p x ↔ t ∈ Icc (varMin x) (varMax x) ∧
        (1-p) * LPM 1 t x = p * UPM 1 t x)) ∧
    cdf (empMeasure ![(0 : ℝ),0,1]) (GapQuantile.quantile (empMeasure ![(0 : ℝ),0,1]) (1/2)) =
      2/3 ∧ (1/2 : ℝ) < 2/3 := by
  refine ⟨?_, ?_, GapQuantile.cdf_quantile_tied_sample, by norm_num⟩
  · intro p hp t
    simpa only [cdf_empMeasure hn, LPMratio, ↓reduceIte] using
      GapQuantile.quantile_le_iff (empMeasure x) hp.1 hp.2 (t := t)
  · intro hx p hp
    have hq := (root_spec (x := x) hn zero_lt_one (d := 1) hp).1
    have hr := ratio_at_expectile hx hp
    refine ⟨fun t => ⟨fun e => e ▸ ⟨hq, hr⟩, fun ⟨ht, he⟩ => ?_⟩, expectile hn hp⟩
    exact eq_of_LPMratio_eq one_ne_zero hn ht hq (he.trans hr.symm)

theorem claim_C336 (hn : 0 < n) (hm : 0 < m) {a b : ℝ}
    (ha : ∀ i, a ≤ x i) (hx : ∀ i, x i ≤ b) (hy : ∀ j, y j ≤ b) :
    ((∀ t, LPM 2 t x ≤ LPM 2 t y) ↔
      (∀ t ∈ Icc a b, LPM 2 t x ≤ LPM 2 t y) ∧ mean y ≤ mean x) ∧
    (sdDomPair x y 3 true = 1 → mean y ≤ mean x) :=
  ⟨⟨fun h => ⟨fun t _ => h t, mean_le_of_LPM_two_le hn hm h⟩,
    fun ⟨h,g⟩ => LPM_two_le_of_le_on_Icc hm ha hx hy g h⟩,
    fun h => mean_le_of_LPM_two_le hn hm ((sdDomPair_tsd_iff hn hm ..).mp h).1⟩

theorem claim_C338 (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ((∀ t, IsKnot x y t → LPM 0 t x ≤ LPM 0 t y) ↔ ∀ t, LPM 0 t x ≤ LPM 0 t y) ∧
    IsKnot ![(0 : ℝ)] ![1] 0 ∧ LPM 0 0 ![(1 : ℝ)] < LPM 0 0 ![(0 : ℝ)] ∧
    nnsFSD ![(0 : ℝ)] ![1] true = .yDom ∧ nnsFSD ![(0 : ℝ)] ![0] true = .noDom := by
  refine ⟨⟨LPM_zero_le_of_le_at_knots, fun h t _ => h t⟩, Or.inl ⟨0, rfl⟩,
    by norm_num [LPM, lowerDev, Fin.sum_univ_one], ?_, ?_⟩
  · apply (nnsFSD_eq_yDom_iff (by decide) (by decide) _ _ _).mpr
    apply (fsdUni_curveDom (by decide) (by decide) _ _ (Or.inl rfl) true).mp
    exact fsdUni_of_forall_lt (by decide) (by intro i; fin_cases i; norm_num)
  · simp [nnsFSD, fsdUni_self (by decide : 0 < 1)]

theorem claim_C339 (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ((∀ t, LPM 1 t x ≤ LPM 1 t y) ↔ ∀ t, IsKnot x y t → LPM 1 t x ≤ LPM 1 t y) ∧
    ((∀ t, LPM 2 t x ≤ LPM 2 t y) ↔ mean y ≤ mean x ∧
      (∀ t, IsKnot x y t → LPM 2 t x ≤ LPM 2 t y) ∧
      ∀ u w, IsKnot x y u → IsKnot x y w → (∀ i, x i ∉ Set.Ioo u w) →
        (∀ j, y j ∉ Set.Ioo u w) → tsdCoefA u x y < 0 →
        tsdCoefB u x y / tsdCoefA u x y ∈ Set.Ioo u w →
        let v := tsdCoefB u x y / tsdCoefA u x y
        LPM 2 v x ≤ LPM 2 v y) ∧
    (sdDomPair x y 2 true = 1 ↔ LPMDominates 1 x y) ∧
    (sdDomPair x y 3 true = 1 ↔ LPMDominates 2 x y) ∧
    tsdKnotOnly tsdVertexCexX tsdVertexCexY = 1 ∧
    sdDomPair tsdVertexCexX tsdVertexCexY 3 true = 0 :=
  ⟨⟨fun h t _ => h t, LPM_one_le_of_le_at_knots hm⟩,
    ⟨fun h => ⟨mean_le_of_LPM_two_le hn hm h, fun t _ => h t, fun _ _ _ _ _ _ _ _ => h _⟩,
      fun ⟨g,h,hv⟩ => LPM_two_le_of_le_at_knots hm g h hv⟩,
    sdDomPair_ssd_iff hn hm x y true, sdDomPair_tsd_iff hn hm x y true,
    tsdKnotOnly_vertexCex, sdDomPair_tsd_vertexCex true⟩

theorem claim_C351 (μ ν : Measure ℝ) [IsProbabilityMeasure μ] [IsProbabilityMeasure ν] :
    ((∀ t, mLPM 0 t μ ≤ mLPM 0 t ν) ↔ ∀ t, cdf μ t ≤ cdf ν t) ∧
    ((∀ t, Integrable (lowerDev 1 t) μ) → (∀ t, Integrable (lowerDev 1 t) ν) →
      ((∀ t, mLPM 1 t μ ≤ mLPM 1 t ν) ↔
       ∀ t, ∫ s in Set.Iic t, cdf μ s ≤ ∫ s in Set.Iic t, cdf ν s)) ∧
    ((∀ t, Integrable (lowerDev 2 t) μ) → (∀ t, Integrable (lowerDev 2 t) ν) →
      ((∀ t, mLPM 2 t μ ≤ mLPM 2 t ν) ↔ ∀ t, ∫ s in Set.Iic t, ∫ u in Set.Iic s, cdf μ u ≤
        ∫ s in Set.Iic t, ∫ u in Set.Iic s, cdf ν u)) :=
  ⟨by simp only [mLPM_zero_eq_cdf], forall_mLPM_one_le_iff, forall_mLPM_two_le_iff⟩

theorem claim_C430 (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    (lower 0 p x = quantile7 p x ∧ upper 0 p x = quantile7 (1-p) x) ∧
    ∀ d : ℝ, 0 < d → (∃ i j, x i ≠ x j) →
      (let t := lower d p x; rLPM d t x / (rLPM d t x + rUPM d t x) = p) ∧
      (let t := upper d p x; rUPM d t x / (rLPM d t x + rUPM d t x) = p) :=
  ⟨zero_pair hp x, fun _ hd hx => ⟨lower_share hd hp hx, upper_share hd hp hx⟩⟩

theorem claim_C440 (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ)
    (k : ℕ) (hk : k = 1 ∨ k = 2 ∨ k = 3) :
    sdDomPairCentred x y k true = 1 ↔
      (∀ t, LPM (k-1) t x ≤ LPM (k-1) t y) ∧ ∃ t, LPM (k-1) t x < LPM (k-1) t y := by
  have h : sdDomPairCentred x y k true = 1 ↔ LPMDominates (k-1) x y := by
    rcases hk with rfl | rfl | rfl <;> exact (fsdUni_curveDom hn hm _ _ (by decide) true)
  exact h.trans ⟨fun ⟨h,t,ht⟩ => ⟨h,t,(h t).lt_of_ne ht⟩, fun ⟨h,t,ht⟩ => ⟨h,t,ht.ne⟩⟩

theorem claim_C441 (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    (∀ t, 0 < LPM 1 t x → LPMratio 1 t x = 1 / (1 + UPM 1 t x / LPM 1 t x)) ∧
    (sdDomPairCentred x y 1 false = 1 ↔
      (∀ p ∈ Set.Ioo (0 : ℝ) 1, lower 1 p y ≤ lower 1 p x) ∧
      ∃ p ∈ Set.Ioo (0 : ℝ) 1, lower 1 p y < lower 1 p x) :=
  ⟨omega_identity x, strict_expectile_order hn hm x y⟩

theorem claim_C442 (hn : 0 < n) (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    lower 0 p x = quantile7 p x ∧ ∀ d : ℝ, 0 < d → ∀ t,
      t = lower d p x ↔ t ∈ Icc (varMin x) (varMax x) ∧
        (1-p) * rLPM d t x = p * rUPM d t x := by
  refine ⟨(zero_pair hp x).1, fun d hd t => ?_⟩
  have h := root_spec (x := x) hn hd hp
  simp only [balance, sub_eq_zero] at h
  exact ⟨fun e => e ▸ ⟨h.1,h.2.1⟩, fun ⟨ht,he⟩ => h.2.2 t ht he⟩

theorem claim_C443 (hn : 0 < n) (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    upper 0 p x = quantile7 (1-p) x ∧ ∀ d : ℝ, 0 < d → ∀ t,
      t = upper d p x ↔ t ∈ Icc (varMin x) (varMax x) ∧
        p * rLPM d t x = (1-p) * rUPM d t x := by
  refine ⟨(zero_pair hp x).2, fun d hd t => ?_⟩
  simpa [upper] using
    (claim_C442 hn (p := 1-p) ⟨sub_nonneg.mpr hp.2, sub_le_self _ hp.1⟩ x).2 d hd t

theorem claim_C445 (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    lower 0 p x = quantile7 p x ∧
    (∀ d : ℝ, 0 < d → (∃ i j, x i ≠ x j) →
      let t := lower d p x; rLPM d t x / (rLPM d t x + rUPM d t x) = p) ∧
    LPMVaR0 (9/10) ![(0 : ℝ),1,2] = 9/5 ∧
    LPMratio 0 (LPMVaR0 (9/10) ![(0 : ℝ),1,2]) ![(0 : ℝ),1,2] < 9/10 :=
  ⟨(zero_pair hp x).1, fun _ hd hx => lower_share hd hp hx,
    LPMVaR0_type7_example.1, LPMVaR0_type7_example.2.2⟩

theorem claim_C446 (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    let h := ((n-1 : ℕ) : ℝ) * p
    LPMVaR0 p x = sortedAt x ⌊h⌋₊ + (h-⌊h⌋₊) * (sortedAt x ⌈h⌉₊ - sortedAt x ⌊h⌋₊) := by
  rw [LPMVaR0, clampP_of_mem hp]
  exact quantile7_eq hp x

theorem claim_C447 (hn : 0 < n) (hp : p ∈ Set.Icc (0 : ℝ) 1) (t : ℝ) :
    t = lower 1 p x ↔ t ∈ Icc (varMin x) (varMax x) ∧
      (1-p) * LPM 1 t x = p * UPM 1 t x := expectile hn hp t

theorem claim_C450 (hn : 0 < n) (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    upper 0 p x = quantile7 (1-p) x ∧ (∀ d, upper d p x = lower d (1-p) x) ∧
    ∀ d : ℝ, 0 < d → upper d p x = -lower d p (fun i => -x i) := by
  refine ⟨(zero_pair hp x).2, fun _ => rfl, fun d hd => ?_⟩
  let q := upper d p x
  have h := ((claim_C443 hn hp x).2 d hd q).mp rfl
  have hneg := (claim_C442 hn hp (fun i => -x i)).2 d hd (-q)
  have hm : -q ∈ Icc (varMin (fun i => -x i)) (varMax (fun i => -x i)) := by
    rw [varMin_neg hn, varMax_neg hn]
    exact ⟨neg_le_neg h.1.2, neg_le_neg h.1.1⟩
  have hU := rUPM_eq_rLPM_neg hd.ne' q x
  have hL := rUPM_eq_rLPM_neg hd.ne' (-q) (fun i => -x i)
  simp only [neg_neg] at hL
  have he := hneg.mpr ⟨hm, by rw [← hU, hL]; exact h.2.symm⟩
  exact neg_eq_iff_eq_neg.mp he

theorem claim_C492 (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    LPMVaR0 p x = quantile7 p x ∧ UPMVaR0 p x = quantile7 (1-p) x := by
  simp [LPMVaR0, UPMVaR0, clampP_of_mem hp]

theorem claim_C493 (x : Fin n → ℝ) (k : Fin 10) :
    let p : ℝ := (2 * k.val + 1) / 20
    LPMVaR0 p x = quantile7 p x := by
  have hk : (k.val : ℝ) ≤ 9 := by exact_mod_cast (by omega : k.val ≤ 9)
  exact (claim_C492 ⟨by positivity, by linarith⟩ x).1

theorem claim_C553 (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    LPMVaR0 p x = quantile7 p x ∧ LPMVaR0 (1/2) ![(0 : ℝ),0,1] = 0 ∧
    LPMratio 0 0 ![(0 : ℝ),0,1] = 2/3 ∧ (1/2 : ℝ) < 2/3 :=
  ⟨(claim_C492 hp x).1, tied_median.1, tied_median.2, by norm_num⟩

theorem claim_C554 (hn : 0 < n) (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    (∀ t, t = lower 1 p x ↔ t ∈ Icc (varMin x) (varMax x) ∧
      (1-p) * LPM 1 t x = p * UPM 1 t x) ∧
    ((∃ i j, x i ≠ x j) → LPMratio 1 (lower 1 p x) x = p) :=
  ⟨expectile hn hp, fun hx => ratio_at_expectile hx hp⟩

theorem claim_C557 (hn : 0 < n) (hd : 0 < d) (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    (∀ t, t = lower d p x ↔ t ∈ Icc (varMin x) (varMax x) ∧
      (1-p) * rLPM d t x = p * rUPM d t x) ∧
    ((∃ i j, x i ≠ x j) → let t := lower d p x; rLPM d t x / (rLPM d t x + rUPM d t x) = p) ∧
    lower d 0 x = varMin x ∧ lower d 1 x = varMax x := by
  have hle := (varMin_le x ⟨0, hn⟩).trans (le_varMax x _)
  refine ⟨(claim_C442 hn hp x).2 d hd, fun hx => lower_share hd hp hx,
    ((claim_C442 hn ⟨le_rfl, zero_le_one⟩ x).2 d hd _ |>.mpr ⟨⟨le_rfl,hle⟩, ?_⟩).symm,
    ((claim_C442 hn ⟨zero_le_one, le_rfl⟩ x).2 d hd _ |>.mpr ⟨⟨hle,le_rfl⟩, ?_⟩).symm⟩
  · simp [A8_VaRDominance.min_lower hd]
  · simp [A8_VaRDominance.max_upper hd]

theorem claim_C559 (hn : 0 < n) (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    UPMVaR0 p x = quantile7 (1-p) x ∧
    (∀ t, t = upper 1 p x ↔ t ∈ Icc (varMin x) (varMax x) ∧
      p * LPM 1 t x = (1-p) * UPM 1 t x) ∧
    ((∃ i j, x i ≠ x j) → UPMratio 1 (upper 1 p x) x = p) := by
  have hc : 1-p ∈ Set.Icc (0 : ℝ) 1 := ⟨sub_nonneg.mpr hp.2, sub_le_self _ hp.1⟩
  refine ⟨(claim_C492 hp x).2, fun t => ?_, fun hx => ?_⟩
  · simpa [upper] using expectile (x := x) hn hc t
  · rw [upper, UPMratio_eq_one_sub_LPMratio hn, ratio_at_expectile hx hc]
    ring

theorem claim_C609 (hp : p ∈ Set.Icc (0 : ℝ) 1) (x : Fin n → ℝ) :
    LPMVaR0 p x = quantile7 p x ∧ LPMVaR0 (9/10) ![(0 : ℝ),1,2] = 9/5 ∧
    LPMratio 0 (9/5) ![(0 : ℝ),1,2] = 2/3 ∧ (2/3 : ℝ) ≠ 9/10 :=
  ⟨(claim_C492 hp x).1, LPMVaR0_type7_example.1, LPMVaR0_type7_example.2.1, by norm_num⟩

theorem claim_C588 (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    let X := fun t => LPM 0 t x
    let Y := fun t => LPM 0 t y
    (nnsFSD x y true = .xDom ↔ CurveDom X Y) ∧
    (nnsFSD x y true = .yDom ↔ CurveDom Y X) ∧
    (nnsFSD x y true = .noDom ↔ ¬ CurveDom X Y ∧ ¬ CurveDom Y X) := by
  dsimp only
  have hx := nnsFSD_eq_xDom_iff hn hm x y true
  have hy := nnsFSD_eq_yDom_iff hn hm x y true
  refine ⟨hx, hy, ?_⟩
  change _ ↔ ¬ CurveDom (domCurve 1 true x) (domCurve 1 true y) ∧
    ¬ CurveDom (domCurve 1 true y) (domCurve 1 true x)
  rw [← hx, ← hy]
  cases nnsFSD x y true <;> simp

theorem claim_C594 (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    sdDomPairCentred x y 1 true = 1 ↔
      (∀ t, IsKnot x y t → (#{i | x i ≤ t} : ℝ) / n ≤ (#{j | y j ≤ t} : ℝ) / m) ∧
        ∃ t, IsKnot x y t ∧ (#{i | x i ≤ t} : ℝ) / n < (#{j | y j ≤ t} : ℝ) / m := by
  simp only [← LPM_zero_eq_ecdf]
  rw [claim_C440 hn hm x y 1 (Or.inl rfl)]
  change ((∀ t, LPM 0 t x ≤ LPM 0 t y) ∧ ∃ t, LPM 0 t x < LPM 0 t y) ↔ _
  constructor
  · rintro ⟨h,t,ht⟩
    refine ⟨fun a _ => h a, ?_⟩
    by_contra! he
    have hr := LPM_zero_le_of_le_at_knots (x := y) (y := x) (fun a ha => he a (Or.symm ha)) t
    exact ht.not_ge hr
  · rintro ⟨h,t,_,ht⟩
    exact ⟨LPM_zero_le_of_le_at_knots h,t,ht⟩

namespace A8_VaRDominance

def bookX : Fin 4 → ℝ := ![0,1,4,5]
def bookY : Fin 4 → ℝ := ![0,2,5/2,3]

private theorem book_cross (t : ℝ) :
    LPM 1 t bookX * UPM 1 t bookY ≤ LPM 1 t bookY * UPM 1 t bookX := by
  simp only [LPM, UPM, Fin.sum_univ_four, bookX, bookY, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.head_cons,
    Matrix.tail_cons, Nat.cast_ofNat]
  have below (v : ℝ) (hv : v ≤ t) := And.intro (lowerDev_one_of_le hv) (upperDev_one_of_le hv)
  have above (v : ℝ) (hv : t ≤ v) := And.intro (lowerDev_one_of_ge hv) (upperDev_one_of_ge hv)
  rcases le_total t 0 with h0 | h0
  · simp only [above 0 (by linarith),
      above 1 (by linarith),
      above 2 (by linarith),
      above (5/2) (by linarith),
      above 3 (by linarith),
      above 4 (by linarith),
      above 5 (by linarith)]
    ring_nf
    exact le_rfl
  rcases le_total t 1 with h1 | h1
  · simp only [below 0 (by linarith),
      above 1 (by linarith),
      above 2 (by linarith),
      above (5/2) (by linarith),
      above 3 (by linarith),
      above 4 (by linarith),
      above 5 (by linarith)]
    linarith
  rcases le_total t 2 with h2 | h2
  · simp only [below 0 (by linarith),
      below 1 (by linarith),
      above 2 (by linarith),
      above (5/2) (by linarith),
      above 3 (by linarith),
      above 4 (by linarith),
      above 5 (by linarith)]
    nlinarith [sq_nonneg (t-9/8)]
  rcases le_total t (5/2) with h3 | h3
  · simp only [below 0 (by linarith),
      below 1 (by linarith),
      below 2 (by linarith),
      above (5/2) (by linarith),
      above 3 (by linarith),
      above 4 (by linarith),
      above 5 (by linarith)]
    linarith
  rcases le_total t 3 with h4 | h4
  · simp only [below 0 (by linarith),
      below 1 (by linarith),
      below 2 (by linarith),
      below (5/2) (by linarith),
      above 3 (by linarith),
      above 4 (by linarith),
      above 5 (by linarith)]
    nlinarith [mul_nonneg (sub_nonneg.mpr h3) (sub_nonneg.mpr h4)]
  rcases le_total t 4 with h5 | h5
  · simp only [below 0 (by linarith),
      below 1 (by linarith),
      below 2 (by linarith),
      below (5/2) (by linarith),
      below 3 (by linarith),
      above 4 (by linarith),
      above 5 (by linarith)]
    nlinarith [mul_nonneg (sub_nonneg.mpr h4) (sub_nonneg.mpr h5)]
  rcases le_total t 5 with h6 | h6
  · simp only [below 0 (by linarith),
      below 1 (by linarith),
      below 2 (by linarith),
      below (5/2) (by linarith),
      below 3 (by linarith),
      below 4 (by linarith),
      above 5 (by linarith)]
    nlinarith [mul_nonneg (sub_nonneg.mpr h5) (sub_nonneg.mpr h6)]
  · simp only [below 0 (by linarith),
      below 1 (by linarith),
      below 2 (by linarith),
      below (5/2) (by linarith),
      below 3 (by linarith),
      below 4 (by linarith),
      below 5 (by linarith)]
    ring_nf
    exact le_rfl

private theorem book_total (s : Fin 4 → ℝ) (h : s 0 ≠ s 1) (t : ℝ) :
    0 < LPM 1 t s + UPM 1 t s := by
  apply ratio_denom_pos
  by_contra! he
  exact h ((he 0).trans (he 1).symm)

theorem book_area : CurveDom (fun t => LPMratio 1 t bookX) (fun t => LPMratio 1 t bookY) := by
  have hx := book_total bookX (by norm_num [bookX])
  have hy := book_total bookY (by norm_num [bookY])
  constructor
  · intro t
    change LPMratio 1 t bookX ≤ LPMratio 1 t bookY
    rw [LPMratio_one_of_pos (hx t), LPMratio_one_of_pos (hy t), div_le_div_iff₀ (hx t) (hy t)]
    linarith [book_cross t]
  · refine ⟨1, ?_⟩
    change LPMratio 1 1 bookX ≠ LPMratio 1 1 bookY
    rw [LPMratio_one_of_pos (hx 1), LPMratio_one_of_pos (hy 1)]
    norm_num [LPM, UPM, lowerDev, upperDev, Fin.sum_univ_four, bookX, bookY]

end A8_VaRDominance

theorem claim_C342 :
    CurveDom (fun t => LPMratio 1 t bookX) (fun t => LPMratio 1 t bookY) ∧
    LPM 0 1 bookY < LPM 0 1 bookX ∧ LPM 0 3 bookX < LPM 0 3 bookY ∧
    Monotone (fun v : ℝ => if 1 < v then (1 : ℝ) else 0) ∧
    mean (fun i => if 1 < bookX i then (1 : ℝ) else 0) = 1/2 ∧
    mean (fun i => if 1 < bookY i then (1 : ℝ) else 0) = 3/4 ∧
    nnsFSD bookX bookY false = .xDom ∧ nnsFSD bookX bookY true = .noDom ∧
    sdDomPairCentred bookX bookY 1 false = 1 ∧ sdDomPairCentred bookX bookY 1 true = 0 ∧
    ∀ {n m : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ) (discrete : Bool),
      sdDomPairCentred x y 1 discrete = 0 ∨ sdDomPairCentred x y 1 discrete = 1 := by
  have h1 : LPM 0 1 bookY < LPM 0 1 bookX := by
    simp only [LPM, lowerDev, Fin.sum_univ_four]; norm_num [bookX, bookY]
  have h3 : LPM 0 3 bookX < LPM 0 3 bookY := by
    simp only [LPM, lowerDev, Fin.sum_univ_four]; norm_num [bookX, bookY]
  refine ⟨book_area,h1,h3,?_,?_,?_,(nnsFSD_eq_xDom_iff (by decide) (by decide) _ _ _).mpr book_area,
    (claim_C588 (by decide) (by decide) _ _).2.2.mpr ⟨fun h => h1.not_ge (h.1 1),
      fun h => h3.not_ge (h.1 3)⟩,
    (fsdUni_curveDom (by decide) (by decide) _ _ (Or.inl rfl) false).mpr book_area, ?_, ?_⟩
  · intro a b hab
    dsimp only
    split_ifs <;> norm_num
    linarith
  · simp only [mean, Fin.sum_univ_four]; norm_num [bookX]
  · simp only [mean, Fin.sum_univ_four]; norm_num [bookY]
  · refine (sdDomPairCentred_fsd_zero_or_one _ _).resolve_right fun h => ?_
    exact h1.not_ge (((fsdUni_curveDom (by decide) (by decide) _ _ (Or.inl rfl) true).mp h).1 1)
  · intro n m x y discrete
    unfold sdDomPairCentred
    split_ifs <;> first | exact Or.inl rfl | exact sdDomPair_eq_zero_or_one _ _ _ _

theorem claim_C341 (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    (∀ t, 0 < LPM 1 t x → LPMratio 1 t x = 1 / (1 + UPM 1 t x / LPM 1 t x)) ∧
    ((∀ t, LPMratio 1 t x ≤ LPMratio 1 t y) ↔
      ∀ p ∈ Set.Ioo (0 : ℝ) 1, lower 1 p y ≤ lower 1 p x) ∧
    sdDomPairCentred c591X c591Y 1 false = 1 ∧ sdDomPairCentred c591X c591Y 1 true = 0 ∧
    ∃ u : ℝ → ℝ, Monotone u ∧ mean (u ∘ c591X) < mean (u ∘ c591Y) :=
  ⟨omega_identity x, expectile_order hn hm x y, c591_omega, c591_not_fsd,
    _, c591_utility.1, c591_utility.2⟩

theorem claim_C590 (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    (sdDomPairCentred x y 1 false = 1 ↔
      CurveDom (fun t => LPMratio 1 t x) (fun t => LPMratio 1 t y)) ∧
    (sdDomPairCentred x y 1 false = 1 ↔
      (∀ p ∈ Set.Ioo (0 : ℝ) 1, lower 1 p y ≤ lower 1 p x) ∧
      ∃ p ∈ Set.Ioo (0 : ℝ) 1, lower 1 p y < lower 1 p x) ∧
    ∀ t, 0 < LPM 1 t x → 0 < LPM 1 t y →
      (LPMratio 1 t x ≤ LPMratio 1 t y ↔ UPM 1 t y / LPM 1 t y ≤ UPM 1 t x / LPM 1 t x) := by
  refine ⟨fsdUni_curveDom hn hm x y (Or.inl rfl) false, strict_expectile_order hn hm x y, ?_⟩
  intro t hx hy
  rw [omega_identity x t hx, omega_identity y t hy, div_le_div_iff₀ (by
    have := UPM_nonneg 1 t x; positivity) (by have := UPM_nonneg 1 t y; positivity)]
  simp

namespace A8_VaRDominance

theorem area_le_of_fsd (hn : 0 < n) (hm : 0 < m)
    (h : ∀ t, LPM 0 t x ≤ LPM 0 t y) (t : ℝ) : LPMratio 1 t x ≤ LPMratio 1 t y := by
  have hu : Monotone (upperDev 1 t) := by
    intro a b hab
    simp only [upperDev_eq_max_pow one_ne_zero, pow_one]
    exact max_le_max_right _ (sub_le_sub_right hab t)
  have hl : Monotone (fun v => -lowerDev 1 t v) := by
    intro a b hab
    simp only [lowerDev_eq_max_pow one_ne_zero, pow_one, neg_le_neg_iff]
    exact max_le_max_right _ (sub_le_sub_left hab t)
  have hU : UPM 1 t y ≤ UPM 1 t x := mean_comp_le_of_LPM_zero_le hn hm h hu
  have hL := mean_comp_le_of_LPM_zero_le hn hm h hl
  simp only [mean, sum_neg_distrib, neg_div, neg_le_neg_iff] at hL
  change LPM 1 t x ≤ LPM 1 t y at hL
  by_cases hy : ∀ i, y i = t
  · rw [LPMratio_of_forall_eq one_ne_zero hm hy]
    exact (LPMratio_mem_Icc 1 t x).2
  push Not at hy
  have hpy := ratio_denom_pos (d := 1) hy
  rw [LPMratio_one_of_pos hpy]
  by_cases hx : ∀ i, x i = t
  · have hx0 : UPM 1 t x = 0 := by simp [UPM, upperDev, hx]
    have hy0 := le_antisymm (hU.trans_eq hx0) (UPM_nonneg 1 t y)
    rw [hy0, add_zero] at hpy ⊢
    rw [div_self hpy.ne']
    exact (LPMratio_mem_Icc 1 t x).2
  push Not at hx
  have hpx := ratio_denom_pos (d := 1) hx
  rw [LPMratio_one_of_pos hpx, div_le_div_iff₀ hpx hpy]
  linarith only [mul_le_mul hL hU (UPM_nonneg 1 t y) (LPM_nonneg 1 t y)]

end A8_VaRDominance

theorem claim_C607 (hn : 0 < n) (hm : 0 < m) (x : Fin n → ℝ) (y : Fin m → ℝ) :
    ((∀ t, LPMratio 1 t x ≤ LPMratio 1 t y) ↔ ∀ p ∈ Set.Ioo (0 : ℝ) 1, lower 1 p y ≤ lower 1 p x) ∧
    (∀ t, 0 < LPM 1 t x → 0 < LPM 1 t y → (LPMratio 1 t x ≤ LPMratio 1 t y ↔
      UPM 1 t y / LPM 1 t y ≤ UPM 1 t x / LPM 1 t x)) ∧
    ((∀ t, LPM 0 t x ≤ LPM 0 t y) → ∀ t, LPMratio 1 t x ≤ LPMratio 1 t y) ∧
    nnsFSD bookX bookY false = .xDom ∧ nnsFSD bookX bookY true = .noDom := by
  obtain ⟨_, _, _, _, _, _, hx, hy, _⟩ := claim_C342
  exact ⟨expectile_order hn hm x y, (claim_C590 hn hm x y).2.2, area_le_of_fsd hn hm, hx, hy⟩

theorem claim_C494 (z : Fin 3000 → ℝ) (k : Fin 10) :
    let p : ℝ := (2 * k.val + 1) / 20
    let j := 150 + 300 * k.val
    LPMratio 0 (LPMVaR0 p z) z = p ↔ sortedAt z (j-1) < sortedAt z j := by
  let p : ℝ := (2 * k.val + 1) / 20
  let j := 150 + 300 * k.val
  have hk : k.val ≤ 9 := by omega
  have hj : j < 3000 := by dsimp [j]; omega
  have hj0 : 1 ≤ j := by dsimp [j]; omega
  have hp : p ∈ Set.Ioo (0 : ℝ) 1 := by
    have hk' : (k.val : ℝ) ≤ 9 := by exact_mod_cast hk
    constructor <;> dsimp [p] <;> first | positivity | linarith
  have hjp : (j : ℝ) = 3000*p := by simp [j,p]; ring
  have hjpred : ((j-1 : ℕ) : ℝ) = (j : ℝ)-1 := by exact_mod_cast Nat.cast_sub (R := ℝ) hj0
  have hfloor : ⌊(2999:ℝ)*p⌋₊ = j-1 := by
    rw [Nat.floor_eq_iff (by exact mul_nonneg (by norm_num) hp.1.le)]
    rw [hjpred]; constructor <;> linarith [hp.1,hp.2]
  have hceil : ⌈(2999:ℝ)*p⌉₊ = j := by
    rw [Nat.ceil_eq_iff (by omega)]
    constructor <;> linarith [hp.1,hp.2]
  let q := LPMVaR0 p z
  have hq : q = sortedAt z (j-1) + (1-p)*(sortedAt z j - sortedAt z (j-1)) := by
    dsimp only [q]
    rw [LPMVaR0, clampP_of_mem ⟨hp.1.le,hp.2.le⟩, quantile7_eq ⟨hp.1.le,hp.2.le⟩]
    norm_num only [Nat.reduceSub, Nat.cast_ofNat]
    rw [hfloor,hceil,hjpred]
    congr 2
    linarith
  have hs := sortedAt_sorted z (j-1) j (by omega) hj
  have hlo : sortedAt z (j-1) ≤ q := by
    rw [hq]; linarith [mul_nonneg (sub_nonneg.mpr hp.2.le) (sub_nonneg.mpr hs)]
  have hlo' := (lt_countLe_iff (sortedAt_sorted z) q (by omega : j-1 < 3000)).mpr hlo
  change LPMratio 0 q z = p ↔ sortedAt z (j-1) < sortedAt z j
  rw [LPMratio_zero_eq_ecdf, ← countLe_sortedAt, Nat.cast_ofNat,
    div_eq_iff (by norm_num : (3000 : ℝ) ≠ 0)]
  constructor
  · intro he
    by_contra! hge
    have heq := le_antisymm hs hge
    have hqj : sortedAt z j ≤ q := heq ▸ hlo
    have hc := (lt_countLe_iff (sortedAt_sorted z) q hj).mpr hqj
    have hcr : (j : ℝ) < (countLe (sortedAt z) 3000 q : ℝ) := by exact_mod_cast hc
    linarith
  · intro hgap
    have hlt : q < sortedAt z j := by
      rw [hq]; linarith only [mul_pos hp.1 (sub_pos.mpr hgap)]
    have hu := mt (lt_countLe_iff (sortedAt_sorted z) q hj).mp hlt.not_ge
    have heq : countLe (sortedAt z) 3000 q = j := by omega
    rw [heq]
    linarith

end
end NNS.ClaimsA
