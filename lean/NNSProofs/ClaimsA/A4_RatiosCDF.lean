-- Copyright (c) 2026 NNS contributors. GPL-3.0-only.
module

public import NNSProofs.GapClaims

namespace NNS.ClaimsA

open Finset MeasureTheory ProbabilityTheory Filter Set Topology
open GapOperators GapQuantile

@[expose] public noncomputable section

variable {n : ℕ}

namespace A4_RatiosCDF

def lowerShare (d t : ℝ) (x : Fin n → ℝ) : ℝ :=
  if d = 0 then LPM 0 t x else if 0 < n ∧ ∀ i, x i = t then 1
  else rLPM d t x / (rLPM d t x + rUPM d t x)

def upperShare (d t : ℝ) (x : Fin n → ℝ) : ℝ :=
  if d = 0 then UPM 0 t x else if 0 < n ∧ ∀ i, x i = t then 0
  else rUPM d t x / (rLPM d t x + rUPM d t x)

theorem lowerShare_nat (d : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    lowerShare d t x = LPMratio d t x := by
  simp [lowerShare, LPMratio, rLPM_natCast, rUPM_natCast]

private theorem total_zero (hn : 0 < n) {d : ℝ} (hd : 0 < d) (t : ℝ) (x : Fin n → ℝ) :
    rLPM d t x + rUPM d t x = 0 ↔ ∀ i, x i = t := by
  rw [rLPM_add_rUPM, div_eq_zero_iff, or_iff_left (Nat.cast_ne_zero.mpr hn.ne'),
    Finset.sum_eq_zero_iff_of_nonneg (fun i _ => Real.rpow_nonneg (abs_nonneg _) d)]
  simp [Real.rpow_eq_zero (abs_nonneg _) hd.ne', sub_eq_zero, eq_comm]

def lowerFrequency (t : ℝ) (x : Fin n → ℝ) : ℝ := (#{i | x i ≤ t} : ℝ) / n

def upperFrequency (t : ℝ) (x : Fin n → ℝ) : ℝ := (#{i | t < x i} : ℝ) / n

def lowerFormula (d t : ℝ) (x : Fin n → ℝ) : ℝ :=
  if rLPM d t x + rUPM d t x = 0 then lowerFrequency t x
  else rLPM d t x / (rLPM d t x + rUPM d t x)

def upperFormula (d t : ℝ) (x : Fin n → ℝ) : ℝ :=
  if rLPM d t x + rUPM d t x = 0 then upperFrequency t x
  else rUPM d t x / (rLPM d t x + rUPM d t x)

private theorem freq_zero (t : ℝ) (x : Fin n → ℝ) :
    lowerShare 0 t x = lowerFrequency t x ∧ upperShare 0 t x = upperFrequency t x :=
  ⟨by simpa [lowerShare, lowerFrequency] using LPM_zero_eq_ecdf t x,
    by simpa [upperShare, UPMratio, upperFrequency] using UPMratio_zero_eq_survival t x⟩

private theorem shares (hn : 0 < n) {d : ℝ} (hd : 0 ≤ d) (t : ℝ) (x : Fin n → ℝ) :
    lowerShare d t x = lowerFormula d t x ∧ upperShare d t x = upperFormula d t x := by
  rcases hd.eq_or_lt with h | hp
  · subst d
    simp [lowerShare, upperShare, lowerFormula, upperFormula,
      rLPM_zero, rUPM_zero, LPM_zero_add_UPM_zero hn]
  · by_cases ht : ∀ i, x i = t
    · have hz := (total_zero hn hp t x).2 ht
      simp [lowerShare, upperShare, lowerFormula, upperFormula, lowerFrequency, upperFrequency,
        hp.ne', hn, ht, hz, Nat.cast_ne_zero.mpr hn.ne']
    · simp [lowerShare, upperShare, lowerFormula, upperFormula, hp.ne', ht,
        (total_zero hn hp t x).not.mpr ht]

private theorem share_partition (hn : 0 < n) {d : ℝ} (hd : 0 ≤ d) (t : ℝ) (x : Fin n → ℝ) :
    lowerShare d t x + upperShare d t x = 1 := by
  rw [(shares hn hd t x).1, (shares hn hd t x).2]
  unfold lowerFormula upperFormula
  split_ifs with h
  · rw [← (freq_zero t x).1, ← (freq_zero t x).2]
    simpa [lowerShare, upperShare] using LPM_zero_add_UPM_zero hn t x
  · rw [← add_div, div_self h]

private theorem induced {F : ℝ → ℝ} (hm : Monotone F)
    (hc : ∀ t, ContinuousWithinAt F (Ici t) t)
    (h0 : Tendsto F atBot (𝓝 0)) (h1 : Tendsto F atTop (𝓝 1)) :
    ∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧ ∀ t, cdf ν t = F t := by
  let f : StieltjesFunction ℝ := ⟨F, hm, hc⟩
  refine ⟨f.measure, ⟨?_⟩, fun t => ?_⟩
  · rw [f.measure_univ h0 h1, sub_zero, ENNReal.ofReal_one]
  · exact congrArg (fun g : StieltjesFunction ℝ => g t) (cdf_measure_stieltjesFunction f h0 h1)

private theorem curve_one (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    cdfCurve 1 t x = if LPM 1 t x + UPM 1 t x = 0 then 1
      else LPM 1 t x / (LPM 1 t x + UPM 1 t x) := by
  simp [cdfCurve, ratio_denom_eq_zero_iff one_ne_zero, hn]

private theorem frequency_differs :
    cdfCurve 1 0 ![(0 : ℝ), 1] = 0 ∧ lowerFrequency 0 ![(0 : ℝ), 1] = 1 / 2 := by
  constructor
  · apply cdfCurve_eq_zero_of_le_min one_ne_zero
    · intro i
      fin_cases i <;> norm_num
    · exact ⟨1, by norm_num⟩
  · unfold lowerFrequency
    rw [← LPMratio_zero_eq_ecdf]
    simp only [LPMratio, ↓reduceIte, LPM, Fin.sum_univ_two, lowerDev]
    norm_num

end A4_RatiosCDF

open A4_RatiosCDF

theorem claim_C018 (hn : 0 < n) {d : ℝ} (hd : 0 ≤ d) (t : ℝ) (x : Fin n → ℝ)
    (ht : 0 < rLPM d t x + rUPM d t x) :
    lowerShare d t x = rLPM d t x / (rLPM d t x + rUPM d t x) ∧
    upperShare d t x = rUPM d t x / (rLPM d t x + rUPM d t x) := by
  simpa only [lowerFormula, upperFormula, ite_eq_right ht.ne'] using shares hn hd t x

theorem claim_C038 (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ) :
    mUPM 0 t μ = μ.real (Ioi t) := mUPM_zero_eq t

theorem claim_C041 (μ : Measure ℝ) [IsProbabilityMeasure μ] (t : ℝ) :
    1 - cdf μ t = μ.real (Ioi t) ∧ μ.real (Ioi t) = mUPM 0 t μ ∧
      cdf μ t + mUPM 0 t μ = 1 ∧
      ∀ {m : ℕ}, 0 < m → ∀ (x : Fin m → ℝ), survivalCurve 0 t x = UPMratio 0 t x := by
  have hp := mLPM_zero_add_mUPM_zero (μ := μ) t
  rw [mLPM_zero_eq_cdf] at hp
  exact ⟨by rw [← mUPM_zero_eq]; linarith, (mUPM_zero_eq t).symm, hp,
    fun hn x => (UPMratio_eq_one_sub_LPMratio hn 0 t x).symm⟩

theorem claim_C262 (μ : Measure ℝ) [IsProbabilityMeasure μ] (h : Integrable (fun v : ℝ => v) μ)
    (hx : ∀ t, μ {t}ᶜ ≠ 0) :
    let F := fun t => mPMRatio 1 t μ
    Continuous F ∧ Monotone F ∧ StrictMonoOn F {t | 0 < F t ∧ F t < 1} ∧
      Tendsto F atBot (𝓝 0) ∧ Tendsto F atTop (𝓝 1) ∧
      (∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧ ∀ t, cdf ν t = F t) ∧
      ¬ StrictMono (fun t => mPMRatio 1 t (empMeasure ![(0 : ℝ), 1])) := by
  have h1 : MemLp (fun v : ℝ => v) (1 : ℕ) μ := by simpa using (memLp_one_iff_integrable.mpr h)
  have hc := continuous_mPMRatio (r := 1) one_ne_zero h1 hx
  have hm := mPMRatio_monotone (r := 1) h1
  have h0 := tendsto_mPMRatio_atBot (r := 1) h1
  have ht := tendsto_mPMRatio_atTop (r := 1) h1
  exact ⟨hc, hm, fun _ hs _ _ hst => mPMRatio_lt_mPMRatio one_ne_zero h1 hst hs.1 hs.2,
    h0, ht, induced hm (fun _ => hc.continuousAt.continuousWithinAt) h0 ht, by
      intro hs
      have h := hs (a := -1) (b := 0) (by norm_num)
      have hnc : ∃ i j, ![(0 : ℝ), 1] i ≠ ![0, 1] j := ⟨0, 1, by norm_num⟩
      simp only [mPMRatio_empMeasure one_ne_zero hnc] at h
      have hlow : LPMratio 1 (-1) ![(0 : ℝ), 1] = 0 :=
        LPMratio_eq_zero_of_forall_lt _ (fun i => by fin_cases i <;> norm_num)
      have hz : LPMratio 1 0 ![(0 : ℝ), 1] = 0 := frequency_differs.1
      linarith⟩

theorem claim_C263 {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j) :
    cdfCurve 1 (mean x) x = 1 / 2 := LPMratio_one_mean hx

theorem claim_C267 (hn : 0 < n) {c : ℝ} {x : Fin n → ℝ} (hx : ∀ i, x i = c) :
    LPM 1 c x = 0 ∧ UPM 1 c x = 0 ∧ ∀ t, cdfCurve 1 t x = if c ≤ t then 1 else 0 := by
  have hz := (ratio_denom_eq_zero_iff one_ne_zero c x).2 hx
  exact ⟨by linarith [UPM_nonneg 1 c x, LPM_nonneg 1 c x],
    by linarith [UPM_nonneg 1 c x, LPM_nonneg 1 c x], cdfCurve_of_const hn 1 hx⟩

theorem claim_C268 (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    cdf (empMeasure x) t = LPM 0 t x ∧ LPM 0 t x = LPMratio 0 t x ∧
      LPMratio 0 t x = (∑ i, if x i ≤ t then (1 : ℝ) else 0) / n :=
  ⟨cdf_empMeasure hn x t, rfl, by rw [Finset.sum_boole]; exact LPMratio_zero_eq_ecdf t x⟩

theorem claim_C269 (hn : 0 < n) (x : Fin n → ℝ) :
    (∀ t, cdf (empMeasure x) t = cdfCurve 0 t x) ∧ empMeasure x (Set.range x)ᶜ = 0 ∧
    (∀ s t, s ≤ t → (∀ i, x i ≤ t → x i ≤ s) → cdfCurve 0 s x = cdfCurve 0 t x) ∧
    (∀ t, Tendsto (fun s => cdfCurve 0 s x) (𝓝[<] t) (𝓝 ((#{i | x i < t} : ℝ) / n)) ∧
      cdfCurve 0 t x = (#{i | x i < t} : ℝ) / n + (#{i | x i = t} : ℝ) / n) ∧
    ∀ i, ¬ ContinuousAt (fun t => cdfCurve 0 t x) (x i) := by
  refine ⟨fun t => cdf_empMeasure hn x t, ?_,
    fun _ _ h hgap => (cdfCurve_zero_eq_of_no_obs x hgap h).symm,
    fun t => ⟨cdfCurve_zero_tendsto_left t x, cdfCurve_zero_eq_lt_add_eq t x⟩, ?_⟩
  · rw [empMeasure, Measure.smul_apply, Measure.coe_finsetSum, Finset.sum_apply]
    have hm := (Set.finite_range x).measurableSet.compl
    simp [Measure.dirac_apply' _ hm]
  · intro i hc
    have h := tendsto_nhds_unique (cdfCurve_zero_tendsto_left (x i) x)
      (hc.tendsto.mono_left nhdsWithin_le_nhds)
    rw [cdfCurve_zero_eq_lt_add_eq] at h
    have hp : 0 < #{j | x j = x i} := Finset.card_pos.mpr ⟨i, by simp⟩
    have : 0 < (#{j | x j = x i} : ℝ) / n := by positivity
    linarith

theorem claim_C273 (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    LPM 1 t x = (∑ i, max (t - x i) 0) / n ∧
    UPM 1 t x = (∑ i, max (x i - t) 0) / n ∧ cdfCurve 1 t x = LPMratio 1 t x ∧
    cdfCurve 1 t x = if LPM 1 t x + UPM 1 t x = 0 then 1
      else LPM 1 t x / (LPM 1 t x + UPM 1 t x) :=
  ⟨by simp [LPM, lowerDev_one], by simp [UPM, upperDev_one], rfl, curve_one hn t x⟩

theorem claim_C276 (hn : 0 < n) (x : Fin n → ℝ) :
    (Continuous (fun t => cdfCurve 1 t x) ↔ ∃ i j, x i ≠ x j) ∧
    ∀ i, ¬ ContinuousAt (fun t => cdfCurve 0 t x) (x i) := by
  refine ⟨⟨fun hc => ?_, cdfCurve_continuous one_ne_zero⟩, (claim_C269 hn x).2.2.2.2⟩
  by_contra h
  push Not at h
  exact not_continuousAt_cdfCurve_of_const hn 1 (fun i => h i ⟨0, hn⟩) hc.continuousAt

theorem claim_C292 {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j) :
    Continuous (fun t => LPMratio 1 t x) ∧ LPMratio 1 (mean x) x = 1 / 2 :=
  ⟨cdfCurve_continuous one_ne_zero hx, LPMratio_one_mean hx⟩

theorem claim_C294 (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    cdfCurve 1 t x = if LPM 1 t x + UPM 1 t x = 0 then 1
      else LPM 1 t x / (LPM 1 t x + UPM 1 t x) := curve_one hn t x

theorem claim_C295 {x : Fin n → ℝ} (hx : ∃ i j, x i ≠ x j) :
    cdfCurve 1 (mean x) x = 1 / 2 := LPMratio_one_mean hx

theorem claim_C360 (μ : Measure ℝ) [IsProbabilityMeasure μ] {d : ℝ} (hd : 1 ≤ d) (t : ℝ)
    (hI : Integrable (fun v : ℝ => |v| ^ d) μ) (ht : μ {t}ᶜ ≠ 0) :
    0 < rlpm d t μ + rupm d t μ ∧ GapClaims.rPMRatio d t μ =
      (∫ v, max (t - v) 0 ^ d ∂μ) /
        ((∫ v, max (t - v) 0 ^ d ∂μ) + ∫ v, max (v - t) 0 ^ d ∂μ) := by
  have hp : 0 < d := lt_of_lt_of_le zero_lt_one hd
  refine ⟨GapClaims.rlpm_add_rupm_pos hp hI ht, ?_⟩
  rw [GapClaims.rPMRatio, rlpm_eq_integral_max_rpow hp, rupm_eq_integral_max_rpow hp.ne']

theorem claim_C427 (hn : 0 < n) {d : ℝ} (hd : 0 ≤ d) (t : ℝ) (x : Fin n → ℝ) :
    lowerShare 0 t x = lowerFrequency t x ∧ lowerShare d t x = lowerFormula d t x ∧
    (let y : Fin 1 → ℝ := fun _ => t
     LPMratio 0 t y = 1 ∧ (#{i | y i < t} : ℝ) / 1 = 0) :=
  ⟨(freq_zero t x).1, (shares hn hd t x).1, by simp [LPMratio, LPM, lowerDev]⟩

theorem claim_C428 (hn : 0 < n) {d : ℝ} (hd : 0 ≤ d) (t : ℝ) (x : Fin n → ℝ)
    (ht : 0 < rLPM d t x + rUPM d t x) :
    upperShare 0 t x = upperFrequency t x ∧
    upperShare d t x = rUPM d t x / (rLPM d t x + rUPM d t x) :=
  ⟨(freq_zero t x).2, (claim_C018 hn hd t x ht).2⟩

theorem claim_C429 (hn : 0 < n) {d : ℝ} (hd : 0 ≤ d) (t : ℝ) (x : Fin n → ℝ) :
    (lowerShare 0 t x = lowerFrequency t x ∧ upperShare 0 t x = upperFrequency t x) ∧
    lowerShare d t x = lowerFormula d t x ∧ upperShare d t x = upperFormula d t x :=
  ⟨freq_zero t x, shares hn hd t x⟩

namespace A4_RatiosCDF

def cdfRow (d : ℝ) (x : Fin n → ℝ) (i : Fin n) : ℝ × ℝ :=
  let t := x (Tuple.sort x i)
  (t, lowerShare d t x)

def jointShare {k : ℕ} (d : ℝ) (X : Matrix (Fin n) (Fin k) ℝ) (t : Fin k → ℝ) : ℝ :=
  if d = 0 ∨ normConstRaw X t d = 0 then clpmRaw X t 0
  else clpmRaw X t d / normConstRaw X t d

end A4_RatiosCDF

open scoped Classical in
theorem claim_C434 {d : ℝ} (hd : 0 ≤ d) (t : ℝ) (x : Fin n → ℝ) (i : Fin n) :
    cdfRow d x i = (x (Tuple.sort x i), lowerFormula d (x (Tuple.sort x i)) x) ∧
    lowerShare d t x = lowerFormula d t x ∧
    (∀ k : ℕ, cdfTable k x i = cdfRow k x i ∧ cdfTarget k t x = LPMratio k t x) ∧
    ∀ {k : ℕ} (X : Matrix (Fin n) (Fin k) ℝ) (u : Fin k → ℝ), jointShare d X u =
      if d = 0 ∨ totalND X u d = 0 then (#{j | ∀ col, X j col ≤ u col} : ℝ) / n
      else clpmND X u d / totalND X u d := by
  refine ⟨Prod.ext rfl (shares (Fin.pos i) hd _ x).1, (shares (Fin.pos i) hd t x).1,
    fun k => ⟨Prod.ext rfl (lowerShare_nat k _ x).symm, rfl⟩, fun X u => ?_⟩
  rw [jointShare, normConstRaw_eq hd, clpmRaw_eq hd]
  split_ifs
  · simpa only [coLPMnD_R_v8898017, ↓reduceIte] using coLPMnD_R_zero_v8898017 X u false
  · rfl

theorem claim_C436 (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) :
    LPMratio 0 t x = cdf (empMeasure x) t ∧ LPMratio 0 t x = lowerFrequency t x :=
  ⟨(cdf_empMeasure hn x t).symm, LPMratio_zero_eq_ecdf t x⟩

theorem claim_C437 (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    (LPMratio 1 t x = if LPM 1 t x + UPM 1 t x = 0 then 1
      else LPM 1 t x / (LPM 1 t x + UPM 1 t x)) ∧
    (cdfCurve 1 0 ![(0 : ℝ), 1] = 0 ∧ lowerFrequency 0 ![(0 : ℝ), 1] = 1 / 2) :=
  ⟨curve_one hn t x, frequency_differs⟩

theorem claim_C473 (x : Fin n → ℝ) (k : Fin n) :
    (let t := x (Tuple.sort x k)
     cdfTable 1 x k = (t, if LPM 1 t x + UPM 1 t x = 0 then 1
       else LPM 1 t x / (LPM 1 t x + UPM 1 t x))) ∧
    (cdfCurve 1 0 ![(0 : ℝ), 1] = 0 ∧ lowerFrequency 0 ![(0 : ℝ), 1] = 1 / 2) :=
  ⟨Prod.ext rfl (curve_one (Fin.pos k) _ _), frequency_differs⟩

theorem claim_C474 (hn : 0 < n) (x : Fin n → ℝ) :
    cdfTarget 1 (mean x) x = LPMratio 1 (mean x) x ∧
    cdfTarget 1 (mean x) x = if ∃ i j, x i ≠ x j then 1 / 2 else 1 := by
  refine ⟨rfl, ?_⟩
  split_ifs with hx
  · exact LPMratio_one_mean hx
  · push Not at hx
    exact LPMratio_mean_of_const hn 1 hx

theorem claim_C475 (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    survivalCurve 1 t x = 1 - cdfCurve 1 t x ∧ survivalCurve 1 t x = UPMratio 1 t x :=
  ⟨rfl, (UPMratio_eq_one_sub_LPMratio hn 1 t x).symm⟩

theorem claim_C487 (t : ℝ) (x : Fin n → ℝ) :
    LPMratio 0 t x = (∑ i, if x i ≤ t then (1 : ℝ) else 0) / n := by
  rw [Finset.sum_boole]
  exact LPMratio_zero_eq_ecdf t x

theorem claim_C491 (hn : 0 < n) {d : ℝ} (hd : 0 ≤ d) (x : Fin n → ℝ) (t : ℝ) :
    lowerShare 0 t x = cdf (empMeasure x) t ∧ upperShare 0 t x = 1 - cdf (empMeasure x) t ∧
    lowerShare d t x = lowerFormula d t x ∧ upperShare d t x = upperFormula d t x := by
  have h := (claim_C436 hn x t).1
  exact ⟨by simpa [lowerShare, LPMratio] using h,
    by simpa [upperShare, UPMratio] using
        (UPMratio_eq_one_sub_LPMratio hn 0 t x).trans (congrArg (1 - ·) h),
    shares hn hd t x⟩

theorem claim_C522 (t : ℝ) (x : Fin n → ℝ) :
    UPM 0 t x = upperFrequency t x ∧ UPM 1 t x = (∑ i, max (x i - t) 0) / n :=
  ⟨UPMratio_zero_eq_survival t x, by simp [UPM, upperDev_one]⟩

theorem claim_C534 (t : ℝ) (x : Fin n → ℝ) :
    LPMratio 0 t x = lowerFrequency t x := LPMratio_zero_eq_ecdf t x

theorem claim_C537 (hn : 0 < n) {d : ℝ} (hd : 0 < d) (t : ℝ) (x : Fin n → ℝ) :
    (rLPM d t x + rUPM d t x = 0 ↔ ∀ i, x i = t) ∧
    (lowerShare d t x = lowerFormula d t x ∧ upperShare d t x = upperFormula d t x) ∧
    ((∀ i, x i = t) → lowerShare d t x = lowerShare 0 t x ∧ lowerShare d t x = 1 ∧
      upperShare d t x = upperShare 0 t x ∧ upperShare d t x = 0) := by
  refine ⟨total_zero hn hd t x, shares hn hd.le t x, fun ht => ?_⟩
  simp [upperShare, lowerShare, hd.ne', hn, ht, LPM, UPM, lowerDev, upperDev,
    Nat.cast_ne_zero.mpr hn.ne']

theorem claim_C539 (hn : 0 < n) {d : ℝ} (hd : 0 ≤ d) (t : ℝ) (x : Fin n → ℝ) :
    upperShare d t x = upperFormula d t x ∧ upperShare 0 t x = upperFrequency t x ∧
    lowerShare d t x + upperShare d t x = 1 :=
  ⟨(shares hn hd t x).2, (freq_zero t x).2, share_partition hn hd t x⟩

theorem claim_C542 (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    cdfCurve 0 t x = cdf (empMeasure x) t ∧
    (cdfCurve 1 t x = if LPM 1 t x + UPM 1 t x = 0 then 1
      else LPM 1 t x / (LPM 1 t x + UPM 1 t x)) ∧
    (cdfCurve 1 0 ![(0 : ℝ), 1] = 0 ∧ lowerFrequency 0 ![(0 : ℝ), 1] = 1 / 2) :=
  ⟨(claim_C436 hn x t).1, curve_one hn t x, frequency_differs⟩

namespace A4_RatiosCDF

def zeroinStop (ε lo hi b c value : ℝ) : Prop :=
  |(c - b) / 2| ≤ 2 * ε * |b| + ε ^ 2 * (hi - lo) / 2 ∨ value = 0

private theorem zeroin_scale {a : ℝ} (ha : 0 < a) (d : ℕ) (ε lo hi b c value : ℝ) :
    zeroinStop ε (a * lo) (a * hi) (a * b) (a * c) (a ^ d * value) ↔
      zeroinStop ε lo hi b c value := by
  unfold zeroinStop
  rw [← mul_sub, mul_div_assoc, abs_mul, abs_mul, abs_of_pos ha, ← mul_sub]
  have h : 2 * ε * (a * |b|) + ε ^ 2 * (a * (hi - lo)) / 2 =
      a * (2 * ε * |b| + ε ^ 2 * (hi - lo) / 2) := by ring
  rw [h, mul_le_mul_iff_right₀ ha, mul_eq_zero, or_iff_right (pow_ne_zero _ ha.ne')]

end A4_RatiosCDF

theorem claim_C602 (d : ℕ) {a : ℝ} (ha : 0 < a) (p lo hi b c ε : ℝ) (x : Fin n → ℝ) :
    let y := fun i => a * x i
    let balance := fun t (z : Fin n → ℝ) => (1 - p) * LPM d t z - p * UPM d t z
    (∀ t, |LPMratio d (a * t) y - p| ≤ rTol ↔ |LPMratio d t x - p| ≤ rTol) ∧
    (zeroinStop ε (a * lo) (a * hi) (a * b) (a * c) (balance (a * b) y) ↔
      zeroinStop ε lo hi b c (balance b x)) := by
  dsimp only
  constructor
  · intro t
    simpa using congrArg (fun z : ℝ => |z - p| ≤ rTol) (LPMratio_affine d ha 0 t x) |>.to_iff
  · have hL := rLPM_mul_const ha (d : ℝ) b x
    have hU := rUPM_mul_const ha (d : ℝ) b x
    simp only [rLPM_natCast, rUPM_natCast, Real.rpow_natCast] at hL hU
    rw [hL, hU]
    convert zeroin_scale ha d ε lo hi b c ((1 - p) * LPM d b x - p * UPM d b x) using 2
    ring

namespace A4_RatiosCDF

private theorem induced_population (μ : Measure ℝ) [IsProbabilityMeasure μ] {d : ℝ}
    (hd : 0 < d) (hI : Integrable (fun v : ℝ => |v| ^ d) μ) (hx : ∀ t, μ {t}ᶜ ≠ 0) :
    ∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧ ∀ t, cdf ν t = GapClaims.rPMRatio d t μ :=
  induced (GapClaims.rPMRatio_monotone hd hI)
    (fun _ => (GapClaims.continuous_rPMRatio hd hI hx).continuousAt.continuousWithinAt)
    (GapClaims.tendsto_rPMRatio_atBot hd hI) (GapClaims.tendsto_rPMRatio_atTop hd hI)

private theorem induced_sample (hn : 0 < n) {d : ℝ} (hd : 0 < d) (x : Fin n → ℝ) :
    ∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧ ∀ t, cdf ν t = lowerShare d t x := by
  by_cases h : ∃ i j, x i ≠ x j
  · have := isProbabilityMeasure_empMeasure hn x
    obtain ⟨ν, hν, he⟩ := induced_population (empMeasure x) hd
      (GapClaims.integrable_empMeasure hn x _) (empMeasure_compl_singleton_ne_zero h)
    refine ⟨ν, hν, fun t => ?_⟩
    have hnall : ¬ (0 < n ∧ ∀ i, x i = t) := by
      obtain ⟨i, j, hij⟩ := h
      exact fun ht => hij ((ht.2 i).trans (ht.2 j).symm)
    rw [he, lowerShare, ite_eq_right hd.ne', ite_eq_right hnall, GapClaims.rPMRatio]
    change (∫ v, rLowerDev d t v ∂empMeasure x) / _ = _
    simp only [rlpm, rupm, integral_empMeasure, rLPM, rUPM]
  · push Not at h
    let c := x ⟨0, hn⟩
    have he : x = fun _ => c := funext fun i => h i _
    rw [he]
    refine ⟨Measure.dirac c, inferInstance, fun t => ?_⟩
    rw [cdf_eq_real]
    rcases lt_trichotomy t c with hlt | heq | hgt
    · have hne : ¬ ∀ _ : Fin n, c = t := fun e => hlt.ne (e ⟨0, hn⟩).symm
      simp [lowerShare, hd.ne', hne, rLPM_const hn, rLowerDev, not_le.mpr hlt,
        measureReal_def, Measure.dirac_apply' _ measurableSet_Iic]
    · subst t
      simp [lowerShare, hd.ne', hn, measureReal_def, Measure.dirac_apply' _ measurableSet_Iic]
    · have hne : ¬ ∀ _ : Fin n, c = t := fun e => hgt.ne (e ⟨0, hn⟩)
      simp [lowerShare, hd.ne', hne, rLPM_const hn, rUPM_const hn, rLowerDev, rUpperDev,
        hgt.le, not_lt.mpr hgt.le, (Real.rpow_pos_of_pos (sub_pos.mpr hgt) d).ne',
        measureReal_def, Measure.dirac_apply' _ measurableSet_Iic]

end A4_RatiosCDF

theorem claim_C618 (hn : 0 < n) {d : ℝ} (hd : 0 < d) (t : ℝ) (x : Fin n → ℝ) :
    (lowerShare 0 t x = lowerFrequency t x ∧ upperShare 0 t x = upperFrequency t x) ∧
    (lowerShare d t x = lowerFormula d t x ∧ upperShare d t x = upperFormula d t x) ∧
    (∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧ ∀ s, cdf ν s = lowerShare d s x) ∧
    (cdfCurve 1 0 ![(0 : ℝ), 1] = 0 ∧ lowerFrequency 0 ![(0 : ℝ), 1] = 1 / 2) ∧
    ∀ (μ : Measure ℝ) [IsProbabilityMeasure μ], Integrable (fun v : ℝ => |v| ^ d) μ →
      (∀ s, μ {s}ᶜ ≠ 0) → ∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧
        ∀ s, cdf ν s = GapClaims.rPMRatio d s μ :=
  ⟨freq_zero t x, shares hn hd.le t x, induced_sample hn hd x, frequency_differs,
    fun μ _ hI hx => induced_population μ hd hI hx⟩

theorem claim_C610 (hn : 0 < n) (x : Fin n → ℝ) (t : ℝ) {p d : ℝ}
    (hp : p ∈ Set.Icc (0 : ℝ) 1) (hd : 0 < d) :
    LPMratio 0 t x = lowerFrequency t x ∧ lowerShare d t x = lowerFormula d t x ∧
    (∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧ ∀ s, cdf ν s = cdfCurve 1 s x) ∧
    ((∃ i j, x i ≠ x j) → LPMratio 1 (mean x) x = 1 / 2 ∧
      (LPMratio 1 t x = p ↔ (1 - p) * ∑ i, lowerDev 1 t (x i) = p * ∑ i, upperDev 1 t (x i)) ∧
      (p ∈ Set.Ioo (0 : ℝ) 1 → LPMratio 1 (threshold 1 (empMeasure x) p) x = p)) ∧
    (∀ (μ : Measure ℝ) [IsProbabilityMeasure μ], Integrable (fun v : ℝ => v) μ →
      (∀ s, μ {s}ᶜ ≠ 0) → mPMRatio 1 (∫ v, v ∂μ) μ = 1 / 2 ∧
        (p ∈ Set.Ioo (0 : ℝ) 1 → (1 - p) * mLPM 1 (threshold 1 μ p) μ =
          p * mUPM 1 (threshold 1 μ p) μ)) ∧
    (LPMVaR 0 p x = quantile7 p x ∧ UPMVaR 0 p x = quantile7 (1 - p) x ∧
      let r := ((n - 1 : ℕ) : ℝ) * p
      quantile7 p x = sortedAt x ⌊r⌋₊ + (r - ⌊r⌋₊) * (sortedAt x ⌈r⌉₊ - sortedAt x ⌊r⌋₊)) ∧
    (0 < p → threshold 0 (empMeasure x) p = sortedAt x (⌈(n : ℝ) * p⌉₊ - 1)) ∧
    (let y := ![(0 : ℝ), 1, 2]
     let band := Icc (LPMVaR 0 (1 / 10) y) (UPMVaR 0 (1 / 10) y)
     band = Set.Icc (1 / 5) (9 / 5) ∧ (empMeasure y).real band = 1 / 3 ∧
       (empMeasure y).real band ≠ 1 - 2 * (1 / 10)) := by
  obtain ⟨ν, hν, he⟩ := induced_sample hn (d := 1) zero_lt_one x
  refine ⟨LPMratio_zero_eq_ecdf t x, (shares hn hd.le t x).1,
    ⟨ν, hν, fun s => (he s).trans (by
      change lowerShare 1 s x = LPMratio 1 s x
      simpa only [Nat.cast_one] using lowerShare_nat 1 s x)⟩,
    fun hx => ⟨LPMratio_one_mean hx, LPMratio_one_eq_iff_expectile hx t p,
      fun h => (threshold_empMeasure one_ne_zero hx h.1 h.2).2⟩, fun μ _ h hμ =>
      ⟨mPMRatio_one_mean h (hμ _), fun h' => threshold_balance one_ne_zero
        (memLp_one_of_integrable h) hμ h'.1 h'.2⟩, ⟨?_, ?_, quantile7_eq hp x⟩, ?_, ?_⟩
  · change quantile7 (clampP p) x = quantile7 p x
    rw [clampP_of_mem hp]
  · change quantile7 (1 - clampP p) x = quantile7 (1 - p) x
    rw [clampP_of_mem hp]
  · intro hp0
    let := isProbabilityMeasure_empMeasure hn x
    rw [threshold_zero, quantile_empMeasure hn x hp0 hp.2, ecdfInv]
  · dsimp only
    rw [UPMVaR_eq_LPMVaR (by decide), clampP_of_mem (by constructor <;> norm_num),
      (by norm_num : (1 : ℝ) - 1 / 10 = 9 / 10)]
    simp only [LPMVaR, ↓reduceIte]
    rw [LPMVaR0_not_minimal.1, LPMVaR0_type7_example.1]
    have hm : (empMeasure ![(0 : ℝ), 1, 2]).real (Icc (1 / 5) (9 / 5)) = 1 / 3 := by
      norm_num [empMeasure, measureReal_def, Fin.sum_univ_three,
        Measure.dirac_apply' _ measurableSet_Icc, ENNReal.toReal_inv]
    exact ⟨rfl, hm, by rw [hm]; norm_num⟩

end

end NNS.ClaimsA
